use std::{collections::{HashMap, HashSet}, sync::atomic::{AtomicU32, AtomicU64, Ordering::{Acquire, Relaxed, Release, SeqCst}}};

use rand_chacha::rand_core::{RngCore, SeedableRng};
use simple_sds_sbwt::{ops::{BitVec, Rank}, raw_vector::AccessRaw};

use crate::{coloring_interface::ColorSetStorage, iterators::VecIterator};

#[derive(Debug, Copy, Clone)]
struct SetElement {
    set_id: usize,
    color: usize,
}

struct MyTransposedColorSetStream<T : Iterator<Item = SetElement>> {
    element_generator: T,
    buf: Vec<usize>,
    leftover_element: Option<SetElement>,
    current_color: usize,
    n_colors: usize,
}


impl<T : Iterator<Item = SetElement>> crate::iterators::USizeIteratorGenerator for MyTransposedColorSetStream<T> {

    type Iter<'a> = VecIterator<'a> where Self: 'a;

    fn next<'b>(&'b mut self) -> Option<Self::Iter<'b>> {
        self.buf.clear();
        dbg!(&self.current_color, &self.n_colors);
        if self.current_color == self.n_colors {
            return None;
        }

        // If there is a leftover element from the previous round, consider that.
        if let Some(x) = &self.leftover_element {
            if x.color < self.current_color {
                // No set ids in this color
                self.current_color += 1;
                return Some(VecIterator::new(&self.buf));
            } else if x.color == self.current_color {
                self.buf.push(x.set_id);
                self.leftover_element = None;
            } else {
                panic!("Programming error: color set element iterator is not in the right order");
            }
        }

        // Read new elements from the generator
        while let Some(x) = self.element_generator.next() {
            eprintln!("Generator gives {:?}", x);
            if x.color == self.current_color {
                self.buf.push(x.set_id);
            } else {
                self.leftover_element = Some(x);
                break;
            }
        }

        eprintln!("Generator returns {:?}", self.buf);
        self.current_color += 1;

        Some(VecIterator::new(&self.buf))
    }
    
}

/// Takes a generator of SetElement structs with set_id in 0..max_n_sets and element in 0..max_n_elements.
/// The element generators must generate the elements in increasing order of element: first all set ids
/// with element 0, then all set ids with element 1, and so on. There must not be duplicate elements
/// in the same set.
/// Returns the CSS and a vector of length n_sets mapping original set ids to new set ids.
fn construct<CSS: ColorSetStorage>(
    element_generator: impl Iterator<Item = SetElement>, 
    element_generator_again: impl Iterator<Item = SetElement>, 
    n_sets: usize, 
    n_colors: usize,
    random_seed: usize)
    -> (CSS, Vec<usize>) {

    // TODO: all the atomic operation here are with SeqCst ordering because I'm not sure
    // if the stores here are guaranteed to happen before loads. There is probably a better
    // way to do this with Acquire/Release ordering, or with fences.

    // Assign a 128-bit fingerprint for each possible element id. 128-bit integers can not be,
    // updated atomically, so instead we use a pair of u64 values which can be updated atomically.
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(random_seed as u64);
    let element_fingerprints: Vec<(u64,u64)> = (0..n_colors).map(|_i| (rng.next_u64(), rng.next_u64())).collect();

    // 128-bit fingerprints for sets of elements. Again we split each u128 into
    // two u64s.
    let mut set_fingerprints = Vec::<(AtomicU64, AtomicU64)>::new();
    set_fingerprints.resize_with(n_sets, || (AtomicU64::new(0), AtomicU64::new(0)));
    let mut set_sizes = Vec::<AtomicU64>::new(); // TODO: could be U32?
    set_sizes.resize_with(n_sets, || AtomicU64::new(0));

    for new in element_generator {
        let (fp1, fp2) = element_fingerprints[new.color];

        set_fingerprints[new.set_id].0.fetch_xor(fp1, SeqCst);
        set_fingerprints[new.set_id].1.fetch_xor(fp2, SeqCst);
        set_sizes[new.set_id].fetch_add(1, SeqCst);
    } 

    // Make set fingeprints not atomic
    let set_fingerprints: Vec<(u64, u64)> = set_fingerprints.into_iter().map(
        |pair| (pair.0.load(SeqCst), pair.1.load(SeqCst))
    ).collect();

    // Make sizes not atomic
    let set_sizes: Vec<usize> = set_sizes.into_iter().map(
        |sz| sz.load(SeqCst) as usize
    ).collect();

    // Mark the lowest set id where each distinct fingerprint occurs 
    let mut distinct_fingerprints = HashMap::<(u64,u64), usize>::new(); // Maps fingerprint to new set id
    let mut marked_sets = simple_sds_sbwt::raw_vector::RawVector::with_len(n_sets, false);
    let mut marked_set_sizes = Vec::<usize>::new();
    let mut n_distinct_set_found = 0_usize;
    for set_id in 0..n_sets {
        let fp1 = set_fingerprints[set_id].0;
        let fp2 = set_fingerprints[set_id].1;
        let fp = (fp1, fp2);

        if let std::collections::hash_map::Entry::Vacant(e) = distinct_fingerprints.entry(fp) {
            e.insert(n_distinct_set_found);
            n_distinct_set_found += 1;
            marked_sets.set_bit(set_id, true);
            marked_set_sizes.push(set_sizes[set_id] as usize);
        }
    }

    // Free memory
    drop(set_sizes);
    drop(element_fingerprints);

    // Build original set id -> new set id vector
    let old_id_to_new_id: Vec<usize> = set_fingerprints.iter().enumerate().map(
        |(_old_id, fingerprint)| distinct_fingerprints[fingerprint])
        .collect();

    // Free memory
    drop(set_fingerprints);

    // Build ran on marked sets
    let mut marked_sets = simple_sds_sbwt::bit_vector::BitVector::from(marked_sets);
    marked_sets.enable_rank();

    // Now we build the ColorSetStorage from the transposed constructor, that is,
    // we need a ColorSetStream that is like an iterator of iterators, where each
    // inner iterator gives is the set ids of all sets that have a given color.
    // Here we make use of the assumption that the element generator generates the
    // elements with increasing order of color id.
    let my_stream = MyTransposedColorSetStream {
        element_generator: element_generator_again.filter(|new| { 
            // Only include sets whose is id sampled
            marked_sets.get(new.set_id)
        }).map(|new| {
            let rank_in_sampled_sets = marked_sets.rank(new.set_id);
            SetElement { set_id: rank_in_sampled_sets, color: new.color }
        }),
        buf: vec![],
        leftover_element: None,
        current_color: 0,
        n_colors
    };

    (*CSS::new_from_transpose(my_stream, n_colors, &marked_set_sizes), old_id_to_new_id)

}

#[cfg(test)]
mod tests{
    use crate::{coloring_interface::ColorSetView, sparse_dense_storage::SparseDenseStorage};

    use super::*;

    #[test]
    fn set_of_sets_test() {
        // Define some sets of sets
        let sets = [
            vec![0, 1, 2],
            vec![2, 3],
            vec![0, 1, 2], // duplicate of first set
            vec![4],
            vec![], // Make sure to test the empty set because it is a special case
            vec![3, 2],    // duplicate of second set
        ];

        let transposed_sets: Vec<Vec<usize>> = vec![
            vec![0, 2], // Color 0
            vec![0, 2], // Color 1
            vec![0, 1, 2, 5], // Color 2
            vec![1, 5], // Color 3
            vec![3], // Color 4
        ];

        // Make an element generator
        let mut elements = vec![];
        transposed_sets.iter().enumerate().for_each(|(color, set_ids)| {
            for set_id in set_ids {
                elements.push(SetElement{set_id: *set_id, color});
            }
        });

        dbg!(&elements);
        
        let distinct_sets = construct::<SparseDenseStorage>(
            elements.iter().copied(),
            elements.iter().copied(),
            sets.len(),
            5,
            123123
        );

        let mut correct_answers = vec![vec![0,1,2], vec![2,3], vec![4], vec![]];
        let mut our_answers: Vec<Vec<usize>> = 
            (0..distinct_sets.n_sets())
            .map(|i| distinct_sets.get_set_view(i).iter().collect::<Vec::<usize>>())
            .collect();

        correct_answers.sort();
        our_answers.sort();
        for i in 0..distinct_sets.n_sets() {
            eprintln!("{:?} {:?}", our_answers[i], correct_answers[i]);
        }
        assert_eq!(correct_answers, our_answers);

    }
}