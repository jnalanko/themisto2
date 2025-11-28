use std::{collections::HashMap, sync::atomic::{AtomicU64, Ordering::{Acquire, Release}}};

use bitvec::order::Lsb0;
use rand_chacha::rand_core::{RngCore, SeedableRng};
use sbwt::{LcsArray, SbwtIndex, SubsetMatrix};
use simple_sds_sbwt::{ops::{BitVec, Rank, Select}, raw_vector::AccessRaw};
use crate::{colex_colored_kmers::mark_key_kmers, coloring_interface::ColorSetStorage, iterators::VecIterator};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct SetElement {
    pub set_id: usize,
    pub color: usize,
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
        // Todo: can we not store all in a buffer and return and iterator instead?
        while let Some(x) = self.element_generator.next() {
            if x.color == self.current_color {
                self.buf.push(x.set_id);
            } else {
                self.leftover_element = Some(x);
                break;
            }
        }

        self.current_color += 1;

        Some(VecIterator::new(&self.buf))
    }
    
}

pub trait ParallelElementGenerator {
    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize);

    // Takes a bit vector with rank support that indicates which sets should be passed
    // to the callback of run. The set ids will be remapped so that the first set that
    // passes the filter gets id 0, the second that passes the filter gets id 1, and so on.
    fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector);
}

/// Takes a generator of SetElement structs with set_id in 0..max_n_sets and element in 0..max_n_elements.
/// There must not be duplicate elements in the same set!
/// Returns the CSS and a vector of length key_kmer_marks.count_ones() that gives the color set id for each key k-mer.
pub fn construct_from_generators_that_do_not_give_duplicates<CSS: ColorSetStorage + Send>(
    mut gen: impl ParallelElementGenerator,
    mut gen_again: impl ParallelElementGenerator,
    key_kmer_marks: bitvec::vec::BitVec,
    n_sets: usize, n_colors: usize, n_threads: usize, random_seed: usize)
    -> (CSS, Vec<usize>) {

    // Build rank support for key k-mer marks
    log::info!("Building rank support for key k-mer marks");
    let key_kmer_marks = crate::util::bitvec_to_simple_sds_raw_bitvec(key_kmer_marks);
    /*let mut new_key_kmer_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(key_kmer_marks.len(), false);
    for b in key_kmer_marks.iter_ones(){
        new_key_kmer_marks.set_bit(b, true);
    }
    drop(key_kmer_marks);
    */
    let mut key_kmer_marks = simple_sds_sbwt::bit_vector::BitVector::from(key_kmer_marks);
    key_kmer_marks.enable_rank();

    log::info!("Building color set fingerprints");
    // Assign a 128-bit fingerprint for each possible element id. 128-bit integers can not be
    // updated atomically, so instead we use a pair of u64 values which can each be updated atomically.
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(random_seed as u64);
    let element_fingerprints: Vec<(u64,u64)> = (0..n_colors).map(|_i| (rng.next_u64(), rng.next_u64())).collect();

    // 128-bit fingerprints for the color set of each key k-mer. Again we split each u128 into
    // two u64s.
    let mut set_fingerprints = Vec::<(AtomicU64, AtomicU64)>::new();
    set_fingerprints.resize_with(key_kmer_marks.count_ones(), || (AtomicU64::new(0), AtomicU64::new(0)));
    let mut set_sizes = Vec::<AtomicU64>::new(); // TODO: could be U32?
    set_sizes.resize_with(key_kmer_marks.count_ones(), || AtomicU64::new(0));
    assert!(set_fingerprints.len() == key_kmer_marks.rank(key_kmer_marks.len()));

    let callback = |e: SetElement| {
        if key_kmer_marks.get(e.set_id) {
            let (fp1, fp2) = element_fingerprints[e.color];
            set_fingerprints[key_kmer_marks.rank(e.set_id)].0.fetch_xor(fp1, Release);
            set_fingerprints[key_kmer_marks.rank(e.set_id)].1.fetch_xor(fp2, Release);
            set_sizes[key_kmer_marks.rank(e.set_id)].fetch_add(1, Release);
        }
    };

    gen.run(callback, n_threads);

    // Make set fingeprints not atomic
    let set_fingerprints: Vec<(u64, u64)> = set_fingerprints.into_iter().map(
        |pair| (pair.0.load(Acquire), pair.1.load(Acquire))
        // The Acquire here means that the Release-writes above must be visible before these loads (I think).
    ).collect();

    // Make sizes not atomic
    let set_sizes: Vec<usize> = set_sizes.into_iter().map(
        |sz| sz.load(Acquire) as usize
    ).collect();

    // Mark the colex-lowest key k-mer where each distinct fingerprint occurs 
    let mut distinct_fingerprints = HashMap::<(u64,u64), usize>::new(); // Maps fingerprint to new set id
    let mut sparsified_key_kmer_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(key_kmer_marks.len(), false);
    let mut sparsified_marked_set_sizes = Vec::<usize>::new();
    let mut n_distinct_sets_found = 0_usize;
    let mut total_set_size = 0_usize;
    for (key_kmer_idx, key_kmer_colex) in key_kmer_marks.one_iter() {
        let fp1 = set_fingerprints[key_kmer_idx].0;
        let fp2 = set_fingerprints[key_kmer_idx].1;
        let fp = (fp1, fp2);

        if let std::collections::hash_map::Entry::Vacant(e) = distinct_fingerprints.entry(fp) {
            e.insert(n_distinct_sets_found);
            n_distinct_sets_found += 1;
            sparsified_key_kmer_marks.set_bit(key_kmer_colex, true);
            sparsified_marked_set_sizes.push(set_sizes[key_kmer_idx]);
            total_set_size += set_sizes[key_kmer_idx];
        }
    }
    sparsified_marked_set_sizes.shrink_to_fit();

    log::info!("{} distinct color sets found", n_distinct_sets_found);
    log::info!("Average color set size: {:.2}", (total_set_size as f64)/(n_distinct_sets_found as f64));

    // Free memory
    drop(set_sizes);
    drop(element_fingerprints);

    // Build original set id -> new set id vector
    let key_kmer_idx_to_new_id: Vec<usize> = set_fingerprints.iter().map(
        |fingerprint| distinct_fingerprints[fingerprint])
        .collect();

    // Free memory
    drop(set_fingerprints);

    // Build rank on marked sets
    let mut sparsified_key_kmer_marks = simple_sds_sbwt::bit_vector::BitVector::from(sparsified_key_kmer_marks);
    sparsified_key_kmer_marks.enable_rank();

    // Filter the second element iterator
    gen_again.set_filter(sparsified_key_kmer_marks);

    (*CSS::new_parallel(gen_again, n_colors, &sparsified_marked_set_sizes, n_threads), key_kmer_idx_to_new_id)
}

#[cfg(test)]
mod tests{
    use simple_sds_sbwt::ops::BitVec;

    use crate::{coloring_interface::ColorSetView, sparse_dense_storage::SparseDenseStorage};

    use super::*;

    struct VecVecGenerator {
        vv: Vec<Vec<usize>>,
        filter: Option<simple_sds_sbwt::bit_vector::BitVector>,
    }

    impl VecVecGenerator {
        pub fn new(vv: Vec<Vec<usize>>) -> VecVecGenerator {
            VecVecGenerator{vv, filter: None}
        }
    }

    impl ParallelElementGenerator for VecVecGenerator {
        fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, _n_threads: usize) {
            for (set_id, set) in self.vv.iter().enumerate() {
                let mut new_id = set_id;
                if let Some(filter) = &self.filter {
                    if !filter.get(set_id) {
                        continue; // This is filtered away
                    } else {
                        // Assign new id
                        new_id = filter.rank(set_id);
                    }
                }
                for &color in set.iter() {
                    callback(SetElement { set_id: new_id, color});
                }
            }
        }

        fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector) {
            self.filter = Some(filter);
        }
    }

    #[test]
    fn set_of_sets_test() {
        // Define some sets of sets
        let sets = vec![
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
        
        let (distinct_sets, old_id_to_new_id) = construct_from_generators_that_do_not_give_duplicates::<SparseDenseStorage>(
            VecVecGenerator::new(sets.clone()),
            VecVecGenerator::new(sets.clone()),
            bitvec::bitvec![1,1,0,1,1,1], // Mark one of the duplicates as non-key
            sets.len(),
            5,
            3,
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
        assert_eq!(old_id_to_new_id, vec![0,1,2,3,1]);

    }
}