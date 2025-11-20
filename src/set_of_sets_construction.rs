use std::{collections::HashSet, sync::atomic::{AtomicU32, AtomicU64, Ordering::Relaxed}};

use rand_chacha::rand_core::{RngCore, SeedableRng};

use crate::coloring_interface::{ColorSetStorage, ColorSetStream};

struct SetElement {
    set_id: usize,
    element: usize,
}

struct MyTransposedColorSetStream<T : Iterator<Item = SetElement>> {
    element_generator: T,
    buf: Vec<usize>,
    leftover_element: Option<SetElement>,
    current_color: usize,
}


impl<T : Iterator<Item = SetElement>> ColorSetStream for MyTransposedColorSetStream<T> {
    fn next(&mut self) -> Option<&[usize]> { // TODO: &[usize] is bad here, should be Iter<Item = usize>
        self.buf.clear();

        // If there is a leftover element from the previous round, consider that.
        if let Some(x) = &self.leftover_element {
            if x.element < self.current_color {
                // No set ids in this color
                self.current_color += 1;
                return Some(&self.buf);
            } else if x.element == self.current_color {
                self.buf.push(x.set_id);
                self.leftover_element = None;
            } else {
                panic!("Programming error: color set element iterator is not in the right order");
            }
        }

        // Read new elements from the generator
        while let Some(x) = self.element_generator.next() {
            if x.element == self.current_color {
                self.buf.push(x.set_id);
            } else {
                self.leftover_element = Some(x);
                self.current_color += 1;
                break;
            }
        }

        return Some(&self.buf);
    }
}

/// Takes a generator of SetElement structs with set_id in 0..max_n_sets and element in 0..max_n_elements.
/// The element generators must generate the elements in increasing order of element: first all set ids
/// with element 0, then all set ids with element 1, and so on.
fn construct<CSS: ColorSetStorage>(
    element_generator: impl Iterator<Item = SetElement>, 
    element_generator_again: impl Iterator<Item = SetElement>, 
    max_n_sets: usize, 
    max_n_elements: usize,
    random_seed: usize)
    -> CSS {


    // Assign a 128-bit fingerprint for each possible element id. 128-bit integers can not be,
    // updated atomically, so instead we use a pair of u64 values which can be updated atomically.
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(random_seed as u64);
    let element_fingerprints: Vec<(u64,u64)> = (0..max_n_elements).map(|_i| (rng.next_u64(), rng.next_u64())).collect();

    // 128-bit fingerprints for sets of elements. Again we split each u128 into
    // two u64s.
    let mut set_fingerprints = Vec::<(AtomicU64, AtomicU64)>::new();
    set_fingerprints.resize_with(max_n_sets, || (AtomicU64::new(0), AtomicU64::new(0)));
    let set_sizes = Vec::<AtomicU64>::new(); // TODO: could be U32?

    for new in element_generator {
        let (fp1, fp2) = element_fingerprints[new.element];

        set_fingerprints[new.set_id].0.fetch_xor(fp1, Relaxed);
        set_fingerprints[new.set_id].1.fetch_xor(fp2, Relaxed);
        set_sizes[new.set_id].fetch_add(1, Relaxed);
    } 

    // From atomic ints to regular ints
    let set_sizes: Vec<usize> =  set_sizes.into_iter().map(|x| x.into_inner() as usize).collect();

    // Mark the lowest set id where each distinct fingerprint occurs 
    let mut distinct_fingerprints = HashSet::<(u64,u64)>::new();
    let mut marked_sets = bitvec::bitvec![0; max_n_sets];
    for set_id in 0..max_n_sets {
        let fp1 = set_fingerprints[set_id].0.load(Relaxed);
        let fp2 = set_fingerprints[set_id].1.load(Relaxed);
        let fp = (fp1, fp2);

        if !distinct_fingerprints.contains(&fp) {
            distinct_fingerprints.insert(fp);
            marked_sets.set(set_id, true);
        }
    }

    // Free memory
    drop(set_fingerprints);
    drop(element_fingerprints);

    // Now we build the ColorSetStorage from the transposed constructor, that is,
    // we need a ColorSetStream that is like an iterator of iterators, where each
    // inner iterator gives is the set ids of all sets that have a given color.
    // Here we make use of the assumption that the element generator generates the
    // elements with increasing order of color id.
    let my_stream = MyTransposedColorSetStream {
        element_generator: element_generator_again.filter(|new| { marked_sets[new.set_id] }),
        buf: vec![],
        leftover_element: None,
        current_color: 0,
    };

    *CSS::new_from_transpose(my_stream, max_n_elements, &set_sizes)

}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn set_of_sets_test() {
        // Define some sets of sets
        let sets = [
            vec![0, 1, 2],
            vec![2, 3],
            vec![0, 1, 2], // duplicate of first set
            vec![4],
            vec![3, 2],    // duplicate of second set
        ];

        // Make an element generator
        let element_generator = sets.iter().enumerate().flat_map(|(set_id, elements)| {
            elements.iter().map(move |&element| SetElement { set_id, element })
        });

        let distinct_sets = construct(
            element_generator.clone(),
            element_generator,
            sets.len(),
            5,
            123123
        );

        dbg!("Distinct sets: {:?}", distinct_sets);
        assert!(false); // TODO
    }
}