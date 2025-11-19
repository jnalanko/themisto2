use std::{collections::HashSet, sync::atomic::{AtomicU32, AtomicU64, Ordering::Relaxed}};

use rand_chacha::rand_core::{RngCore, SeedableRng};

struct SetElement {
    set_id: usize,
    element: usize,
}

/// Takes a generator of SetElement structs with set_id in 0..max_n_sets and element in 0..max_n_elements
fn construct(
    element_generator: impl Iterator<Item = SetElement>, 
    element_generator_again: impl Iterator<Item = SetElement>, 
    max_n_sets: usize, 
    max_n_elements: usize,
    random_seed: usize,)
    -> Vec<Vec<usize>> {


    // Assign a 128-bit fingerprint for each possible element id. 128-bit integers can not be,
    // we instead use a pair of u64-bit values which can be updated atomically.
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(random_seed as u64);
    let element_fingerprints: Vec<(u64,u64)> = (0..max_n_elements).map(|_i| (rng.next_u64(), rng.next_u64())).collect();

    // 128-bit fingerprints for sets of elements. Again we split each u128 into
    // two u64s.
    let mut set_fingerprints = Vec::<(AtomicU64, AtomicU64)>::new();
    set_fingerprints.resize_with(max_n_sets, || (AtomicU64::new(0), AtomicU64::new(0)));

    for new in element_generator {
        let (fp1, fp2) = element_fingerprints[new.element];

        set_fingerprints[new.set_id].0.fetch_xor(fp1, Relaxed);
        set_fingerprints[new.set_id].1.fetch_xor(fp2, Relaxed);
    } 

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

    // Iterate sets again and store the marked sets
    let mut distinct_sets: Vec<Vec<usize>> = vec![vec![]; distinct_fingerprints.len()];
    for new in element_generator_again {
        if marked_sets[new.set_id] {
            let distinct_id = marked_sets[..new.set_id].count_ones(); // TODO: rank query
            distinct_sets[distinct_id].push(new.element); 
        }
    }

    distinct_sets
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
        assert!(false);
    }
}