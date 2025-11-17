use std::collections::HashSet;

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
    fingerprint_modulus: u64,
    fingerprint_base: u64)
    -> Vec<Vec<usize>> {

    let p = fingerprint_modulus; // TODO: assert that this is prime
    let b = fingerprint_base;

    // Build fingerprints
    let mut A: Vec<u64> = vec![0; max_n_sets]; // Fingerprints
    for new in element_generator {
        let c = new.element as u64;
        A[new.set_id] = (A[new.set_id] + mod_pow(b, c, p)) % p;
    } 

    // Mark the lowest set id where each distinct fingerprint occurs 
    let mut distinct_fingerprints = HashSet::<u64>::new();
    let mut marked_sets = bitvec::bitvec![0; max_n_sets];
    for (set_id, fingerprint) in A.into_iter().enumerate() {
        if !distinct_fingerprints.contains(&fingerprint) {
            distinct_fingerprints.insert(fingerprint);
            marked_sets.set(set_id, true);
        }
    }

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

fn mod_pow(mut base: u64, mut exp: u64, modulo: u64) -> u64 {
    let mut result = 1u64;
    base %= modulo;

    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % modulo;
        }
        base = base * base % modulo;
        exp >>= 1;
    }

    result
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
            1_000_000_007,
              123_456_789
        );

        dbg!("Distinct sets: {:?}", distinct_sets);
        assert!(false);
    }
}