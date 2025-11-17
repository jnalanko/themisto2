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

    // Biggest prime under 2^90:
    let p: u128 = 1237940039285380274899124191;

    // Draw the base of the hash from 0..p. We can only draw 64 bits at a time,
    // so we'll combine two 64-bit numbers and take mod p of that.
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(random_seed as u64);
    let b = (rng.next_u64() as u128 + ((rng.next_u64() as u128) << 64)) % p;

    // Build fingerprints. Initialize a conceptual array of max_n_set
    // 128-bit integers, but since we can't update 128-bit integers
    // atomically, instead we simulate the array with twice the amount
    // of 64-bit integers.
    let mut A = Vec::<AtomicU64>::new();
    A.resize_with(max_n_sets*2, || AtomicU64::new(0));

    // Least significant word carries. Up to 2^32 carries supported, which means
    // sets of size up to 2^32.
    let mut lsw_carries = Vec::<AtomicU32>::new();
    lsw_carries.resize_with(max_n_sets, || AtomicU32::new(0));

    for new in element_generator {
        let c = new.element as u128;
        let to_add = mod_pow(b, c, p);
        let to_add_lsw = to_add as u64; // Least significant word
        let to_add_msw = (to_add >> 64) as u64; // Least significant word

        let lsw_before = A[2*new.set_id].fetch_add(to_add_lsw, Relaxed);
        let _msw_before = A[2*new.set_id + 1].fetch_add(to_add_msw, Relaxed);

        // Carry happened if lsw_before + to_add_lsw >= 2^64
        // Which is the same as lsw_before >= 2^64 - to_add_lsw.
        // This is same as:
        if lsw_before >= to_add_lsw.wrapping_neg() {
            lsw_carries[new.set_id].fetch_add(1, Relaxed); // Record the carry
        }
    } 

    // Mark the lowest set id where each distinct fingerprint occurs 
    let mut distinct_fingerprints = HashSet::<u128>::new();
    let mut marked_sets = bitvec::bitvec![0; max_n_sets];
    for set_id in 0..max_n_sets {
        let lsw = A[2*set_id].load(Relaxed) as u128;
        let msw = A[2*set_id + 1].load(Relaxed) as u128;
        let carries = lsw_carries[set_id].load(Relaxed) as u128;
        let mut total: u128 = lsw + (msw << 64);
        total %= p;
        total += carries << 64;
        total %= p;

    }
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

fn mod_pow(mut base: u128, mut exp: u128, modulo: u128) -> u128 {
    let mut result = 1_u128;
    base %= modulo;

    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, modulo);
        }
        base = mul_mod(base, base, modulo);
        exp >>= 1;
    }

    result
}

/// Compute (a * b) % m over u128.
/// Uses the ancient egyptian multiplication algorithm:
/// https://en.wikipedia.org/wiki/Ancient_Egyptian_multiplication 
/// IMPORTANT: assumes a < m and b < m (and of course m > 0)
pub fn mul_mod(mut a: u128, mut b: u128, m: u128) -> u128 {
    debug_assert!(m > 0);
    debug_assert!(a < m);
    debug_assert!(b < m);

    let mut res = 0u128;

    while b > 0 {
        if b & 1 == 1 {
            // res = (res + a) modulo m
            res = add_mod(res, a, m)
        }

        // a = 2a modulo b
        a = add_mod(a, a, m); 

        // b = floor(b/2)
        b >>= 1;
    }

    res
}

fn add_mod(r: u128, a: u128, m: u128) -> u128 {
    debug_assert!(m > 0);
    debug_assert!(r < m);
    debug_assert!(a < m);

    let t = m - a;
    if r >= t {
        // r + a >= m, so result is r + a - m == r - (m - a)
        r - t
    } else {
        // r + a < m, safe to add directly
        r + a
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn modular_arithmetic() {
        // Prime close to 2^128
        let p: u128 = 340282366920938463463374607431768211297;
        assert_eq!(mod_pow(5, 128, p), 124204064441846203354619122274609400273);
    }

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
            340282366920938463463374607431768211297,
             12312312312312312312312312312312312312
        );

        dbg!("Distinct sets: {:?}", distinct_sets);
        assert!(false);
    }
}