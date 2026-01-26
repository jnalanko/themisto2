use std::{collections::HashMap, ops::Range, path::Path, sync::{Arc, atomic::{AtomicU64, Ordering::{Acquire, Release}}}};

use rand_chacha::rand_core::{RngCore, SeedableRng};
use rayon::slice::ParallelSliceMut;
use simple_sds_sbwt::ops::{BitVec, Rank, Select};
use crate::{coloring_interface::ColorSetStorage, int_vec::CompactIntVec, util};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct SetElement {
    pub set_id: usize,
    pub color: usize,
}

pub trait ParallelElementGenerator {
    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize);

    // Takes a bit vector with rank support that indicates which sets should be passed
    // to the callback of run. The set ids will be remapped so that the first set that
    // passes the filter gets id 0, the second that passes the filter gets id 1, and so on.
    // It's Arc because if we have multiple generators, we need them to be able to share
    // the filter.
    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>);

    fn rewind(&mut self); // Rewind to the start for another run
}

/// Takes a generator of SetElement structs with set_id in 0..max_n_sets and element in 0..max_n_elements.
/// There must not be duplicate elements in the same set! The callback only has to return ids for
/// sets that are marked in the key_kmer_marks bitvector.
/// Returns three things:
/// * A bit vector marking a subset of the key-kmers such that every marked k-mer has a distinct color
///   set.
/// * The sizes of the color sets of the marked k-mers, in colex order
/// * A vector of length key_kmer_marks.count_ones() that gives the color set id for each key k-mer.
///   The color set id is the rank of the 1-bit of the representative k-mer of the color set in the returned marks.
pub fn find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(
    mut gen: impl ParallelElementGenerator,
    key_kmer_marks: bitvec::vec::BitVec,
    n_colors: usize, n_threads: usize, random_seed: usize)
    -> (bitvec::vec::BitVec, Vec<usize>, CompactIntVec) {

    // Build rank support for key k-mer marks
    log::info!("Building rank support for key k-mer marks");
    let key_kmer_marks = crate::util::bitvec_to_simple_sds_raw_bitvec(key_kmer_marks);
    let mut key_kmer_marks = simple_sds_sbwt::bit_vector::BitVector::from(key_kmer_marks);
    key_kmer_marks.enable_rank();
    let n_key_kmers = key_kmer_marks.count_ones();

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

    drop(element_fingerprints); // Free memory

    log::info!("Sorting by fingerprint");
    // Make set fingeprints not atomic and add colex positions as the third element and the set
    // size as the fourth element
    let mut set_quadruples: Vec<(u64, u64, usize, usize)> = set_fingerprints.into_iter().map(
        |pair| (pair.0.load(Acquire), pair.1.load(Acquire), 0, 0) // Colex and set size will be filled in next
        // The Acquire here means that the Release-writes above must be visible before these loads (I think).
    ).collect();
    for (fp_idx, colex) in key_kmer_marks.one_iter() {
        set_quadruples[fp_idx].2 = colex;
    }
    // Make sizes not atomic and add to the quadruples
    for (idx, sz) in set_sizes.into_iter().enumerate() {
        set_quadruples[idx].3 = sz.load(Acquire) as usize;
    }
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| set_quadruples.par_sort_unstable());

    log::info!("Assigning representative k-mers");
    let mut sufficient_kmer_marks = bitvec::bitvec![0; key_kmer_marks.len()];
    util::for_each_run_with_key_mut(&mut set_quadruples, |x| (x.0, x.1), |run| {
        let min_colex = run.iter().map(|x| x.2).min().unwrap(); // Min colex in the run. Unwrap is okay because there can not be empty runs
        for x in run {
            x.0 = min_colex as u64; // Store the min colex over the first fingerprint
        }
        assert!(!sufficient_kmer_marks[min_colex]);
        sufficient_kmer_marks.set(min_colex, true);
    });

    log::info!("Sorting by representative id");
    thread_pool.install(|| set_quadruples.par_sort_unstable()); // Sort by min colex
    let mut key_kmer_idx_to_new_id: Vec<usize> = vec![0; n_key_kmers]; 
    let mut set_sizes: Vec<usize> = vec![0; sufficient_kmer_marks.count_ones()];
    let mut set_id = 0_usize;
    log::info!("Collecting set ids and sizes");
    util::for_each_run_with_key(&set_quadruples, |x| x.0, |run| { // Run with the same min colex
        let mut set_size = 0;
        for (_, _, key_kmer_colex, size) in &set_quadruples[run.clone()] {
            let key_kmer_idx = key_kmer_marks.rank(*key_kmer_colex);
            key_kmer_idx_to_new_id[key_kmer_idx] = set_id;
            set_size = *size; // The size should be the same for all lets since this is a run
        }
        set_sizes[set_id] = set_size;
        set_id += 1;
    });

    let n_sets = set_id;
    log::info!("{} distinct color sets found", n_sets);
    log::info!("Average color set size: {:.2}", (set_sizes.iter().sum::<usize>() as f64) / (n_sets as f64));

    (sufficient_kmer_marks, set_sizes, CompactIntVec::from_vec(key_kmer_idx_to_new_id))

}

// The generator must not provide duplicate elements! Otherwise the final data structure will be corrupted because
// the sampled_set_sizes is not correct.
pub fn build_color_set_storage<CSS: ColorSetStorage + Send>(n_colors: usize, colex_sample_marks: bitvec::vec::BitVec, sampled_set_sizes: Vec<usize>, mut gen: impl ParallelElementGenerator, n_threads: usize) -> CSS {
    let mut colex_sample_marks = crate::util::bitvec_to_simple_sds_bitvec(colex_sample_marks);
    colex_sample_marks.enable_rank();

    gen.set_filter(Arc::new(colex_sample_marks));

    *CSS::new_parallel(gen, n_colors, &sampled_set_sizes, n_threads)
}

// The generator must not provide duplicate elements! Otherwise the final data structure will be corrupted because
// the sampled_set_sizes is not correct.
pub fn build_color_set_storage_to_disk<CSS: ColorSetStorage + Send>(n_colors: usize, colex_sample_marks: bitvec::vec::BitVec, sampled_set_sizes: Vec<usize>, mut element_gens: Vec<(impl crate::set_of_sets_construction::ParallelElementGenerator, Range<usize>)>, out_prefix: &Path, n_threads: usize) {
    let mut colex_sample_marks = crate::util::bitvec_to_simple_sds_bitvec(colex_sample_marks);
    colex_sample_marks.enable_rank();

    let colex_sample_marks = Arc::new(colex_sample_marks);

    for (gen, _) in element_gens.iter_mut() {
        gen.set_filter(colex_sample_marks.clone());
    }

    CSS::new_parallel_to_disk(element_gens, sampled_set_sizes, out_prefix, n_threads);
}

pub fn compute_set_sizes_assuming_no_duplicate_elements(element_gen: &mut impl ParallelElementGenerator, n_set_ids: usize, n_threads: usize) -> Vec<usize> {
    let mut sizes = Vec::<AtomicU64>::new(); // TODO: could be U32?
    sizes.resize_with(n_set_ids, || AtomicU64::new(0));

    let callback = |e: SetElement| {
        sizes[e.set_id].fetch_add(1, Release);
    };

    element_gen.run(callback, n_threads);

    let sizes = sizes.into_iter().map(
        |x| x.load(Acquire) as usize
    ).collect();
    sizes
}

#[cfg(test)]
mod tests{
    use simple_sds_sbwt::ops::BitVec;

    use crate::{coloring_interface::ColorSetView, sparse_dense_storage::SparseDenseStorage};
    use bitvec::prelude::*;

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

        fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
            self.filter = Some((*filter).clone());
        }

        fn rewind(&mut self) {
            // Nothing needs to done, calling run() again already works
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
        
        let (new_marks, set_sizes, key_kmer_to_set_id) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(
            VecVecGenerator::new(sets.clone()),
            bitvec::bitvec![1,1,0,1,1,1], // Mark one of the duplicates as non-key
            sets.len(),
            3,
            123123
        );
        let distinct_sets = build_color_set_storage::<SparseDenseStorage>(5, new_marks, set_sizes, VecVecGenerator::new(sets.clone()), 3);

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
        assert_eq!(key_kmer_to_set_id.to_vec(), vec![0,1,2,3,1]);

    }
}