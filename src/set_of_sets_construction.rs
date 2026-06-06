use std::{ops::Range, path::Path, sync::{Arc, atomic::{AtomicU64, Ordering::{Acquire, Release}}}};

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


// A generator the wraps a generator and adds an offset to all color ids
#[allow(dead_code)] // Could be useful in the future
pub struct GenWithColorIdOffset<ParallelElementGenerator> {
    pub inner: ParallelElementGenerator,
    pub offset: usize,
}

impl<T: ParallelElementGenerator> ParallelElementGenerator for GenWithColorIdOffset<T> {
    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        self.inner.run(|x| {
            callback(SetElement { set_id: x.set_id, color: x.color + self.offset});
        }, n_threads);
    }

    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.inner.set_filter(filter);
    }

    fn rewind(&mut self) {
        self.inner.rewind();
    }
}

fn log_memory_usage() {
    if let Some(stats) = memory_stats::memory_stats() {
        log::info!("Current memory usage: {} ", human_bytes::human_bytes(stats.physical_mem as f64));
    }
}

pub static SET_OF_SETS_CONSTRUCTION_MAX_SBWT_LEN: usize = 1_usize << 40; // We are bit-packing to 40 bits
/// Takes a generator of SetElement structs with set_id in 0..max_n_sets and element in 0..max_n_elements.
/// There must not be duplicate elements in the same set! The callback only has to return ids for
/// sets that are marked in the key_kmer_marks bitvector.
/// Returns four things:
/// * A bit vector marking a subset of the key-kmers such that every marked k-mer has a distinct color
///   set.
/// * The sizes of the color sets of the marked k-mers, in colex order
/// * A vector of length key_kmer_marks.count_ones() that gives the color set id for each key k-mer.
///   The color set id is the rank of the 1-bit of the representative k-mer of the color set in the returned marks.
/// * The provided key_kmer_marks as a simple_sds_sbwt bitvector, now with rank support
pub fn find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(
    mut gen: impl ParallelElementGenerator,
    key_kmer_marks: bitvec::vec::BitVec,
    n_colors: u32, n_threads: usize, random_seed: usize)
    -> (bitvec::vec::BitVec, Vec<usize>, CompactIntVec, simple_sds_sbwt::bit_vector::BitVector) {

    assert!(key_kmer_marks.len() <= SET_OF_SETS_CONSTRUCTION_MAX_SBWT_LEN);

    log_memory_usage();
    // Build rank support for key k-mer marks
    log::info!("Building rank support for key k-mer marks");
    let key_kmer_marks = crate::util::bitvec_to_simple_sds_raw_bitvec(key_kmer_marks);
    let mut key_kmer_marks = simple_sds_sbwt::bit_vector::BitVector::from(key_kmer_marks);
    key_kmer_marks.enable_rank();
    let n_key_kmers = key_kmer_marks.count_ones();

    log_memory_usage();
    log::info!("Building color set fingerprints");
    // Assign a 128-bit fingerprint for each possible element id. 128-bit integers can not be
    // updated atomically, so instead we use a pair of u64 values which can each be updated atomically.
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(random_seed as u64);
    let element_fingerprints: Vec<(u64,u64)> = (0..n_colors).map(|_i| (rng.next_u64(), rng.next_u64())).collect();

    // 128-bit fingerprints for the color set of each key k-mer. Again we split each u128 into
    // two u64s. We also store the color set size as a third AtomicU64. So we have a triple of
    // u64 for each key k-mer. We store this as a flat vector of u64.
    let mut set_fingerprints_and_sizes = Vec::<AtomicU64>::new();
    set_fingerprints_and_sizes.resize_with(key_kmer_marks.count_ones()*3, || AtomicU64::new(0));
    assert!(set_fingerprints_and_sizes.len() == 3 * key_kmer_marks.rank(key_kmer_marks.len()));

    log_memory_usage();

    let callback = |e: SetElement| {
        if key_kmer_marks.get(e.set_id) {
            let (fp1, fp2) = element_fingerprints[e.color];
            let key_kmer_idx = key_kmer_marks.rank(e.set_id);
            set_fingerprints_and_sizes[key_kmer_idx*3 + 0].fetch_xor(fp1, Release);
            set_fingerprints_and_sizes[key_kmer_idx*3 + 1].fetch_xor(fp2, Release);
            set_fingerprints_and_sizes[key_kmer_idx*3 + 2].fetch_add(1, Release);
        }
    };

    gen.run(callback, n_threads);

    drop(element_fingerprints); // Free memory

    // Turn atomic integers into regular ones. I hope this compiles into just a re-interpretation
    // of the data, not a copy. Also abbreviate the "sfs" for set fingerprints and sizes.
    let mut sfs: Vec<u64> = set_fingerprints_and_sizes.into_iter().map(
        |x| x.load(Acquire)
        // The Acquire here means that the Release-writes above must be visible before these loads (I think).
    ).collect();

    // Pack things. Before:
    // [fp1: u64][fp2: u64][size: u64]
    // After
    // [[fp1: u64][[fp2 : u56][size_high_bits: u8]] [[colex: u40][size_low_bits: u24]]
    log::info!("Bit-packing things");
    for (fp_idx, colex) in key_kmer_marks.one_iter() {
        let size = sfs[fp_idx*3 + 2]; // Assumed to be at most 32 bits
        let size_MSB = size >> 24;

        sfs[fp_idx*3 + 1] &= !(0xFF_u64); // Clear be bits
        sfs[fp_idx*3 + 1] |= size_MSB; // Write the new bits

        sfs[fp_idx*3 + 2] &= !0xFFFFFFFFFF000000_u64; // Clear the 40 most significant bits
        sfs[fp_idx*3 + 2] |= (colex as u64) << 24; // Write the new bits
    }

    log::info!("Sorting by fingerprint");
    assert!(sfs.len() % 3 == 0);
    let triples: &mut [[u64; 3]] = unsafe {
        std::slice::from_raw_parts_mut(
            sfs.as_mut_ptr().cast::<[u64; 3]>(),
            sfs.len() / 3,
        )
    };

    // We use an 8 MiB stack for the thread pool because Rayon seems to sometimes run into a 
    // stack overflow with huge sorts.
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).stack_size(1 << 23).build().unwrap();
    thread_pool.install(|| triples.par_sort_unstable()); // Sorts by fingerprint

    // Now we go from this:
    // [[fp1: u64][[fp2 : u56][size_high_bits: u8]] [[colex: u40][size_low_bits: u24]]
    // To pairs of u64:
    // [[colex: u56][size_high_bits: u8]][[color_set_id: u40][size_low_bits: u24]]
    log::info!("Marking sufficient k-mers");
    let mut sufficient_kmer_marks = bitvec::bitvec![0; key_kmer_marks.len()];
    log_memory_usage();
    let mut color_set_id = 0_usize;
    let mut fp_run_start = 0_usize;
    for fp_idx in 0..n_key_kmers {
        let cur_fp = (sfs[fp_idx*3 + 0], sfs[fp_idx*3 + 1]); // this also includes the 8 MSB bits of the set size -> is ok.
        let size = ((sfs[fp_idx*3 + 1] & 0xFF) << 24) | (sfs[fp_idx*3 + 2] & 0xFFFFFF);
        let colex = sfs[fp_idx*3 + 2] >> 24; // 40 bits

        // Now this is tricky! We're writing a pair over the existing triples.
        // But we're never going read those words again in this loop,
        // so this is ok.
        sfs[fp_idx*2 + 0] = colex << 8; // 56 bits
        sfs[fp_idx*2 + 0] |= size >> 24; // 8 most significant bits of size
        sfs[fp_idx*2 + 1] = (color_set_id as u64) << 24 ; // 40 bits
        sfs[fp_idx*2 + 1] |= size & 0xFFFFFF ; // 24 bits

        if fp_idx + 1 == n_key_kmers || cur_fp != (sfs[(fp_idx+1)*3 + 0], sfs[(fp_idx+1)*3 + 1]) {
            // End of fingerprint run
            let fp_run_end = fp_idx + 1; // Exclusive
            let mut min_colex = usize::MAX;
            for i in fp_run_start..fp_run_end {
                min_colex = std::cmp::min(min_colex, sfs[i*2 + 0] as usize >> 8);
            }
            sufficient_kmer_marks.set(min_colex, true);

            color_set_id += 1;
            fp_run_start = fp_run_end;
        }
    }
    let n_distinct_sets = color_set_id;
    log::info!("{} distinct color sets found", n_distinct_sets);

    // Free memory
    sfs.truncate(n_key_kmers*2); // It has n_key_kmers pairs now
    sfs.shrink_to_fit();
    log_memory_usage();

    log::info!("Sorting by colex");
    let mut colex_set_id_sizes = sfs; // Rename to reflect what it now has
    assert!(colex_set_id_sizes.len() % 2 == 0, "flat pair vector must have even length");
    let pairs: &mut [[u64; 2]] = unsafe {
        std::slice::from_raw_parts_mut(
            colex_set_id_sizes.as_mut_ptr().cast::<[u64; 2]>(),
            colex_set_id_sizes.len() / 2,
        )
    };
    thread_pool.install(|| pairs.par_sort_unstable()); // Sorts by colex because colex is at the MSB bits of the first word

    log::info!("Storing final color set ids and sizes");

    let mut sufficient_set_sizes: Vec<usize> = vec![0; n_distinct_sets];
    // Maps fingerprint-sort-order color set id to colex-order id (matching rank in sufficient_kmer_marks).
    // Populated lazily: representatives always have the minimum colex in their class, so each
    // representative is encountered before any non-representative of the same class.
    let mut fp_id_to_colex_id: Vec<u64> = vec![0; n_distinct_sets];
    log_memory_usage();
    let mut sufficient_kmer_idx = 0_usize;
    for pairs_idx in 0..n_key_kmers {
        let w1 = colex_set_id_sizes[pairs_idx*2 + 0];
        let w2 = colex_set_id_sizes[pairs_idx*2 + 1];
        let colex = w1 >> 8;
        let set_size = ((w1 & 0xFF) << 24) | (w2 & 0x00FFFFFF);
        let fp_set_id = (w2 >> 24) as usize;
        if sufficient_kmer_marks[colex as usize] {
            fp_id_to_colex_id[fp_set_id] = sufficient_kmer_idx as u64;
            sufficient_set_sizes[sufficient_kmer_idx] = set_size as usize;
            sufficient_kmer_idx += 1;
        }
        colex_set_id_sizes[pairs_idx] = fp_id_to_colex_id[fp_set_id]; // Overwriting in-place
    }
    colex_set_id_sizes.truncate(n_key_kmers);
    colex_set_id_sizes.shrink_to_fit();
    log_memory_usage();

    let key_kmer_idx_to_new_id = colex_set_id_sizes; // Rename to reflect what we now have

    // Interpret as usize. If this makes a copy it's ok, it's not the space peak
    let key_kmer_idx_to_new_id = key_kmer_idx_to_new_id.into_iter().map(|x| x as usize).collect();

    // Bit-pack
    log::info!("Bit packing color set ids");
    let key_kmer_idx_to_new_id = CompactIntVec::from_vec(key_kmer_idx_to_new_id);
    log_memory_usage();

    log::info!("Average color set size: {:.2}", (sufficient_set_sizes.iter().sum::<usize>() as f64) / (n_distinct_sets as f64));

    (sufficient_kmer_marks, sufficient_set_sizes, key_kmer_idx_to_new_id, key_kmer_marks)

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
pub fn build_color_set_storage_to_disk<CSS: ColorSetStorage + Send>(colex_sample_marks: bitvec::vec::BitVec, sampled_set_sizes: Vec<usize>, mut element_gens: Vec<(impl crate::set_of_sets_construction::ParallelElementGenerator, Range<usize>)>, out_prefix: &Path, n_threads: usize) {
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

    sizes.into_iter().map(
        |x| x.load(Acquire) as usize
    ).collect()
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
        
        let (new_marks, set_sizes, key_kmer_to_set_id, _key_kmer_marks) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(
            VecVecGenerator::new(sets.clone()),
            bitvec::bitvec![1,1,0,1,1,1], // Mark one of the duplicates as non-key
            sets.len() as u32,
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