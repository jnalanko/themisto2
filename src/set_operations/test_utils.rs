//! Helpers for the tests of set operations: building small colored indexes from named sequences.

use std::{collections::HashMap, hash::BuildHasherDefault};

use jseqio::seq_db::SeqDB;
use rustc_hash::FxHasher;
use sbwt::{BitPackedKmerSortingMem, LcsArray, SbwtIndex, SubsetMatrix};
use simple_sds_sbwt::ops::Rank;

use crate::{colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers, mark_key_kmers}, coloring_interface::{ColorSetStorage, ColorSetView}, int_vec::CompactIntVec, io::RewindableSeqStreamGenerator, sparse_dense_storage::SparseDenseStorage, util::VecVecRewindableGen};

// A sample distance larger than any unitig in the tests, so that only the structurally
// required key k-mers get marked. Any k-mer that a set operation forgets to mark then results
// in a wrong color set.
pub(crate) const NO_SAMPLING: usize = 1000;

/// Output:
/// - Distinct color sets encoded as ColorSetStorage
/// - HashMap from color set to its index in ColorSets
fn hash_and_encode_distinct_sets<'a, CSS: ColorSetStorage>(colex_to_set: &'a CSS, n_colors: usize) -> (CSS, HashMap::<CSS::SetView<'a>, usize, BuildHasherDefault::<FxHasher>>) {
    let n_sets = colex_to_set.n_sets();

    log::info!("Hashing distinct color sets");

    let mut distinct_sets = HashMap::<CSS::SetView<'a>, usize, BuildHasherDefault::<FxHasher>>::default(); // Set -> id
    let mut distinct_set_colex_ranks = Vec::<usize>::new();
    let bar = indicatif::ProgressBar::new(n_sets as u64);
    for colex in 0..n_sets {
        let key = colex_to_set.get_set_view(colex);
        if !distinct_sets.contains_key(&key) {
            distinct_sets.insert(key, distinct_sets.len());
            distinct_set_colex_ranks.push(colex);
        }
        if colex % 1000 == 0 {
            bar.inc(1000);
        }
    }
    bar.finish();

    log::info!("{} distinct color sets found", distinct_sets.len());

    // Create an iterator of iterators, each inner iterator iterating over one color set
    let color_sets_iterator = distinct_set_colex_ranks.into_iter().map(|colex| {
        colex_to_set.get_set_view(colex).iter()
    });

    let colorsets = CSS::new_from_iter_of_iters(color_sets_iterator, n_colors);

    (*colorsets, distinct_sets)

}

pub(crate) fn gen_random_dna_string(len: usize, seed: u64) -> Vec<u8> {
    use rand_chacha::rand_core::{RngCore, SeedableRng};

    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(seed);
    (0..len).map(|_| { 
        match rng.next_u64() % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("Impossible")
        }
    }).collect()
}

fn build_color_sets<CSS: ColorSetStorage>(sbwt1: &SbwtIndex<SubsetMatrix>, lcs1: &LcsArray, dbs1: Vec<SeqDB>, n_threads: usize) 
-> (Vec<usize>, CSS){
    let n_colors_1 = dbs1.len();
    let bms1 = crate::bitmap_storage::build_from_seq_dbs(dbs1, sbwt1, lcs1, n_threads);

    let iter_of_iters_1 = (0..sbwt1.n_sets()).map(|colex| bms1.get_set_view(colex).iter());
    let colex_to_css_1 = *CSS::new_from_iter_of_iters(iter_of_iters_1, n_colors_1);

    let (distinct_css_1, set_to_id_1) = hash_and_encode_distinct_sets(&colex_to_css_1, n_colors_1);
    let colex_to_id: Vec<usize> = (0..sbwt1.n_sets()).map(|colex| {
        set_to_id_1[&colex_to_css_1.get_set_view(colex)]
    }).collect(); 

    (colex_to_id, distinct_css_1)
}

/// Turns (color name, sequence) pairs into (color id, sequence) pairs, where the color id
/// is the index of the name in `names`. Names that are not yet in `names` are appended, so
/// a name that appears several times is one color with several sequences.
pub(crate) fn assign_color_ids(names: &mut Vec<String>, colors: &[(&str, Vec<u8>)]) -> Vec<(usize, Vec<u8>)> {
    colors.iter().map(|(name, seq)| {
        let color = names.iter().position(|x| x == name).unwrap_or_else(|| {
            names.push(name.to_string());
            names.len() - 1
        });
        (color, seq.clone())
    }).collect()
}

fn build_sbwt(k: usize, seqs: &[Vec<u8>]) -> (SbwtIndex<SubsetMatrix>, LcsArray) {
    let (mut sbwt, lcs) = BitPackedKmerSortingMem::new_from_vecs(seqs, k)
        .add_rev_comp(false)
        .build_lcs(true)
        .n_threads(3)
        .precalc_length(5)
        .dedup_batches(true)
        .run();
    sbwt.build_select();
    (sbwt, lcs.unwrap())
}

/// Builds a CompactColexKmers where color i is named names[i] and consists of the sequences
/// s with (i, s) in seqs.
pub(crate) fn build_named_coloring(k: usize, names: &[String], seqs: &[(usize, Vec<u8>)], sample_distance: usize, n_threads: usize) -> CompactColexKmers<SparseDenseStorage> {
    let mut dbs: Vec<SeqDB> = names.iter().map(|_| SeqDB::new()).collect();
    for (color, seq) in seqs {
        dbs[*color].push_seq(seq);
    }
    let all_seqs: Vec<Vec<u8>> = seqs.iter().map(|(_, seq)| seq.clone()).collect();
    let (sbwt, lcs) = build_sbwt(k, &all_seqs);
    let (colex_to_id, storage) = build_color_sets::<SparseDenseStorage>(&sbwt, &lcs, dbs, n_threads);

    let mut gen: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(VecVecRewindableGen::new(all_seqs));
    let dbg = sbwt::dbg::Dbg::new(&sbwt, Some(&lcs), n_threads);
    let key_kmers = mark_key_kmers(&sbwt, &dbg, sample_distance, &mut gen, n_threads, 1, false);
    let sampled_ids: Vec<usize> = colex_to_id.iter().enumerate().filter(|(i, _)| key_kmers[*i]).map(|(_,x)| *x).collect();
    assert!(key_kmers.count_ones() == sampled_ids.len());

    let mut key_kmers = crate::util::bitvec_to_simple_sds_bitvec(key_kmers);
    key_kmers.enable_rank();

    let colex_map = ColexToColorSetMap{
        sampling: key_kmers,
        color_set_ids: CompactIntVec::from_vec(sampled_ids),
    };

    CompactColexKmers::new(sbwt, lcs, colex_map, storage, Some(names))
}
