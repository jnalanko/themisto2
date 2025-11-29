use std::{cmp::{max, min}, collections::HashMap, hash::{Hash, Hasher}, sync::Arc};

use bitvec::slice::BitSlice;
use rustc_hash::FxHasher;
use sbwt::{dbg::Node, LcsArray, MergeInterleaving, SbwtIndex, SubsetMatrix, SubsetSeq};

use bitvec::prelude::*;
use simple_sds_sbwt::{ops::{BitVec, Rank}, raw_vector::AccessRaw};

use crate::{colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView}, int_vec::CompactIntVec};

#[derive(Debug, Eq, PartialEq)]
pub struct BitKey<'a> { // Bitslice with a custom hash function
    pub bits: &'a BitSlice,
}

impl std::hash::Hash for BitKey<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash 64 bits at a time
        let len = self.bits.len();
        let n_words = len.div_ceil(64);
        for i in 0..n_words {
            let start = 64*i;
            let end = min(64*(i+1), len);
            let word: u64 = self.bits[start..end].load();
            word.hash(state);
        }
        len.hash(state);  // include length to distinguish e.g. 0b1 from 0b10
    }
}


fn figure_out_if_we_need_to_sample_nonsampled_vs_absent(
    absent_sbwt: &SbwtIndex<SubsetMatrix>, 
    mut absent_colex: usize, // Position in the absent sbwt where k-mer would be inserted
    merged_colex: usize,
    merged_leader_marks: &bitvec::vec::BitVec<u64, Lsb0>,
    absent_merge_marks: &bitvec::vec::BitVec<u64, Lsb0>) -> bool {

    // This node may become the end of a colored unitig in the merged graph. So we may need
    // to sample it. 
    // 
    // This happens if any of the following happen:
    //   (i)   The merged graph has a new outneighbor for this node (unitig ends).
    //   (ii)  The current outneighbor gets a new in-neighbor (unitig ends).
    //   (iii) There will be an edge from the node in the present SBWT to a node
    //         in the absent SBWT. Then the node from the absent SBWT may introduce 
    //         a new color, in which case the colored unitig ends.
    //
    // We assume that all color sets are non-empty, which means that if there is an
    // outedge into the absent sbwt, then this always introduces a new color in case (iii).
    // Under this assumption, if case (i) or (ii) happens, case (iii) also happens, so it's enough
    // to check only for case (iii). If our assumption that all color sets are nonempty
    // does not hold, it only means that we may sample a node unnecessarily, but the
    // color set structure is still correct. 

    assert!(!absent_merge_marks[merged_colex]); // Should be absent
    let mut s = merged_colex;
    while !merged_leader_marks[s] {
        // merged_leader_marks[0] is always set so s > 0 if we are here
        s -= 1;
        if absent_merge_marks[s] {
            absent_colex -= 1;
        }
    }
    let mut e = merged_colex+1;
    while e < merged_leader_marks.len() && !merged_leader_marks[e] {
        e += 1;
    }

    // [s..e) is the suffix group of the present k-mer in the merged sbwt.
    for i in s..e {
        if absent_merge_marks[i] {
            // Suffix group leader in the absent sbwt
            for c_idx in 0..absent_sbwt.alphabet().len() {
                if absent_sbwt.sbwt().set_contains(absent_colex, c_idx as u8) {
                    return true; // Sample x
                }
            }
            return false; // Suffix group leader did not have any edge
        }
    }
    false
}

struct PartitionedIdMap {
    #[allow(clippy::type_complexity)]
    hashmaps: Vec<HashMap::<(Option::<usize>, Option::<usize>), usize>>,
}

struct PartitionedReadOnlyIdMap {
    #[allow(clippy::type_complexity)]
    hashmaps: Vec<HashMap::<(Option::<usize>, Option::<usize>), usize>>,
    cumul_sizes: Vec<usize> // index i contains total length of hash maps [0..i)
}

impl PartitionedIdMap {
    fn hash_pair(x: (Option<usize>, Option<usize>)) -> u64 {
        let mut hasher = FxHasher::default();
        x.hash(&mut hasher);
        hasher.finish()
    }

    fn insert_pair(&mut self, x: (Option<usize>, Option<usize>)) {
        let r = Self::hash_pair(x);
        let hash_map_idx = (r / (u64::MAX / self.hashmaps.len() as u64)) as usize;
        let H = &mut self.hashmaps[hash_map_idx];
        if !H.contains_key(&x) {
            H.insert(x, H.len());
        }
    }

    // There is not method to get a pair. For that, first convert the struct
    // into PartitionedReadOnlyIdMap, which does some precalc to make the
    // lookup faster.
}

impl PartitionedReadOnlyIdMap {
    fn new(old: PartitionedIdMap) -> Self {
        let mut cumul_sizes = Vec::<usize>::with_capacity(old.hashmaps.len() + 1); 
        cumul_sizes.push(0);
        old.hashmaps.iter().for_each(|H| {
            cumul_sizes.push(cumul_sizes.last().unwrap() + H.len());
        });

        PartitionedReadOnlyIdMap{hashmaps: old.hashmaps, cumul_sizes}
    }

    fn total_len(&self) -> usize {
        self.cumul_sizes[self.hashmaps.len()]
    }

    fn get_new_id_of_pair(&self, x: (Option<usize>, Option<usize>)) -> usize {
        let r = PartitionedIdMap::hash_pair(x);
        let hash_map_idx = (r / (u64::MAX / self.hashmaps.len() as u64)) as usize;
        self.cumul_sizes[hash_map_idx] + self.hashmaps[hash_map_idx][&x]
    }

    #[allow(clippy::type_complexity)]
    fn get_old_ids_sorted_by_new_id(&self) -> Vec<(usize, (Option::<usize>, Option::<usize>))> {
        // Collect all elements (new id, old id pair) from the hash maps
        let n_pairs_total = self.total_len();
        let mut id_pairs_in_new_id_order = self.hashmaps.iter().fold(
            Vec::<(usize, (Option::<usize>, Option::<usize>))>::with_capacity(n_pairs_total),
            |mut acc, H| {
                let len_before = acc.len();
                acc.extend(
                    H.iter().map(|(pair, new_id)| (*new_id + len_before, *pair))
                );
                acc
            }
        );
        id_pairs_in_new_id_order.sort();
        id_pairs_in_new_id_order
    }

}


fn compute_color_id_pairs_and_merged_unitig_sampling<CSS: ColorSetStorage>(coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, lcs1: &LcsArray, lcs2: &LcsArray, merge_plan: &MergeInterleaving, n_threads: usize) -> (PartitionedReadOnlyIdMap, simple_sds_sbwt::raw_vector::RawVector) {
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();

    // Distinct color id pairs, inserted into key-disjoint hash maps. The values are
    // color set ids within the hash map.
    let hashmaps = vec![HashMap::<(Option::<usize>, Option::<usize>), usize>::new(); n_threads];
    let mut new_id_map = PartitionedIdMap{hashmaps};

    let mut colex1 = 0_usize;
    let mut colex2 = 0_usize;

    let mut color_set_sample_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(merged_len, false);
    log::info!("Building DBG support");
    let dbg1 = sbwt::dbg::Dbg::new(&(*coloring1.sbwt()), Some(lcs1), n_threads);
    let dbg2 = sbwt::dbg::Dbg::new(&(*coloring2.sbwt()), Some(lcs2), n_threads);
    let mut outlabel_buf_1 = Vec::<u8>::new();
    let mut outlabel_buf_2 = Vec::<u8>::new();

    log::info!("Computing new color set id pairs and merged unitig sampling");
    #[derive(Debug)]
    enum Case { // Three cases in a loop below
        Sampled(usize),
        NotSampled,
        Absent,
    }
    for merged_colex in 0..merged_len {
        if !merge_plan.is_dummy[merged_colex] {
            let c1 = if !merge_plan.s1[merged_colex] {
                Case::Absent
            } else if coloring1.get_map().sampling.get(colex1) {
                Case::Sampled(coloring1.colex_to_set_id(colex1))
            } else {
                Case::NotSampled
            };

            let c2 = if !merge_plan.s2[merged_colex] {
                Case::Absent
            } else if coloring2.get_map().sampling.get(colex2) {
                Case::Sampled(coloring2.colex_to_set_id(colex2))
            } else {
                Case::NotSampled
            };

            // Ok, this is going to get a bit verbose but bear with me. We have
            // 3 * 3 = 9 cases. There are two symmetric pairs of cases and three unique cases. We could
            // reduce code duplication by making symmetric cases call a common function,
            // but it's so few lines of code anyway so let's just go with this.
            match (c1, c2) {
                (Case::Sampled(id1), Case::Sampled(id2)) => {
                    new_id_map.insert_pair((Some(id1), Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::Sampled(id1), Case::NotSampled) => {
                    let id2 = coloring2.colex_to_set_id(colex2);
                    new_id_map.insert_pair((Some(id1), Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::Sampled(id1), Case::Absent) => {
                    new_id_map.insert_pair((Some(id1), None));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::NotSampled, Case::Sampled(id2)) => {
                    //eprintln!("Case 3");
                    let id1 = coloring1.colex_to_set_id(colex1);
                    new_id_map.insert_pair((Some(id1), Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::NotSampled, Case::NotSampled) => {
                    // K-mer is in both SBWTs but its not sampled in either one.
                    // Since it is not sampled in either SBWT, the outdegree of this k-mer
                    // is 1 in both. But we might still need to sample it in the merged graph.
                    // There are two cases:
                    // 1) The outneighbor k-mers are the same k-mer. Then the outdegree in the merged graph 
                    //    will be 1, and that outneighbor will have the same color set id pair as this
                    //    one -> this node does not need to be sampled
                    // 2) The outneighbor k-mers are different. Now we have a new outgoing branch at this 
                    //    node. Which means this node needs to be sampled.

                    outlabel_buf_1.clear();
                    outlabel_buf_2.clear();
                    dbg1.push_outlabels(Node{id: colex1}, &mut outlabel_buf_1);
                    dbg2.push_outlabels(Node{id: colex2}, &mut outlabel_buf_2);
                    assert_eq!(outlabel_buf_1.len(), 1);
                    assert_eq!(outlabel_buf_2.len(), 1);
                    //eprintln!("{} {}", *outlabel_buf_1.first().unwrap() as char, *outlabel_buf_2.first().unwrap() as char);
                    match (outlabel_buf_1.first(), outlabel_buf_2.first()) {
                        (Some(a), Some(b)) => {
                            if a != b { // Case 2 in the comment above
                                color_set_sample_marks.set_bit(merged_colex, true);
                                let id1 = coloring1.colex_to_set_id(colex1);
                                let id2 = coloring2.colex_to_set_id(colex2);
                                new_id_map.insert_pair((Some(id1), Some(id2)));
                            } else { // The else-branch would be case 1 but then there is nothing to do

                            }
                        }
                        _ => panic!("Bug at computing color set samples bit vector in merge") // Both should have outdegree > 0
                    }
                },
                (Case::NotSampled, Case::Absent) => {
                    let id1 = coloring1.colex_to_set_id(colex1);
                    new_id_map.insert_pair((Some(id1), None));
                    if figure_out_if_we_need_to_sample_nonsampled_vs_absent(&coloring2.sbwt(), colex2, merged_colex, &merge_plan.is_leader, &merge_plan.s2) {
                        color_set_sample_marks.set_bit(merged_colex, true);
                    }
                },
                (Case::Absent, Case::Sampled(id2)) => {
                    new_id_map.insert_pair((None, Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::Absent, Case::NotSampled) => {
                    let id2 = coloring2.colex_to_set_id(colex2);
                    new_id_map.insert_pair((None, Some(id2)));
                    if figure_out_if_we_need_to_sample_nonsampled_vs_absent(&coloring1.sbwt(), colex1, merged_colex, &merge_plan.is_leader, &merge_plan.s1) {
                        color_set_sample_marks.set_bit(merged_colex, true);
                    }
                },
                (Case::Absent, Case::Absent) => panic!("Nonexisting merged kmer") // Impossible
            }
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s2[merged_colex] as usize;
    }

    (PartitionedReadOnlyIdMap::new(new_id_map), color_set_sample_marks)

}

struct TwoSetMerger<L: Iterator<Item = usize>, R: Iterator<Item = usize>> {
    left: Option<L>,
    right: Option<R>,
    left_n_colors: usize, // The left set get colors 0..left_n_colors, the right set gets left_n_colors..
}

impl<L: Iterator<Item = usize>, R: Iterator<Item = usize>> Iterator for TwoSetMerger<L,R> {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        // Terrible branch city. TODO: do better.

        // Try to take from left
        if let Some(l) = &mut self.left {
            if let Some(x) = l.next() {
                return Some(x);
            } 
        }

        // Could not take from left -> take from right
        if let Some(r) = &mut self.right {
            r.next().map(|x| self.left_n_colors + x)
        } else {
            None // Finished
        }
    }
}


fn encode_merged_color_sets<CSS: ColorSetStorage>(new_id_map: &PartitionedReadOnlyIdMap, coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>) -> CSS {

    let n_colors_1 = coloring1.get_set_storage().get_full_set().iter().count();
    let n_colors_2 = coloring2.get_set_storage().get_full_set().iter().count();
    let id_pairs_in_new_id_order = new_id_map.get_old_ids_sorted_by_new_id();

    // Create an iterator of combined sets
    let mut pair_id = 0_usize;
    let n_pairs = id_pairs_in_new_id_order.len();
    let n_colors_1_ref = &n_colors_1; // Reference to move by reference into the closure 
    let iter_of_iters = std::iter::from_fn(move || {
        if pair_id == n_pairs {
            None
        } else {
            let (_, (left, right)) = id_pairs_in_new_id_order[pair_id]; 
            pair_id += 1;

            match (left,right) {
                (Some(x), Some(y)) => {
                    let set1 = coloring1.set_id_to_set(x);
                    let set2 = coloring2.set_id_to_set(y);
                    Some(TwoSetMerger{left: Some(set1.iter()), right: Some(set2.iter()), left_n_colors: *n_colors_1_ref})
                },
                (Some(x), None) => {
                    let set1 = coloring1.set_id_to_set(x);
                    Some(TwoSetMerger{left: Some(set1.iter()), right: None, left_n_colors: *n_colors_1_ref})
                },
                (None, Some(y)) => {
                    let set2 = coloring2.set_id_to_set(y);
                    Some(TwoSetMerger{left: None, right: Some(set2.iter()), left_n_colors: *n_colors_1_ref})
                }
                (None, None) => panic!("Nonexisting color set id pair")
            }
        }
    });

    *CSS::new_from_iter_of_iters(iter_of_iters, n_colors_1 + n_colors_2)

}

fn store_new_sampled_color_ids<CSS: ColorSetStorage>(n_distinct_color_sets: usize, merge_plan: &MergeInterleaving, color_set_sample_marks: &simple_sds_sbwt::bit_vector::BitVector, coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, pair_to_new_id_maps: &PartitionedReadOnlyIdMap) -> CompactIntVec {
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();

    let bits_per_color_set_id = n_distinct_color_sets.next_power_of_two().trailing_zeros() as usize;
    let mut sampled_ids = CompactIntVec::new(color_set_sample_marks.count_ones(), bits_per_color_set_id);
    let mut n_items_pushed = 0_usize;
    let mut colex1 = 0_usize;
    let mut colex2 = 0_usize;
    for merged_colex in 0..merged_len {
        if color_set_sample_marks.get(merged_colex) {
            let color_set_id_1 = if merge_plan.s1[merged_colex] {
                Some(coloring1.colex_to_set_id(colex1))
            } else {
                None
            };

            let color_set_id_2 = if merge_plan.s2[merged_colex] {
                Some(coloring2.colex_to_set_id(colex2))
            } else {
                None
            };

            // The merge plan should not have a zero-bit at the same position in s1 and s2
            assert!(color_set_id_1.is_some() || color_set_id_2.is_some());
            let id = pair_to_new_id_maps.get_new_id_of_pair((color_set_id_1, color_set_id_2));
            sampled_ids.set(n_items_pushed, id);
            n_items_pushed += 1;
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s2[merged_colex] as usize;
    }

    sampled_ids
}

pub fn merge_compact_colorings<CSS: ColorSetStorage>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, optimize_peak_ram: bool, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::MergeInterleaving::new(&(*coloring1.sbwt()), &(*coloring2.sbwt()), optimize_peak_ram, n_threads);

    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();    

    log::info!("Computing color id pairs and merged sampling");
    let (new_id_map, color_set_sample_marks) = compute_color_id_pairs_and_merged_unitig_sampling(&coloring1, &coloring2, &coloring1.lcs(), &coloring2.lcs(), &merge_plan, n_threads);

    let mut color_set_sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(color_set_sample_marks);
    color_set_sample_marks.enable_rank();
    let n_sampled = color_set_sample_marks.rank(color_set_sample_marks.len());
    log::info!("Sampled {} out of {} SBWT nodes ({:.2}%)", n_sampled, merged_len, n_sampled as f64 / merged_len as f64 * 100.0);

    log::info!("Encoding distinct merged color sets");
    let css = encode_merged_color_sets(&new_id_map, &coloring1, &coloring2);

    log::info!("Storing new sampled color set ids");
    let n_distinct_color_sets = new_id_map.total_len(); 
    let sampled_ids = store_new_sampled_color_ids(n_distinct_color_sets, &merge_plan, &color_set_sample_marks, &coloring1, &coloring2, &new_id_map);

    log::info!("Interleaving SBWTs");
    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);

    // Collect old color names before dropping the structs
    let mut new_color_names = coloring1.get_color_names().clone();
    new_color_names.extend(coloring2.get_color_names().clone());

    let sbwt1 = (*coloring1.sbwt()).clone(); // Todo: avoid clone. Currently unavoidable because we have just a reference to the SBWT, but the merge needs an owned value.
    drop(coloring1);

    let sbwt2 = (*coloring2.sbwt()).clone(); // Todo: avoid clone
    drop(coloring2);

    log::info!("Interleaving SBWTs");
    let merged_sbwt = Arc::new(sbwt::merge(sbwt1, sbwt2, merge_plan, precalc_len, n_threads));

    log::info!("Computing the merged LCS array"); // Todo: could we do this during the interleave?
    let merged_lcs = LcsArray::from_sbwt(&merged_sbwt, n_threads);

    let new_coloring = CompactColexKmers::new(
        merged_sbwt.clone(),
        merged_lcs,
        ColexToColorSetMap {
            sbwt: merged_sbwt.clone(),
            sampling: color_set_sample_marks,
            color_set_ids: sampled_ids
        },
        css,
        Some(&new_color_names)
    );

    log::info!("Color merge finished");
    new_coloring

}


#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use jseqio::seq_db::SeqDB;
    use sbwt::{BitPackedKmerSortingMem, LcsArray, SbwtIndex, SubsetMatrix};
    use simple_sds_sbwt::ops::{BitVec, Rank};

    use crate::{bitmap_storage::build_from_seq_dbs, colex_colored_kmers::{ColexToColorSetMap, hash_and_encode_distinct_sets, mark_key_kmers}, coloring_interface::{ColorSetStorage, ColorSetView}, int_vec::CompactIntVec, sparse_dense_storage::SparseDenseStorage, util::VecVecSeqStream};

    use super::CompactColexKmers;


    #[cfg(test)]
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
        let bms1 = build_from_seq_dbs(dbs1, &sbwt1, &lcs1, n_threads);

        let iter_of_iters_1 = (0..sbwt1.n_sets()).into_iter().map(|colex| bms1.get_set_view(colex).iter());
        let colex_to_css_1 = *CSS::new_from_iter_of_iters(iter_of_iters_1, n_colors_1);

        let (distinct_css_1, set_to_id_1) = hash_and_encode_distinct_sets(&colex_to_css_1, n_colors_1);
        let colex_to_id: Vec<usize> = (0..sbwt1.n_sets()).into_iter().map(|colex| {
            set_to_id_1[&colex_to_css_1.get_set_view(colex)]
        }).collect(); 

        (colex_to_id, distinct_css_1)
    }

    #[test]
    fn test_merge() {

        if std::env::var("RUST_LOG").is_err() {
            std::env::set_var("RUST_LOG", "info")
        }
        env_logger::init();

        let n_threads = 3;

        for k in 3_usize..10_usize { // k < 3 does not work because construction uses 3-mer binning.

            let input_seqs_1: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (i + k.pow(4)) as u64)).collect();
            let input_seqs_2: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (123456 + i + k.pow(4)) as u64)).collect();

            let mut all_input_seq_slices = Vec::<&[u8]>::new();
            all_input_seq_slices.extend(input_seqs_1.iter().map(|s| s.as_slice()));
            all_input_seq_slices.extend(input_seqs_2.iter().map(|s| s.as_slice()));

            let mut all_input_seqs: Vec<Vec<u8>> = all_input_seq_slices.iter().map(|s| s.to_vec()).collect();

            let mut dbs1 = Vec::<SeqDB>::new();
            let mut dbs2 = Vec::<SeqDB>::new();
            let mut dbs_both = Vec::<SeqDB>::new();
            for seq in input_seqs_1.iter() {
                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs1.push(db);

                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs_both.push(db);
            }
            for seq in input_seqs_2.iter() {
                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs2.push(db);

                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs_both.push(db);
            }

            let (mut sbwt1, lcs1) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(false)
                .k(k)
                .build_lcs(true)
                .n_threads(3)
                .precalc_length(5)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
            .run_from_vecs(&input_seqs_1);

            let (mut sbwt2, lcs2) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(false)
                .k(k)
                .build_lcs(true)
                .n_threads(3)
                .precalc_length(5)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
            .run_from_vecs(&input_seqs_2);

            let (mut sbwt_both, lcs_both) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(false)
                .k(k)
                .build_lcs(true)
                .n_threads(3)
                .precalc_length(5)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
            .run_from_slices(&all_input_seq_slices);

            sbwt1.build_select();
            sbwt2.build_select();
            sbwt_both.build_select();

            let sbwt1 = Arc::new(sbwt1);
            let sbwt2 = Arc::new(sbwt2);
            let sbwt_both = Arc::new(sbwt_both);

            let lcs1 = lcs1.unwrap();
            let lcs2 = lcs2.unwrap();
            let lcs_both = lcs_both.unwrap();


            let sample_distance = 3;

            let (colex_to_id_1, storage_1) = build_color_sets::<SparseDenseStorage>(&sbwt1, &lcs1, dbs1, n_threads); 
            let (colex_to_id_2, storage_2) = build_color_sets::<SparseDenseStorage>(&sbwt2, &lcs2, dbs2, n_threads); 
            let (colex_to_id_both, storage_both)= build_color_sets::<SparseDenseStorage>(&sbwt_both, &lcs_both, dbs_both, n_threads); 
            
            let key_kmers_1 = mark_key_kmers(&sbwt1, &lcs1, sample_distance, VecVecSeqStream::new(input_seqs_1.clone()), n_threads);
            let key_kmers_2 = mark_key_kmers(&sbwt2, &lcs2, sample_distance, VecVecSeqStream::new(input_seqs_2.clone()), n_threads);
            let key_kmers_both = mark_key_kmers(&sbwt_both, &lcs_both, sample_distance, VecVecSeqStream::new(all_input_seqs.clone()), n_threads);

            let sampled_ids_1: Vec<usize> = colex_to_id_1.iter().enumerate().filter(|(i, _)| key_kmers_1[*i]).map(|(_,x)| *x).collect();
            let sampled_ids_2: Vec<usize> = colex_to_id_2.iter().enumerate().filter(|(i, _)| key_kmers_2[*i]).map(|(_,x)| *x).collect();
            let sampled_ids_both: Vec<usize> = colex_to_id_both.iter().enumerate().filter(|(i, _)| key_kmers_both[*i]).map(|(_,x)| *x).collect();

            assert!(key_kmers_1.count_ones() == sampled_ids_1.len());
            let mut key_kmers_1 = crate::util::bitvec_to_simple_sds_bitvec(key_kmers_1);
            let mut key_kmers_2 = crate::util::bitvec_to_simple_sds_bitvec(key_kmers_2);
            let mut key_kmers_both = crate::util::bitvec_to_simple_sds_bitvec(key_kmers_both);

            key_kmers_1.enable_rank();
            key_kmers_2.enable_rank();
            key_kmers_both.enable_rank();

            let colex_map_1 = ColexToColorSetMap{
                sbwt: sbwt1.clone(),
                sampling: key_kmers_1,
                color_set_ids: CompactIntVec::from_vec(sampled_ids_1),
            };

            let colex_map_2 = ColexToColorSetMap{
                sbwt: sbwt2.clone(),
                sampling: key_kmers_2,
                color_set_ids: CompactIntVec::from_vec(sampled_ids_2),
            };

            let colex_map_both = ColexToColorSetMap{
                sbwt: sbwt_both.clone(),
                sampling: key_kmers_both,
                color_set_ids: CompactIntVec::from_vec(sampled_ids_both),
            };

            let ccc1 = CompactColexKmers::new(sbwt1, lcs1, colex_map_1, storage_1, None);
            let ccc2 = CompactColexKmers::new(sbwt2, lcs2, colex_map_2, storage_2, None);
            let ccc_both = CompactColexKmers::new(sbwt_both, lcs_both, colex_map_both, storage_both, None);

            let ccc_merged = super::merge_compact_colorings(ccc1, ccc2, true, n_threads);
            let sbwt_merged = &ccc_merged.sbwt();

            for colex in 0..ccc_both.sbwt().n_sets() {
                let kmer = ccc_both.sbwt().access_kmer(colex);

                if kmer.iter().all(|c| *c != b'$') { // Not a dummy k-mer
                    let true_colors: Vec<usize> = ccc_both.colex_to_set(colex).iter().collect();
                    let range = sbwt_merged.search(&kmer).unwrap();
                    assert_eq!(range.len(), 1);
                    let colex_merged = range.start;
                    //let merged_colors = ccc_merged.colex_to_set(colex_merged).as_bitvec(ccc_both.n_colors);
                    let merged_colors: Vec<usize> = ccc_merged.colex_to_set(colex_merged).iter().collect();

                    eprintln!("{} {} {:?} {:?} {} {}", colex, String::from_utf8_lossy(&kmer), true_colors, sbwt_merged.search(&kmer), ccc_merged.get_map().sampling.get(colex_merged), ccc_merged.colex_to_set_id(colex_merged));
                    assert_eq!(true_colors, merged_colors);
                }

            }
        }
    }
}