use std::cmp::max;

use sbwt::{dbg::{Dbg, Node}, merge, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

fn mark_key_kmers_for<CSS: ColorSetStorage + Send + Sync>(coloring: &CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, n_threads: usize) {

    let k = merged_sbwt.k();
    assert_eq!(k, coloring.get_k());

    log::info!("Initializing DBG");
    let dbg = Dbg::new(&coloring.sbwt(), Some(&coloring.lcs()), n_threads);

    log::info!("Iterating unitigs");
    dbg.iter_unitigs_with_callback(|nodes, unitig|{
        assert!(unitig.len() >= k);
        let unitig_colex_ranks = nodes.iter().map(|v| v.id).collect::<Vec<usize>>(); // TODO: avoid this allocation
        let (_, subunitig_ranges) = coloring.break_to_colored_subunitigs(&unitig_colex_ranks, unitig);
        let mut in_neighbor_buf = Vec::<(Node, u8)>::new(); // TODO: avoid this allocation

        for subunitig_range in subunitig_ranges {
            // (s,e) = (start of first k-mer, start of the k-mer after the last k-mer)
            let (s,e) = (subunitig_range.start, subunitig_range.end); 
            assert!(s < e);
            let last_kmer = &unitig[e-1..e-1+k];
            let merged_colex_last = merged_sbwt.search(last_kmer);
            let merged_colex_last = merged_colex_last.unwrap_or_else(|| panic!("k-mer not found in merged SBWT: {:?}", String::from_utf8_lossy(last_kmer)));
            assert!(merged_colex_last.len() == 1);
            key_kmer_marks.set(merged_colex_last.start, true);

            let first_kmer = &unitig[s..s+k];
            let merged_colex_first = merged_sbwt.search(first_kmer);
            let merged_colex_first = merged_colex_first.unwrap_or_else(|| panic!("k-mer not found in merged SBWT: {:?}", String::from_utf8_lossy(first_kmer)));
            assert!(merged_colex_first.len() == 1);
            let merged_colex_first = merged_colex_first.start;

            in_neighbor_buf.clear();
            merged_dbg.push_in_neighbors(Node{id: merged_colex_first}, &mut in_neighbor_buf);
            for (in_node, _) in in_neighbor_buf.iter() {
               key_kmer_marks.set(in_node.id, true);
            }
        }
    }, n_threads);
}

fn mark_new_key_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_dbg: &Dbg<'_, SubsetMatrix>, n_threads: usize) -> bitvec::vec::BitVec {
    let k = merged_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert_eq!(k, coloring2.get_k());

    let key_kmer_marks = AtomicBitmap::new(merged_sbwt.n_sets());
    mark_key_kmers_for(&coloring1, &merged_sbwt, &merged_dbg, &key_kmer_marks, n_threads);
    mark_key_kmers_for(&coloring2, &merged_sbwt, &merged_dbg, &key_kmer_marks, n_threads);

    key_kmer_marks.into_bitvec()
}

pub fn new_merge<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, optimize_peak_ram: bool, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::MergeInterleaving::new(&(*coloring1.sbwt()), &(*coloring2.sbwt()), optimize_peak_ram, n_threads);
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());

    log::info!("Merging SBWTs");
    let sbwt1 = (*coloring1.sbwt()).clone(); // Todo: avoid clone. Currently unavoidable because we have just an Arc to the SBWT, but the merge needs an owned value.
    let sbwt2 = (*coloring2.sbwt()).clone(); // Todo: avoid clone. Currently unavoidable because we have just an Arc to the SBWT, but the merge needs an owned value.
    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);
    let merged_sbwt = sbwt::merge(sbwt1, sbwt2, merge_plan, precalc_len, n_threads);

    log::info!("Building the LCS array for the merged SBWT");
    let merged_sbwt_lcs = LcsArray::from_sbwt(&merged_sbwt, n_threads);

    log::info!("Initializing DBG for the merged SBWT");
    let dbg = Dbg::new(&merged_sbwt, Some(&merged_sbwt_lcs), n_threads);

    let new_key_kmer_marks = mark_new_key_kmers(&coloring1, &coloring2, &merged_sbwt, n_threads);

    todo!();
}