use std::cmp::max;

use sbwt::{dbg::Dbg, merge, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

fn mark_key_kmers_for<CSS: ColorSetStorage + Send + Sync>(coloring: &CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, key_kmer_marks: &AtomicBitmap, n_threads: usize) {

    let k = merged_sbwt.k();
    assert_eq!(k, coloring.get_k());

    log::info!("Initializing DBG");
    let dbg = Dbg::new(&coloring.sbwt(), Some(&coloring.lcs()), n_threads);

    log::info!("Iterating unitigs");
    dbg.iter_unitigs_with_callback(|nodes, unitig|{
        assert!(unitig.len() >= k);
        let unitig_colex_ranks = nodes.iter().map(|v| v.id).collect::<Vec<usize>>(); // TODO: avoid this allocation
        let (_, subunitig_ranges) = coloring.break_to_colored_subunitigs(&unitig_colex_ranks, unitig);

        for subunitig_range in subunitig_ranges {
            // (s,e) = (start of first k-mer, start of the k-mer after the last k-mer)
            let (s,e) = (subunitig_range.start, subunitig_range.end); 
            assert!(s < e);
            let last_kmer = &unitig[e-1..e-1+k];
            let merged_colex = merged_sbwt.search(last_kmer);
            let merged_colex = merged_colex.unwrap_or_else(|| panic!("k-mer from DBG1 not found in merged SBWT: {:?}", String::from_utf8_lossy(last_kmer)));
            assert!(merged_colex.len() == 1);
            key_kmer_marks.set(merged_colex.start, true);
        }
    }, n_threads);
}

fn mark_new_key_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, n_threads: usize) -> bitvec::vec::BitVec {
    let k = merged_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert_eq!(k, coloring2.get_k());

    let key_kmer_marks = AtomicBitmap::new(merged_sbwt.n_sets());
    mark_key_kmers_for(&coloring1, &merged_sbwt, &key_kmer_marks, n_threads);
    mark_key_kmers_for(&coloring2, &merged_sbwt, &key_kmer_marks, n_threads);

    key_kmer_marks.into_bitvec()
}

pub fn new_merge<CSS: ColorSetStorage>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, optimize_peak_ram: bool, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::MergeInterleaving::new(&(*coloring1.sbwt()), &(*coloring2.sbwt()), optimize_peak_ram, n_threads);
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());

    log::info!("Merging SBWTs");
    let sbwt1 = (*coloring1.sbwt()).clone(); // Todo: avoid clone. Currently unavoidable because we have just an Arc to the SBWT, but the merge needs an owned value.
    let sbwt2 = (*coloring2.sbwt()).clone(); // Todo: avoid clone. Currently unavoidable because we have just an Arc to the SBWT, but the merge needs an owned value.
    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);
    let merged_sbwt = sbwt::merge(sbwt1, sbwt2, merge_plan, precalc_len, n_threads);

    let merged_sbwt_len = LcsArray::from_sbwt(&merged_sbwt, n_threads);

    todo!();
}