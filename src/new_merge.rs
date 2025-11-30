use std::cmp::max;

use sbwt::{dbg::Dbg, merge, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

fn mark_new_key_kmers<CSS: ColorSetStorage>(coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_sbwt_lcs: &LcsArray, n_threads: usize){
    let si = StreamingIndex::new(&merged_sbwt, &merged_sbwt_lcs);

    let k = merged_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert_eq!(k, coloring2.get_k());

    log::info!("Initializing DBG 1");
    let dbg1 = Dbg::new(&coloring1.sbwt(), Some(&coloring1.lcs()), n_threads);

    log::info!("Initializing DBG 2");
    let dbg2 = Dbg::new(&coloring2.sbwt(), Some(&coloring2.lcs()), n_threads);

    let key_kmer_marks = AtomicBitmap::new(merged_sbwt.n_sets());

    log::info!("Iterating unitigs 1");
    dbg1.iter_unitigs_with_callback(|_nodes, unitig|{
        assert!(unitig.len() >= k);
        let last_kmer = &unitig[unitig.len()-k..];
        let merged_colex = merged_sbwt.search(last_kmer);
        let merged_colex = merged_colex.unwrap_or_else(|| panic!("k-mer from DBG1 not found in merged SBWT: {:?}", String::from_utf8_lossy(last_kmer)));
        assert!(merged_colex.len() == 1);
        key_kmer_marks.set(merged_colex.start, true);
    }, n_threads);

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