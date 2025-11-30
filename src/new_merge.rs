use std::{cmp::max, sync::Arc};

use sbwt::{dbg::{Dbg, Node}, merge, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::ColorSetStorage, parallel_ms_iteration::{MergedElementGenerator, MsElementGenerator}, set_of_sets_construction::{build_color_set_storage, find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates}};

fn mark_key_kmers_for<'a, CSS: ColorSetStorage + Send + Sync>(coloring: &'a CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, n_threads: usize) -> Dbg<'a, SubsetMatrix> {

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

    dbg

}

fn mark_new_key_kmers<'a, 'b, CSS: ColorSetStorage + Send + Sync>(coloring1: &'a CompactColexKmers<CSS>, coloring2: &'b CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_dbg: &Dbg<'_, SubsetMatrix>, n_threads: usize) -> (bitvec::vec::BitVec, Dbg<'a, SubsetMatrix>, Dbg<'b, SubsetMatrix>) {
    let k = merged_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert_eq!(k, coloring2.get_k());

    let key_kmer_marks = AtomicBitmap::new(merged_sbwt.n_sets());
    let dbg1 = mark_key_kmers_for(&coloring1, &merged_sbwt, &merged_dbg, &key_kmer_marks, n_threads);
    let dbg2 = mark_key_kmers_for(&coloring2, &merged_sbwt, &merged_dbg, &key_kmer_marks, n_threads);

    (key_kmer_marks.into_bitvec(), dbg1, dbg2)
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
    let merged_dbg = Dbg::new(&merged_sbwt, Some(&merged_sbwt_lcs), n_threads);

    log::info!("=== Phase 1/3: marking new key k-mers ===");
    let (new_key_kmer_marks, dbg1, dbg2) = mark_new_key_kmers(&coloring1, &coloring2, &merged_sbwt, &merged_dbg, n_threads);

    log::info!("=== PHASE 2/3: Building color set finperprints for key k-mers ===");
    let random_seed = 123123; // Todo: be more random
    let gen = MergedElementGenerator {
        merged_sbwt: &merged_sbwt,
        coloring1: &coloring1,
        coloring2: &coloring2,
        dbg1: &dbg1,
        dbg2: &dbg1,
        filter: None,
    };

    let n_colors = coloring1.get_set_storage().n_colors() + coloring2.get_set_storage().n_colors();
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen, new_key_kmer_marks.clone(), merged_sbwt.n_sets(), n_colors, n_threads, random_seed);

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    let gen = MergedElementGenerator {
        merged_sbwt: &merged_sbwt,
        coloring1: &coloring1,
        coloring2: &coloring2,
        dbg1: &dbg1,
        dbg2: &dbg1,
        filter: None,
    };
        
    let css = build_color_set_storage(n_colors, repr_kmer_marks, distinct_set_sizes, gen, n_threads);

    log::info!("Building rank support for key k-mer marks");
    let mut key_kmer_marks = crate::util::bitvec_to_simple_sds_bitvec(new_key_kmer_marks);
    key_kmer_marks.enable_rank();
    assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));

    let merged_sbwt = Arc::new(merged_sbwt);
    let colex_map = ColexToColorSetMap {
        sbwt: merged_sbwt.clone(), // Clones just the Arc
        sampling: key_kmer_marks, 
        color_set_ids: key_kmer_idx_to_set_id,
    };

    let mut color_names = Vec::<String>::new();
    color_names.extend(coloring1.get_color_names().iter().cloned());
    color_names.extend(coloring2.get_color_names().iter().cloned());

    CompactColexKmers::<CSS>::new(merged_sbwt, merged_sbwt_lcs, colex_map, css, Some(&color_names))
}