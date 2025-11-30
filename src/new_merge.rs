use std::{cmp::max, sync::Arc};

use sbwt::{dbg::{Dbg, Node}, merge, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::ColorSetStorage, parallel_ms_iteration::{MergedElementGenerator, MsElementGenerator}, set_of_sets_construction::{build_color_set_storage, find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates}};

fn mark_kmer(kmer: &[u8], sbwt: &SbwtIndex<SubsetMatrix>, marks: &AtomicBitmap) {
    let colex = sbwt.search(kmer);
    let colex = colex.unwrap_or_else(|| panic!("k-mer not found in merged SBWT: {:?}", String::from_utf8_lossy(kmer)));
    assert!(colex.len() == 1);
    marks.set(colex.start, true);
}

fn mark_in_neighbors<'a>(kmer: &[u8], sbwt: &'a SbwtIndex<SubsetMatrix>, dbg: &Dbg<'a, SubsetMatrix>, marks: &AtomicBitmap) {
    let colex = sbwt.search(kmer);
    let colex = colex.unwrap_or_else(|| panic!("k-mer not found in merged SBWT: {:?}", String::from_utf8_lossy(kmer)));
    assert!(colex.len() == 1);
    let colex = colex.start;

    let mut in_neighbor_buf = Vec::<(Node, u8)>::new(); // TODO: avoid this allocation
    dbg.push_in_neighbors(Node{id: colex}, &mut in_neighbor_buf);
    for (in_node, _) in in_neighbor_buf.iter() {
        marks.set(in_node.id, true);
    }
}

fn mark_key_kmers_for<'a, CSS: ColorSetStorage + Send + Sync>(coloring: &'a CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, visited_kmers_before: &AtomicBitmap, n_threads: usize) -> Dbg<'a, SubsetMatrix> {

    let k = merged_sbwt.k();
    assert_eq!(k, coloring.get_k());

    log::info!("Initializing DBG");
    let dbg = Dbg::new(&coloring.sbwt(), Some(&coloring.lcs()), n_threads);
    let visited_kmers_now = AtomicBitmap::new(merged_sbwt.n_sets());

    log::info!("Iterating unitigs");
    dbg.iter_unitigs_with_callback(|nodes, unitig|{
        assert!(unitig.len() >= k);
        let unitig_colex_ranks = nodes.iter().map(|v| v.id).collect::<Vec<usize>>(); // TODO: avoid this allocation
        let (_, subunitig_ranges) = coloring.break_to_colored_subunitigs(&unitig_colex_ranks, unitig);

        let mut prev_was_visited = false;

        for kmer_start in 0..nodes.len() {
            let kmer_colex = nodes[kmer_start].id;
            let visited = visited_kmers_before.get(kmer_colex);
            if visited {
                if !prev_was_visited {
                    // Start of a new colored subunitig
                    // -> Mark all in-neighbors for sampling
                    todo!();
                } else {
                    // Extending the colored subunitig -> no need to mark
                }
            } else { // Not visited
                if prev_was_visited {
                    // One past the end of a color subunitig
                    // -> mark previous node for sampling
                    todo!();
                } else{
                    // Extending a colored subunitig without interference 
                    // from previous colors -> no need to mark
                }
            }
            prev_was_visited = visited;
        }

        for subunitig_range in subunitig_ranges {
            // (s,e) = (start of first k-mer, start of the k-mer after the last k-mer)
            let (s,e) = (subunitig_range.start, subunitig_range.end); 
            assert!(s < e);
            let last_kmer = &unitig[e-1..e-1+k];

            let first_kmer = &unitig[s..s+k];
        }
    }, n_threads);

    todo!(); // Update visited marks

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
    let mut merged_sbwt = sbwt::merge(sbwt1, sbwt2, merge_plan, precalc_len, n_threads);
    merged_sbwt.build_select();

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
        merged_lcs: &merged_sbwt_lcs,
        coloring1: &coloring1,
        coloring2: &coloring2,
        dbg1: &dbg1,
        dbg2: &dbg2,
        filter: None,
    };

    let n_colors = coloring1.get_set_storage().n_colors() + coloring2.get_set_storage().n_colors();
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen, new_key_kmer_marks.clone(), merged_sbwt.n_sets(), n_colors, n_threads, random_seed);

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    let gen = MergedElementGenerator {
        merged_sbwt: &merged_sbwt,
        merged_lcs: &merged_sbwt_lcs,
        coloring1: &coloring1,
        coloring2: &coloring2,
        dbg1: &dbg1,
        dbg2: &dbg2,
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
            //let sample_distance = 1;

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

            let ccc_merged = super::new_merge(ccc1, ccc2, true, n_threads);
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