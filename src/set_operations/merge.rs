use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};
use crate::set_operations::{color_sets::ColorCombination, colors::merged_color_mapping, pipeline::{set_operation, SbwtOperation}};

/// Merges two colored indexes: the result has the k-mers of both indexes, and the color set of
/// each k-mer is the union of its color sets in the two indexes. If merge_shared_colors is true,
/// colors with the same name in the two indexes are the same color of the result.
///
/// If keep_redundant_dummies is true, the dummy nodes that become redundant in the SBWT merge are
/// kept. The result has the same k-mers and colors, but possibly more dummy nodes.
pub fn merge_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, merge_shared_colors: bool, optimize_peak_ram: bool, keep_redundant_dummies: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {
    let (color2_to_merged, merged_color_names) = merged_color_mapping(coloring1.get_color_names(), coloring2.get_color_names(), merge_shared_colors).unwrap_or_else(|e| {
        log::error!("{}", e);
        panic!("{}", e);
    });
    let colors = ColorCombination::Union { color2_to_result: &color2_to_merged };
    set_operation(coloring1, coloring2, SbwtOperation::Merge { keep_redundant_dummies }, colors, &merged_color_names, optimize_peak_ram, sample_distance, n_threads)
}

#[cfg(test)]
mod tests {
    use crate::{coloring_interface::{ColorSetStorage, ColorSetView}, set_operations::test_utils::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING}};

    /// Merges the colorings of input_seqs_1 and input_seqs_2 (one color per sequence, all colors
    /// distinct) and checks the result against a coloring built directly from all sequences.
    fn check_merge(k: usize, input_seqs_1: &[Vec<u8>], input_seqs_2: &[Vec<u8>], input_sample_distance: usize, merge_sample_distance: usize, n_threads: usize) {
        let names: Vec<String> = (0..input_seqs_1.len() + input_seqs_2.len()).map(|i| format!("s{}", i)).collect();
        let (names1, names2) = names.split_at(input_seqs_1.len());
        let colors1: Vec<(&str, Vec<u8>)> = names1.iter().zip(input_seqs_1).map(|(name, seq)| (name.as_str(), seq.clone())).collect();
        let colors2: Vec<(&str, Vec<u8>)> = names2.iter().zip(input_seqs_2).map(|(name, seq)| (name.as_str(), seq.clone())).collect();

        // With distinct names, merging shared colors must not change anything
        for merge_shared_colors in [false, true] {
            check_named_merge(k, &colors1, &colors2, merge_shared_colors, input_sample_distance, merge_sample_distance, n_threads);
        }
    }

    /// Merges two colorings given as (color name, sequence) pairs (see assign_color_ids) and checks
    /// that every k-mer gets the same color set as in a coloring built directly from the expected
    /// merged colors. Colors with the same name are shared if merge_shared_colors is true.
    fn check_named_merge(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], merge_shared_colors: bool, input_sample_distance: usize, merge_sample_distance: usize, n_threads: usize) {
        let mut names1 = Vec::<String>::new();
        let mut names2 = Vec::<String>::new();
        let input_seqs_1 = assign_color_ids(&mut names1, colors1);
        let input_seqs_2 = assign_color_ids(&mut names2, colors2);

        // Expected merged colors: the colors of the first coloring keep their ids, and the new
        // colors of the second coloring come after them.
        let mut expected_names = names1.clone();
        let expected_seqs_2 = if merge_shared_colors {
            // A color of the second coloring with a name of the first coloring gets that id
            assign_color_ids(&mut expected_names, colors2)
        } else {
            // Every color of the second coloring is new, even if its name is not
            expected_names.extend(names2.iter().cloned());
            input_seqs_2.iter().map(|(color, seq)| (names1.len() + color, seq.clone())).collect()
        };
        let expected_seqs = [input_seqs_1.clone(), expected_seqs_2].concat();

        let ccc_both = build_named_coloring(k, &expected_names, &expected_seqs, input_sample_distance, n_threads);

        for keep_redundant_dummies in [false, true] {
            let ccc1 = build_named_coloring(k, &names1, &input_seqs_1, input_sample_distance, n_threads);
            let ccc2 = build_named_coloring(k, &names2, &input_seqs_2, input_sample_distance, n_threads);

            let ccc_merged = super::merge_compact_colex_kmers(ccc1, ccc2, merge_shared_colors, true, keep_redundant_dummies, merge_sample_distance, n_threads);
            assert_eq!(ccc_merged.get_color_names(), &expected_names);
            assert_eq!(ccc_merged.get_set_storage().n_colors(), expected_names.len());
            let sbwt_merged = &ccc_merged.sbwt();

            assert_eq!(sbwt_merged.n_kmers(), ccc_both.sbwt().n_kmers());
            if keep_redundant_dummies {
                assert!(sbwt_merged.n_sets() >= ccc_both.sbwt().n_sets());
            } else {
                // With the redundant dummies removed, we get the same nodes as when building from scratch
                assert_eq!(sbwt_merged.n_sets(), ccc_both.sbwt().n_sets());
            }

            for colex in 0..ccc_both.sbwt().n_sets() {
                let kmer = ccc_both.sbwt().access_kmer(colex);

                if kmer.iter().all(|c| *c != b'$') { // Not a dummy k-mer
                    let true_colors: Vec<usize> = ccc_both.colex_to_set(colex).iter().collect();
                    let range = sbwt_merged.search(&kmer).unwrap();
                    assert_eq!(range.len(), 1);
                    let colex_merged = range.start;
                    let merged_colors: Vec<usize> = ccc_merged.colex_to_set(colex_merged).iter().collect();

                    assert_eq!(true_colors, merged_colors, "k = {}, k-mer {}, keep_redundant_dummies = {}", k, String::from_utf8_lossy(&kmer), keep_redundant_dummies);
                }
            }
        }
    }

    // =====================================================================
    // Tests
    // =====================================================================

    #[test]
    fn test_merge() {

        // Opt-in logging: set RUST_LOG=info to see output; silent by default.
        let _ = env_logger::try_init();

        let n_threads = 3;

        for k in 3_usize..10_usize { // k < 3 does not work because construction uses 3-mer binning.

            let input_seqs_1: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (i + k.pow(4)) as u64)).collect();
            let input_seqs_2: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (123456 + i + k.pow(4)) as u64)).collect();

            check_merge(k, &input_seqs_1, &input_seqs_2, 3, 5, n_threads);
        }
    }

    #[test]
    fn test_merge_multiple_colored_subunitigs_in_one_unitig() {
        let _ = env_logger::try_init();
        let k = 11;
        let genome = gen_random_dna_string(80, 1);

        // Color 1 covers the middle of the unitig of color 0, so the single unitig
        // of the first coloring breaks into three colored subunitigs.
        let input_seqs_1 = vec![genome.clone(), genome[20..50].to_vec()];
        let input_seqs_2 = vec![gen_random_dna_string(80, 2)];

        check_merge(k, &input_seqs_1, &input_seqs_2, NO_SAMPLING, NO_SAMPLING, 3);
        check_merge(k, &input_seqs_2, &input_seqs_1, NO_SAMPLING, NO_SAMPLING, 3);
    }

    #[test]
    fn test_merge_shared_kmers_in_middle_of_unitig() {
        let _ = env_logger::try_init();
        let k = 11;
        let genome = gen_random_dna_string(80, 3);

        // The k-mers of the first coloring form a run in the middle of the unitig
        // of the second coloring, so processing the second coloring enters and
        // leaves a run of previously visited k-mers.
        let inner = vec![genome[20..50].to_vec()];
        let outer = vec![genome.clone()];
        check_merge(k, &inner, &outer, NO_SAMPLING, NO_SAMPLING, 3);
        check_merge(k, &outer, &inner, NO_SAMPLING, NO_SAMPLING, 3);

        // Partial overlap: a shared run at the start of one unitig and at the end of the other.
        let left = vec![genome[0..50].to_vec()];
        let right = vec![genome[30..80].to_vec()];
        check_merge(k, &left, &right, NO_SAMPLING, NO_SAMPLING, 3);
        check_merge(k, &right, &left, NO_SAMPLING, NO_SAMPLING, 3);
    }

    #[test]
    fn test_merge_overlapping_random_substrings() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};

        let _ = env_logger::try_init();
        let n_threads = 3;

        // Both colorings consist of random substrings of the same genome, so the
        // colorings share many k-mers and unitigs of the colorings partially overlap
        // in all kinds of ways.
        for k in [5_usize, 8, 11, 15] {
            for (sample_distance, seed) in [(NO_SAMPLING, 0_u64), (3, 1)] {
                let genome = gen_random_dna_string(300, 1000 + k as u64 + seed);
                let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(2000 + k as u64 * 100 + seed);
                let mut random_substrings = |n: usize| -> Vec<Vec<u8>> {
                    (0..n).map(|_| {
                        let start = rng.next_u64() as usize % (genome.len() - k);
                        let len = k + rng.next_u64() as usize % 80;
                        genome[start..(start + len).min(genome.len())].to_vec()
                    }).collect()
                };
                let input_seqs_1 = random_substrings(6);
                let input_seqs_2 = random_substrings(6);
                check_merge(k, &input_seqs_1, &input_seqs_2, sample_distance, sample_distance, n_threads);
            }
        }
    }

    #[test]
    fn test_merge_shared_colors() {
        let _ = env_logger::try_init();
        let n_threads = 3;

        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 5000 + k as u64);
                let seed = 6000 + 10 * k as u64;

                // A and D are private to one index. B has unrelated sequences in the two indexes.
                // C has overlapping substrings of the same genome in both indexes, so some k-mers
                // have C on both sides, and it must be reported only once.
                let colors1 = vec![
                    ("A", gen_random_dna_string(60, seed)),
                    ("B", gen_random_dna_string(60, seed + 1)),
                    ("B", genome[150..200].to_vec()),
                    ("C", genome[0..100].to_vec()),
                ];
                let colors2 = vec![
                    ("B", gen_random_dna_string(60, seed + 2)),
                    ("D", genome[40..160].to_vec()),
                    ("C", genome[50..130].to_vec()),
                    ("C", genome[140..190].to_vec()),
                ];

                check_named_merge(k, &colors1, &colors2, true, sample_distance, sample_distance, n_threads);
                check_named_merge(k, &colors2, &colors1, true, sample_distance, sample_distance, n_threads);

                // Without the flag, equally named colors stay distinct
                check_named_merge(k, &colors1, &colors2, false, sample_distance, sample_distance, n_threads);
            }
        }
    }

    #[test]
    fn test_merge_all_colors_shared() {
        let _ = env_logger::try_init();
        let k = 9;
        let genome = gen_random_dna_string(200, 77);
        let colors1 = vec![("X", genome[0..120].to_vec()), ("Y", genome[60..200].to_vec())];
        let colors2 = vec![("Y", genome[0..90].to_vec()), ("X", genome[100..200].to_vec())];
        check_named_merge(k, &colors1, &colors2, true, NO_SAMPLING, NO_SAMPLING, 3);
        check_named_merge(k, &colors1, &colors2, true, 3, 3, 3);
    }


}
