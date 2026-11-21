use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};
use crate::set_operations::{color_sets::ColorCombination, pipeline::{set_operation, SbwtOperation}};

/// Subtracts the second colored index from the first: the result has the k-mers of the first index
/// that are not in the second. No k-mer of the result is in the second index, so each k-mer keeps
/// its color set of the first index, and the result has the colors of the first index, with the
/// same ids (also colors that no k-mer of the result has anymore).
///
/// The SBWT of the first index needs select support.
pub fn difference_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {
    let result_color_names = coloring1.get_color_names().to_vec();
    set_operation(coloring1, coloring2, SbwtOperation::Difference, ColorCombination::FromFirst { color1_to_result: None }, &result_color_names, optimize_peak_ram, sample_distance, n_threads)
}

#[cfg(test)]
mod tests {
    use crate::{coloring_interface::{ColorSetStorage, ColorSetView}, set_operations::test_utils::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING}};

    /// Subtracts two colorings given as (color name, sequence) pairs (see assign_color_ids) and
    /// checks that the result has exactly the k-mers of the first coloring that are not in the
    /// second, each with its color set in the first coloring.
    fn check_difference(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], input_sample_distance: usize, sample_distance: usize) {
        let n_threads = 3;
        let mut names1 = Vec::<String>::new();
        let mut names2 = Vec::<String>::new();
        let seqs1 = assign_color_ids(&mut names1, colors1);
        let seqs2 = assign_color_ids(&mut names2, colors2);

        let ccc1 = build_named_coloring(k, &names1, &seqs1, input_sample_distance, n_threads);
        let ccc2 = build_named_coloring(k, &names2, &seqs2, input_sample_distance, n_threads);

        // The expected k-mers with their expected color sets
        let mut expected = Vec::<(Vec<u8>, Vec<usize>)>::new();
        for colex1 in 0..ccc1.sbwt().n_sets() {
            let kmer = ccc1.sbwt().access_kmer(colex1);
            if kmer.contains(&b'$') { continue; } // Dummy
            if ccc2.sbwt().search(&kmer).is_some_and(|range| !range.is_empty()) { continue; }
            expected.push((kmer, ccc1.colex_to_set(colex1).iter().collect()));
        }

        let result = super::difference_compact_colex_kmers(ccc1, ccc2, true, sample_distance, n_threads);
        assert_eq!(result.get_color_names(), &names1);
        assert_eq!(result.get_set_storage().n_colors(), names1.len());
        assert_eq!(result.sbwt().n_kmers(), expected.len());

        for (kmer, expected_set) in expected {
            let range = result.sbwt().search(&kmer).unwrap();
            assert_eq!(range.len(), 1);
            let set: Vec<usize> = result.colex_to_set(range.start).iter().collect();
            assert_eq!(set, expected_set, "k = {}, k-mer {}", k, String::from_utf8_lossy(&kmer));
        }
    }

    fn check_both_ways(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], sample_distance: usize) {
        check_difference(k, colors1, colors2, sample_distance, sample_distance);
        check_difference(k, colors2, colors1, sample_distance, sample_distance);
    }

    #[test]
    fn test_difference_overlapping_substrings() {
        let _ = env_logger::try_init();
        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 7000 + k as u64);
                let seed = 8000 + 10 * k as u64;

                // Overlapping substrings of the same genome, so the difference cuts unitigs and
                // colored subunitigs of the first coloring in the middle
                let colors1 = vec![
                    ("A", genome[30..120].to_vec()),
                    ("B", gen_random_dna_string(60, seed)),
                    ("B", genome[150..200].to_vec()),
                    ("C", genome[0..100].to_vec()),
                ];
                let colors2 = vec![
                    ("B", gen_random_dna_string(60, seed + 1)),
                    ("B", genome[160..190].to_vec()),
                    ("D", genome[40..160].to_vec()),
                    ("C", genome[50..130].to_vec()),
                ];
                check_both_ways(k, &colors1, &colors2, sample_distance);
            }
        }
    }

    #[test]
    fn test_difference_random_substrings() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};
        let _ = env_logger::try_init();

        // Both indexes consist of random substrings of the same genome, so the difference cuts the
        // unitigs and colored subunitigs of the first index in all kinds of ways.
        for k in [5_usize, 8, 11, 15] {
            for (sample_distance, seed) in [(NO_SAMPLING, 0_u64), (3, 1)] {
                let genome = gen_random_dna_string(300, 9000 + k as u64 + seed);
                let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(9500 + k as u64 * 100 + seed);
                let names = ["P", "Q", "R", "S"];
                let mut random_colors = |n: usize, n_names: usize| -> Vec<(&str, Vec<u8>)> {
                    (0..n).map(|_| {
                        let start = rng.next_u64() as usize % (genome.len() - k);
                        let len = k + rng.next_u64() as usize % 120;
                        (names[rng.next_u64() as usize % n_names], genome[start..(start + len).min(genome.len())].to_vec())
                    }).collect()
                };
                let colors1 = random_colors(8, 3);
                let colors2 = random_colors(8, 4);
                check_both_ways(k, &colors1, &colors2, sample_distance);
            }
        }
    }

    #[test]
    fn test_difference_no_overlap_and_full_overlap() {
        let _ = env_logger::try_init();
        let k = 9;
        let genome = gen_random_dna_string(150, 11);
        let colors1 = vec![("X", genome[0..100].to_vec()), ("Y", genome[20..60].to_vec())];
        let other = vec![("Z", gen_random_dna_string(100, 12))];
        // Nothing in common: the result is the first index
        check_both_ways(k, &colors1, &other, NO_SAMPLING);
        // Everything in common: the result is empty
        check_difference(k, &colors1, &colors1, NO_SAMPLING, NO_SAMPLING);
    }

    #[test]
    fn test_difference_cycle() {
        // A difference whose k-mers form a cycle needs no dummy chains, which takes the direct
        // pass of sbwt::difference
        let _ = env_logger::try_init();
        let k = 5;
        let cycle: Vec<u8> = b"ACG".iter().cycle().take(k + 32).copied().collect();
        check_difference(k, &[("X", cycle)], &[("Y", b"TTTTTTTT".to_vec())], NO_SAMPLING, NO_SAMPLING);
    }
}
