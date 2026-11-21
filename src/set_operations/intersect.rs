use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};
use crate::set_operations::{color_sets::ColorCombination, colors::{intersected_color_mapping, merged_color_mapping}, pipeline::{set_operation, SbwtOperation}};

/// How the color set of a k-mer of the intersection is computed from its color sets in the two
/// input indexes.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntersectColors {
    /// The union of the two color sets, as when merging the indexes
    Union,
    /// The intersection of the two color sets, matching colors by name. The result has only the
    /// colors whose names are in both indexes.
    Intersect,
}

/// Intersects two colored indexes: the result has the k-mers that are in both indexes. The color
/// set of each k-mer is the union or the intersection of its color sets in the two indexes,
/// depending on `colors`. With the union, merge_shared_colors is as in
/// merge::merge_compact_colex_kmers. With the intersection, colors are always matched by name, and
/// a k-mer whose color sets do not intersect has an empty color set.
///
/// Both SBWTs need select support.
#[allow(clippy::too_many_arguments)]
pub fn intersect_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, colors: IntersectColors, merge_shared_colors: bool, optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {
    let fail = |e: String| -> ! {
        log::error!("{}", e);
        panic!("{}", e);
    };
    let (names1, names2) = (coloring1.get_color_names(), coloring2.get_color_names());
    match colors {
        IntersectColors::Union => {
            let (color2_to_result, result_color_names) = merged_color_mapping(names1, names2, merge_shared_colors).unwrap_or_else(|e| fail(e));
            let colors = ColorCombination::Union { color2_to_result: &color2_to_result };
            set_operation(coloring1, coloring2, SbwtOperation::Intersect, colors, &result_color_names, optimize_peak_ram, sample_distance, n_threads)
        },
        IntersectColors::Intersect => {
            let (color1_to_result, color2_to_result, result_color_names) = intersected_color_mapping(names1, names2).unwrap_or_else(|e| fail(e));
            let colors = ColorCombination::Intersection { color1_to_result: &color1_to_result, color2_to_result: &color2_to_result };
            set_operation(coloring1, coloring2, SbwtOperation::Intersect, colors, &result_color_names, optimize_peak_ram, sample_distance, n_threads)
        },
    }
}

#[cfg(test)]
mod tests {
    use crate::{coloring_interface::{ColorSetStorage, ColorSetView}, set_operations::test_utils::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING}};

    use super::IntersectColors;

    /// Intersects two colorings given as (color name, sequence) pairs (see assign_color_ids) and
    /// checks that the result has exactly the k-mers of both colorings, and that each k-mer has the
    /// union or intersection of its color sets in the inputs, computed here from the color names.
    fn check_intersect(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], colors: IntersectColors, merge_shared_colors: bool, input_sample_distance: usize, sample_distance: usize) {
        let n_threads = 3;
        let mut names1 = Vec::<String>::new();
        let mut names2 = Vec::<String>::new();
        let seqs1 = assign_color_ids(&mut names1, colors1);
        let seqs2 = assign_color_ids(&mut names2, colors2);

        // Expected result color names, and the result color id of each color of both inputs
        let (expected_names, color1_to_result, color2_to_result): (Vec<String>, Vec<Option<usize>>, Vec<Option<usize>>) = match (colors, merge_shared_colors) {
            (IntersectColors::Union, false) => {
                let names = [names1.clone(), names2.clone()].concat();
                (names, (0..names1.len()).map(Some).collect(), (names1.len()..names1.len() + names2.len()).map(Some).collect())
            },
            (IntersectColors::Union, true) => {
                let mut names = names1.clone();
                for name in names2.iter() {
                    if !names.contains(name) { names.push(name.clone()); }
                }
                let id = |name: &String| names.iter().position(|x| x == name);
                let (map1, map2) = (names1.iter().map(id).collect(), names2.iter().map(id).collect());
                (names, map1, map2)
            },
            (IntersectColors::Intersect, _) => {
                let names: Vec<String> = names1.iter().filter(|name| names2.contains(name)).cloned().collect();
                let id = |name: &String| names.iter().position(|x| x == name);
                let (map1, map2) = (names1.iter().map(id).collect(), names2.iter().map(id).collect());
                (names, map1, map2)
            },
        };

        let ccc1 = build_named_coloring(k, &names1, &seqs1, input_sample_distance, n_threads);
        let ccc2 = build_named_coloring(k, &names2, &seqs2, input_sample_distance, n_threads);

        // The expected k-mers with their expected color sets, from searching the k-mers of ccc1 in ccc2
        let mut expected = Vec::<(Vec<u8>, Vec<usize>)>::new();
        for colex1 in 0..ccc1.sbwt().n_sets() {
            let kmer = ccc1.sbwt().access_kmer(colex1);
            if kmer.contains(&b'$') { continue; } // Dummy
            let Some(range2) = ccc2.sbwt().search(&kmer) else { continue; };
            if range2.is_empty() { continue; }
            assert_eq!(range2.len(), 1);
            let set1: Vec<usize> = ccc1.colex_to_set(colex1).iter().filter_map(|c| color1_to_result[c]).collect();
            let set2: Vec<usize> = ccc2.colex_to_set(range2.start).iter().filter_map(|c| color2_to_result[c]).collect();
            let mut set: Vec<usize> = match colors {
                IntersectColors::Union => [set1, set2].concat(),
                IntersectColors::Intersect => set1.into_iter().filter(|c| set2.contains(c)).collect(),
            };
            set.sort_unstable();
            set.dedup();
            expected.push((kmer, set));
        }

        let result = super::intersect_compact_colex_kmers(ccc1, ccc2, colors, merge_shared_colors, true, sample_distance, n_threads);
        assert_eq!(result.get_color_names(), &expected_names);
        assert_eq!(result.get_set_storage().n_colors(), expected_names.len());
        assert_eq!(result.sbwt().n_kmers(), expected.len());

        for (kmer, expected_set) in expected {
            let range = result.sbwt().search(&kmer).unwrap();
            assert_eq!(range.len(), 1);
            let mut set: Vec<usize> = result.colex_to_set(range.start).iter().collect();
            set.sort_unstable();
            assert_eq!(set, expected_set, "k = {}, k-mer {}, colors = {:?}, merge_shared_colors = {}", k, String::from_utf8_lossy(&kmer), colors, merge_shared_colors);
        }
    }

    fn check_all_modes(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], sample_distance: usize) {
        for (colors, merge_shared_colors) in [(IntersectColors::Union, false), (IntersectColors::Union, true), (IntersectColors::Intersect, false)] {
            check_intersect(k, colors1, colors2, colors, merge_shared_colors, sample_distance, sample_distance);
            check_intersect(k, colors2, colors1, colors, merge_shared_colors, sample_distance, sample_distance);
        }
    }

    #[test]
    fn test_intersect_overlapping_substrings() {
        let _ = env_logger::try_init();
        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 7000 + k as u64);
                let seed = 8000 + 10 * k as u64;

                // A and D are private to one index. B has unrelated sequences in the two indexes
                // apart from one shared stretch. C has overlapping substrings of the genome in both.
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
                    ("C", genome[140..190].to_vec()),
                ];
                check_all_modes(k, &colors1, &colors2, sample_distance);
            }
        }
    }

    #[test]
    fn test_intersect_random_substrings() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};
        let _ = env_logger::try_init();

        // Both indexes consist of random substrings of the same genome, with color names drawn
        // from a small pool, so the intersection cuts the unitigs and colored subunitigs of the
        // inputs in all kinds of ways.
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
                check_all_modes(k, &colors1, &colors2, sample_distance);
            }
        }
    }

    #[test]
    fn test_intersect_no_shared_colors() {
        let _ = env_logger::try_init();
        let k = 9;
        let genome = gen_random_dna_string(150, 11);
        let colors1 = vec![("X", genome[0..100].to_vec())];
        let colors2 = vec![("Y", genome[50..150].to_vec())];
        // With the intersection, every k-mer gets an empty color set
        check_all_modes(k, &colors1, &colors2, NO_SAMPLING);
    }


}
