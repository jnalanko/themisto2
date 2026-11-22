//! Subsets of colored indexes: the k-mers and colors of a subset of the colors (e.g. of the species).

use std::sync::Arc;

use rand_chacha::rand_core::{RngCore, SeedableRng};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use sbwt::{dbg::{Dbg, Node}, LcsArray, SbwtIndex, SubsetMatrix};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{Coloring, DbgMarks}, coloring_interface::{ColorSetStorage, ColorSetView}, unitig_export::break_to_colored_subunitigs};
use crate::set_operations::{color_sets::ColorCombination, key_kmers::SubunitigBoundaries, pipeline::{colored_result, FirstInputSubunitigs}};

type Bits = bitvec::vec::BitVec<u64, bitvec::order::Lsb0>;

const KEPT_COLORS_FINGERPRINT_SEED: u64 = 2741;

/// The subset of a colored index for the colors marked in `kept_colors` (one entry per color): the
/// k-mers whose color set contains at least one of these colors, each with the kept colors of its
/// color set. The result has only the kept colors, with new ids in their original order: the same
/// index as one built from scratch from the sequences of the kept colors.
///
/// With `keep_all_colors`, the k-mers are the same, but each keeps its whole color set, and the
/// result has all colors of the index with their ids: the same k-mers and color sets as the
/// intersection of the index with the subset, with the union of the colors and shared colors
/// merged (which has more key k-mers than needed).
///
/// The index is given as its SBWT `sbwt`, which needs select support, the marks of its DBG and its
/// coloring. Returns these for the result, and the LCS array of its SBWT.
#[allow(clippy::too_many_arguments)]
pub fn subset_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(sbwt: SbwtIndex<SubsetMatrix>, marks: DbgMarks, coloring: &Coloring<CSS>, kept_colors: &[bool], keep_all_colors: bool, optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> (SbwtIndex<SubsetMatrix>, DbgMarks, Coloring<CSS>, LcsArray) {
    assert_eq!(kept_colors.len(), coloring.get_color_names().len());
    let mut result_color_names = Vec::<String>::new();
    let color_to_result: Vec<Option<usize>> = kept_colors.iter().zip(coloring.get_color_names()).map(|(&kept, name)| {
        (kept || keep_all_colors).then(|| {
            result_color_names.push(name.clone());
            result_color_names.len() - 1
        })
    }).collect();

    // The DBG takes over the marks, so they are freed with it once the kept k-mers are marked
    let dbg = Dbg::from_marks(&sbwt, marks.0, marks.1);
    log::info!("Marking the k-mers with the kept colors");
    let (keep, boundaries) = kmers_with_colors(&dbg, coloring, kept_colors, keep_all_colors, n_threads);
    log::info!("Keeping {} of {} k-mers", keep.count_ones(), sbwt.n_kmers());
    let group_starts = to_u64_words(dbg.k_minus_1_marks());
    drop(dbg);

    log::info!("Filtering SBWT");
    let precalc_len = sbwt.get_lookup_table().prefix_length;
    // The SBWT goes into an Arc for sbwt::filter, as in pipeline::set_operation
    let sbwt = Arc::new(sbwt);
    let (result_sbwt, result_merge_plan, removed) = sbwt::filter_with_removed_positions(sbwt.clone(), keep, Some(group_starts), precalc_len, optimize_peak_ram, n_threads)
        .unwrap_or_else(|_| panic!("The SBWT of the index has no select support"));
    // Of the merge plan, only the positions of the input (s1) are needed while dummy nodes have no
    // colors, and the starts of its (k-1)-suffix groups, for the groups of the result
    let sbwt::ThreeWayInterleaving { s1, s2, is_leader, .. } = result_merge_plan;
    // sbwt::filter has dropped its clone by now (see the unwraps in pipeline::set_operation)
    let sbwt = Arc::try_unwrap(sbwt).unwrap();

    // With all colors, they keep their ids
    let colors = ColorCombination::FromFirst { color1_to_result: (!keep_all_colors).then_some(&color_to_result[..]) };
    // The unitigs of the index were walked above, so the key k-mers come from their boundaries
    colored_result(FirstInputSubunitigs::Known(&sbwt, boundaries), coloring, None, result_sbwt, &s1, &s2, is_leader, Some(&removed), "subset", colors, &result_color_names, sample_distance, false, n_threads)
}

/// A bit vector with usize words, as sbwt computes them, with u64 words.
fn to_u64_words(bits: &bitvec::vec::BitVec) -> Bits {
    let mut result = Bits::from_vec(bits.as_raw_slice().iter().map(|&w| w as u64).collect());
    result.truncate(bits.len());
    result
}

/// Marks the k-mers of `coloring` whose color set contains at least one of the colors marked in
/// `kept_colors` (one entry per color): one bit per node of the SBWT, set for these k-mers. Dummy
/// nodes are not marked. `dbg` is the DBG of the SBWT of `coloring`.
///
/// A colored subunitig (see break_to_colored_subunitigs) has one color set, so the condition is
/// evaluated once per distinct color set, and the k-mers are marked one colored subunitig at a time.
/// Also returns the boundaries of the kept colored subunitigs, found in the same walk. Neighboring
/// kept subunitigs of a unitig with the same kept colors have the same color set in the result, so
/// they count as one there: the boundary between them is left out, and needs no key k-mers. With
/// `keep_all_colors`, the result keeps the whole color sets, so this holds only for the same color
/// set.
pub(crate) fn kmers_with_colors<CSS: ColorSetStorage + Send + Sync>(dbg: &Dbg<'_, SubsetMatrix>, coloring: &Coloring<CSS>, kept_colors: &[bool], keep_all_colors: bool, n_threads: usize) -> (Bits, SubunitigBoundaries) {
    let sets = coloring.get_set_storage();
    assert_eq!(kept_colors.len(), sets.n_colors());
    // The colors of each color set that are in the result, as a 128-bit fingerprint: the XOR of a
    // random value per such color, as set_of_sets_construction identifies the color sets of the
    // result
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(KEPT_COLORS_FINGERPRINT_SEED);
    let color_fingerprints: Vec<u128> = kept_colors.iter().map(|&kept| {
        let fingerprint = (rng.next_u64() as u128) << 64 | rng.next_u64() as u128;
        if kept || keep_all_colors { fingerprint } else { 0 }
    }).collect();
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let (set_is_kept, kept_fingerprints): (Vec<bool>, Vec<u128>) = thread_pool.install(|| {
        (0..sets.n_sets()).into_par_iter().map(|id| {
            let (mut is_kept, mut fingerprint) = (false, 0_u128);
            for color in sets.get_set_view(id).iter() {
                is_kept |= kept_colors[color];
                fingerprint ^= color_fingerprints[color];
            }
            (is_kept, fingerprint)
        }).unzip()
    });
    drop(color_fingerprints);

    let sbwt = dbg.sbwt();
    let marks = AtomicBitmap::new(sbwt.n_sets());
    let firsts = AtomicBitmap::new(sbwt.n_sets());
    let lasts = AtomicBitmap::new(sbwt.n_sets());
    // The kept k-mers at the ends of the unitigs, whose boundaries are decided once all kept
    // k-mers are known (below)
    let kept_unitig_firsts = AtomicBitmap::new(sbwt.n_sets());
    let kept_unitig_lasts = AtomicBitmap::new(sbwt.n_sets());
    dbg.iter_unitigs_with_callback(|nodes| {
        let unitig_colex_ranks: Vec<usize> = nodes.iter().map(|v| v.id).collect();
        let (set_ids, ranges) = break_to_colored_subunitigs(&unitig_colex_ranks, &[], coloring.get_map(), sbwt);
        // Whether the subunitig at j is kept with the same kept colors as the one at i
        let same_kept_colors = |i: usize, j: usize| set_is_kept[set_ids[j]] && kept_fingerprints[set_ids[j]] == kept_fingerprints[set_ids[i]];
        for (i, range) in ranges.iter().enumerate() {
            // A subunitig has one color set, so it is kept or removed as a whole. The k-mers of a
            // removed one are not in the result, so its boundaries would mark nothing there.
            if set_is_kept[set_ids[i]] {
                let (first, last) = (unitig_colex_ranks[range.start], unitig_colex_ranks[range.end - 1]);
                if i == 0 {
                    kept_unitig_firsts.set(first, true);
                } else if !same_kept_colors(i, i - 1) {
                    firsts.set(first, true);
                }
                if i + 1 == ranges.len() {
                    kept_unitig_lasts.set(last, true);
                } else if !same_kept_colors(i, i + 1) {
                    lasts.set(last, true);
                }
                for &colex in &unitig_colex_ranks[range.clone()] {
                    marks.set(colex, true);
                }
            }
        }
    }, n_threads);

    // At an end of a unitig, the k-mers that made the branch there may all be gone: if the kept
    // k-mer v has one kept neighbor u on that side, and u has one kept neighbor, v, the two are
    // inside a unitig of the result, and with the same kept colors that is no boundary there.
    // Otherwise the DBG of the result has its own boundary there (a branch or the end of a unitig).
    // Most ends fail the first test, which needs only the marks; the kept colors (a walk to the
    // next sampled k-mer for the first k-mer of a unitig) are looked up only for the others.
    let kept_colors_of = |colex: usize| kept_fingerprints[coloring.colex_to_set_id(colex, sbwt)];
    let single_kept_neighbor = |v: usize, out: bool, buf: &mut Vec<(Node, u8)>| -> Option<usize> {
        buf.clear();
        if out { dbg.push_out_neighbors(Node { id: v }, buf) } else { dbg.push_in_neighbors(Node { id: v }, buf) }
        let mut kept = buf.iter().map(|(u, _)| u.id).filter(|&u| marks.get(u));
        let first = kept.next()?;
        kept.next().is_none().then_some(first)
    };
    let joined_in_result = |v: usize, out: bool, buf: &mut Vec<(Node, u8)>| {
        single_kept_neighbor(v, out, buf).is_some_and(|u| {
            single_kept_neighbor(u, !out, buf) == Some(v) && kept_colors_of(u) == kept_colors_of(v)
        })
    };
    for (ends, boundaries, out) in [(kept_unitig_firsts, &firsts, false), (kept_unitig_lasts, &lasts, true)] {
        let ends = ends.into_bitvec();
        thread_pool.install(|| {
            ends.as_raw_slice().par_iter().enumerate().for_each_init(Vec::<(Node, u8)>::new, |buf, (word_index, &word)| {
                let mut bits = word;
                while bits != 0 {
                    let v = word_index * usize::BITS as usize + bits.trailing_zeros() as usize;
                    if !joined_in_result(v, out, buf) {
                        boundaries.set(v, true);
                    }
                    bits &= bits - 1;
                }
            });
        });
    }

    // Converted one after the other, so that only one u64 copy exists at a time (both bitmaps
    // exist from the walk until each is converted)
    let firsts = SubunitigBoundaries::sparse(&to_u64_words(&firsts.into_bitvec()));
    let lasts = SubunitigBoundaries::sparse(&to_u64_words(&lasts.into_bitvec()));
    let boundaries = SubunitigBoundaries { firsts, lasts };
    (to_u64_words(&marks.into_bitvec()), boundaries)
}

#[cfg(test)]
mod tests {
    use crate::{coloring_interface::ColorSetView, set_operations::test_utils::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING}};

    /// Checks kmers_with_colors on a coloring of (color name, sequence) pairs (see
    /// assign_color_ids) against the color set of every k-mer, for the colors named `kept`.
    fn check_kmers_with_colors(k: usize, colors: &[(&str, Vec<u8>)], kept: &[&str], sample_distance: usize) {
        let mut names = Vec::<String>::new();
        let seqs = assign_color_ids(&mut names, colors);
        let (sbwt, (sufgroups, dummy_marks), coloring) = build_named_coloring(k, &names, &seqs, sample_distance, 3);
        let kept_colors: Vec<bool> = names.iter().map(|name| kept.contains(&name.as_str())).collect();

        let dbg = sbwt::dbg::Dbg::from_marks(&sbwt, sufgroups, dummy_marks);
        let (marks, boundaries) = super::kmers_with_colors(&dbg, &coloring, &kept_colors, false, 3);
        // Every run of kept colored subunitigs with the same kept colors has a first and a last
        // k-mer, which are kept
        use simple_sds_sbwt::ops::{BitVec, Select};
        assert_eq!(boundaries.firsts.count_ones(), boundaries.lasts.count_ones());
        assert!(boundaries.firsts.one_iter().chain(boundaries.lasts.one_iter()).all(|(_, colex)| marks[colex]));
        assert_eq!(marks.len(), sbwt.n_sets());
        for colex in 0..sbwt.n_sets() {
            let kmer = sbwt.access_kmer(colex);
            let expected = !kmer.contains(&b'$') && coloring.colex_to_set(colex, &sbwt).iter().any(|color| kept_colors[color]);
            assert_eq!(marks[colex], expected, "k = {}, k-mer {}, kept {:?}, sample distance {}", k, String::from_utf8_lossy(&kmer), kept, sample_distance);
        }
    }

    /// Takes the subset of a coloring of (color name, sequence) pairs (see assign_color_ids) for the
    /// colors named `kept`, and checks that it is the same as a coloring built from scratch from the
    /// sequences of those colors: the same nodes, and the same color set for every k-mer.
    fn check_subset(k: usize, colors: &[(&str, Vec<u8>)], kept: &[&str], input_sample_distance: usize, sample_distance: usize) {
        let mut names = Vec::<String>::new();
        let seqs = assign_color_ids(&mut names, colors);
        let (sbwt, marks, coloring) = build_named_coloring(k, &names, &seqs, input_sample_distance, 3);
        let kept_colors: Vec<bool> = names.iter().map(|name| kept.contains(&name.as_str())).collect();
        let (result_sbwt, _, result_coloring, _) = super::subset_compact_colex_kmers(sbwt, marks, &coloring, &kept_colors, false, true, sample_distance, 3);
        let case = format!("k = {k}, kept {kept:?}, sample distances {input_sample_distance} and {sample_distance}");

        let mut expected_names = Vec::<String>::new();
        let expected_colors: Vec<(&str, Vec<u8>)> = colors.iter().filter(|(name, _)| kept.contains(name)).cloned().collect();
        let expected_seqs = assign_color_ids(&mut expected_names, &expected_colors);
        assert_eq!(result_coloring.get_color_names(), &expected_names, "{case}");
        if expected_seqs.is_empty() {
            assert_eq!(result_sbwt.n_kmers(), 0, "{case}");
            return;
        }
        let (expected_sbwt, _, expected_coloring) = build_named_coloring(k, &expected_names, &expected_seqs, input_sample_distance, 3);
        assert_eq!(result_sbwt.n_kmers(), expected_sbwt.n_kmers(), "{case}");
        assert_eq!(result_sbwt.n_sets(), expected_sbwt.n_sets(), "{case}");
        for colex in 0..expected_sbwt.n_sets() {
            let kmer = expected_sbwt.access_kmer(colex);
            if kmer.contains(&b'$') { continue; }
            let range = result_sbwt.search(&kmer).unwrap();
            assert_eq!(range.len(), 1, "{case}");
            let result_set: Vec<usize> = result_coloring.colex_to_set(range.start, &result_sbwt).iter().collect();
            let expected_set: Vec<usize> = expected_coloring.colex_to_set(colex, &expected_sbwt).iter().collect();
            assert_eq!(result_set, expected_set, "{case}, k-mer {}", String::from_utf8_lossy(&kmer));
        }
    }

    /// Takes the subset of a coloring of (color name, sequence) pairs (see assign_color_ids) for the
    /// colors named `kept` keeping all colors, and checks that it is the same as the intersection of
    /// the coloring with its subset, with the union of the colors and shared colors merged: the same
    /// colors, nodes and color set for every k-mer; and, without sampling, exactly the key k-mers
    /// that the color sets need.
    fn check_subset_keeping_all_colors(k: usize, colors: &[(&str, Vec<u8>)], kept: &[&str], sample_distance: usize) {
        let mut names = Vec::<String>::new();
        let seqs = assign_color_ids(&mut names, colors);
        let kept_colors: Vec<bool> = names.iter().map(|name| kept.contains(&name.as_str())).collect();
        let build = || build_named_coloring(k, &names, &seqs, sample_distance, 3);
        let (sbwt, marks, coloring) = build();
        let (result_sbwt, (result_sufgroups, result_dummy_marks), result_coloring, _) = super::subset_compact_colex_kmers(sbwt, marks, &coloring, &kept_colors, true, true, sample_distance, 3);
        let case = format!("k = {k}, kept {kept:?}, sample distance {sample_distance}");
        assert_eq!(result_coloring.get_color_names(), &names, "{case}");
        if !kept_colors.contains(&true) {
            assert_eq!(result_sbwt.n_kmers(), 0, "{case}");
            return;
        }

        let (sbwt, marks, coloring) = build();
        let (subset_sbwt, subset_marks, subset_coloring, _) = super::subset_compact_colex_kmers(sbwt, marks, &coloring, &kept_colors, false, true, sample_distance, 3);
        let (sbwt, marks, coloring) = build();
        let (expected_sbwt, _, expected_coloring, _) = crate::set_operations::intersect::intersect_compact_colex_kmers(sbwt, marks, &coloring, subset_sbwt, subset_marks, &subset_coloring, crate::set_operations::intersect::IntersectColors::Union, true, true, sample_distance, 3);
        assert_eq!(expected_coloring.get_color_names(), &names, "{case}");
        assert_eq!(result_sbwt.n_kmers(), expected_sbwt.n_kmers(), "{case}");
        assert_eq!(result_sbwt.n_sets(), expected_sbwt.n_sets(), "{case}");
        for colex in 0..expected_sbwt.n_sets() {
            let kmer = expected_sbwt.access_kmer(colex);
            if kmer.contains(&b'$') { continue; }
            let range = result_sbwt.search(&kmer).unwrap();
            assert_eq!(range.len(), 1, "{case}");
            let result_set: Vec<usize> = result_coloring.colex_to_set(range.start, &result_sbwt).iter().collect();
            let expected_set: Vec<usize> = expected_coloring.colex_to_set(colex, &expected_sbwt).iter().collect();
            assert_eq!(result_set, expected_set, "{case}, k-mer {}", String::from_utf8_lossy(&kmer));
        }
        // Without sampling, the key k-mers are exactly the last k-mers of the runs of k-mers with one
        // color set along the unitigs of the result (the in-neighbor of the first k-mer of a run is
        // the last of the run before it, or the last k-mer of another unitig)
        if sample_distance == NO_SAMPLING {
            use simple_sds_sbwt::ops::BitVec;
            let dbg = sbwt::dbg::Dbg::from_marks(&result_sbwt, result_sufgroups, result_dummy_marks);
            let key = crate::atomic_bitmap::AtomicBitmap::new(result_sbwt.n_sets());
            dbg.iter_unitigs_with_callback(|nodes| {
                let sets: Vec<Vec<usize>> = nodes.iter().map(|v| result_coloring.colex_to_set(v.id, &result_sbwt).iter().collect()).collect();
                for i in 0..nodes.len() {
                    if i + 1 == nodes.len() || sets[i] != sets[i + 1] { key.set(nodes[i].id, true); }
                }
            }, 1);
            let minimal = key.into_bitvec().count_ones();
            assert_eq!(result_coloring.get_map().sampling.count_ones(), minimal, "{case}");
        }
    }

    #[test]
    fn test_subset_keeping_all_colors() {
        let _ = env_logger::try_init();
        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 1300 + k as u64);
                // As in test_subset: the kept k-mers keep color sets with removed colors, which
                // change within unitigs; the branching case makes the subset cut unitigs
                let colors = vec![
                    ("A", genome[0..120].to_vec()),
                    ("B", genome[60..200].to_vec()),
                    ("C", genome[90..150].to_vec()),
                    ("D", gen_random_dna_string(80, 1400 + k as u64)),
                    ("B", gen_random_dna_string(50, 1450 + k as u64)),
                ];
                let mut b = genome[40..80].to_vec();
                b.extend(gen_random_dna_string(30, 1500 + k as u64));
                let mut c = gen_random_dna_string(30, 1600 + k as u64);
                c.extend(&genome[100..150]);
                let branching = vec![("A", genome[0..160].to_vec()), ("B", b), ("C", c)];
                for kept in [&["A"][..], &["B", "D"], &["C"], &["A", "C"], &["A", "B", "C", "D"], &[]] {
                    check_subset_keeping_all_colors(k, &colors, kept, sample_distance);
                }
                for kept in [&["A"][..], &["B"], &["C"], &["B", "C"]] {
                    check_subset_keeping_all_colors(k, &branching, kept, sample_distance);
                }
            }
        }
    }

    #[test]
    fn test_subset() {
        let _ = env_logger::try_init();
        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 800 + k as u64);
                // Overlapping substrings of one genome, so that the subset cuts unitigs and
                // colored subunitigs, and a color with two sequences
                let colors = vec![
                    ("A", genome[0..120].to_vec()),
                    ("B", genome[60..200].to_vec()),
                    ("C", genome[90..150].to_vec()),
                    ("D", gen_random_dna_string(80, 900 + k as u64)),
                    ("B", gen_random_dna_string(50, 950 + k as u64)),
                ];
                for kept in [&["A"][..], &["B", "D"], &["C"], &["A", "C"], &["A", "B", "C", "D"], &[]] {
                    check_subset(k, &colors, kept, sample_distance, sample_distance);
                }
            }
        }
    }

    /// Dropping colors can give neighboring colored subunitigs of the index the same kept colors,
    /// and remove the branches between its unitigs. These boundaries are no boundaries in the
    /// subset, so they must not leave key k-mers there: with one kept color of one sequence, the
    /// subset has exactly the key k-mers of an index built from it.
    #[test]
    fn test_subset_merges_subunitigs() {
        use simple_sds_sbwt::ops::BitVec;
        for k in [8_usize, 11] {
            let genome = gen_random_dna_string(200, 1000 + k as u64);
            // Within A, the color sets change where B and C start and end
            let overlapping = vec![
                ("A", genome[0..120].to_vec()),
                ("B", genome[60..200].to_vec()),
                ("C", genome[90..150].to_vec()),
            ];
            // B leaves A and C joins it, so that A spans several unitigs
            let mut b = genome[40..80].to_vec();
            b.extend(gen_random_dna_string(30, 1100 + k as u64));
            let mut c = gen_random_dna_string(30, 1200 + k as u64);
            c.extend(&genome[100..150]);
            let branching = vec![("A", genome[0..160].to_vec()), ("B", b), ("C", c)];
            for (case, colors) in [("overlapping", &overlapping), ("branching", &branching)] {
                for sample_distance in [NO_SAMPLING, 3] {
                    let mut names = Vec::<String>::new();
                    let seqs = assign_color_ids(&mut names, colors);
                    let (sbwt, marks, coloring) = build_named_coloring(k, &names, &seqs, sample_distance, 3);
                    let kept_colors: Vec<bool> = names.iter().map(|name| name == "A").collect();
                    let (_, _, result_coloring, _) = super::subset_compact_colex_kmers(sbwt, marks, &coloring, &kept_colors, false, true, sample_distance, 3);

                    let mut expected_names = Vec::<String>::new();
                    let expected_seqs = assign_color_ids(&mut expected_names, &colors[..1]);
                    let (_, _, expected_coloring) = build_named_coloring(k, &expected_names, &expected_seqs, sample_distance, 3);
                    assert_eq!(result_coloring.get_map().sampling.count_ones(), expected_coloring.get_map().sampling.count_ones(), "k = {k}, sample distance {sample_distance}, {case}");
                    for kept in [&["A"][..], &["B"], &["C"], &["B", "C"], &["A", "B", "C"]] {
                        check_subset(k, colors, kept, sample_distance, sample_distance);
                    }
                }
            }
        }
    }

    #[test]
    fn test_kmers_with_colors() {
        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 600 + k as u64);
                // Overlapping substrings of one genome, so that colored subunitigs with different
                // color sets follow each other within unitigs
                let colors = vec![
                    ("A", genome[0..120].to_vec()),
                    ("B", genome[60..200].to_vec()),
                    ("C", genome[90..150].to_vec()),
                    ("D", gen_random_dna_string(80, 700 + k as u64)),
                ];
                for kept in [&["A"][..], &["B", "D"], &["C"], &["A", "B", "C", "D"], &[]] {
                    check_kmers_with_colors(k, &colors, kept, sample_distance);
                }
            }
        }
    }
}
