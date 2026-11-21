//! The steps that all set operations on colored indexes share: interleaving the two SBWTs,
//! computing the SBWT of the result from the interleaving, and building the colors of the result
//! from the colors of the inputs.

use std::{cmp::max, sync::Arc};

use sbwt::{dbg::Dbg, LcsArray, SbwtIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::ColorSetStorage, util::u64_bitvec_to_simple_sds};
use crate::set_operations::{color_sets::{ColorCombination, ElementGeneratorFromInterleaving}, key_kmers::{mark_key_kmers_for, mark_key_kmers_from_boundaries, mark_structural_key_kmers, SubunitigBoundaries}, result_positions::InputToResult};
use crate::set_of_sets_construction::{build_color_set_storage, find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates};

/// The k-mers of the result of a set operation.
#[derive(Clone, Copy, Debug)]
pub(crate) enum SbwtOperation {
    /// The k-mers of both inputs. If keep_redundant_dummies is true, the dummy nodes that become
    /// redundant in the merge are kept: the result has the same k-mers and colors, but possibly
    /// more dummy nodes.
    Merge { keep_redundant_dummies: bool },
    /// The k-mers that are in both inputs.
    Intersect,
    /// The k-mers of the first input that are not in the second input.
    Difference,
}

impl SbwtOperation {
    /// How the result is called in the log
    fn result_name(&self) -> &'static str {
        match self {
            SbwtOperation::Merge { .. } => "merged",
            SbwtOperation::Intersect => "intersection",
            SbwtOperation::Difference => "difference",
        }
    }
}

type Bits = bitvec::vec::BitVec<u64, bitvec::order::Lsb0>;

/// Computes the result of `operation` on the k-mers of coloring1 and coloring2. The color set of
/// each k-mer of the result combines its color sets in the inputs as given by `colors`, and the
/// result has the colors `result_color_names`.
///
/// For the intersection and the difference, the SBWT of coloring1 needs select support.
#[allow(clippy::too_many_arguments)]
pub(crate) fn set_operation<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, operation: SbwtOperation, colors: ColorCombination, result_color_names: &[String], optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let interleaving = sbwt::MergeInterleaving::new(coloring1.sbwt(), coloring2.sbwt(), optimize_peak_ram, n_threads);
    assert_eq!(interleaving.s1.len(), interleaving.s2.len());

    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);
    // Temporarily destructure the colorings into parts in order to be able to put the SBWTs into
    // Arcs for the set operations of sbwt, which take Arcs by design (see the comment at sbwt::merge).
    let (sbwt1, lcs1, map1, sets1, color_names_1) = coloring1.into_parts();
    let (sbwt2, lcs2, map2, sets2, color_names_2) = coloring2.into_parts();
    // No k-mer of a difference is in coloring2, and its colors are those of coloring1, so the
    // difference needs only the SBWT of coloring2, for the set operation
    let is_difference = matches!(operation, SbwtOperation::Difference);
    let colors2 = if is_difference {
        drop((lcs2, map2, sets2, color_names_2));
        None
    } else {
        Some((lcs2, map2, sets2, color_names_2))
    };
    let sbwt1 = Arc::new(sbwt1);
    let sbwt2 = Arc::new(sbwt2);
    let interleaving = Arc::new(interleaving);

    // The result consists of the positions of an interleaving of the inputs that were not removed,
    // in order. Of that interleaving, only the positions of the inputs (s1 and s2) are needed
    // while dummy nodes have no colors, and the starts of its (k-1)-suffix groups, for the groups
    // of the result. The clones here clone just the Arcs.
    let (result_sbwt, s1, s2, is_leader, removed): (SbwtIndex<SubsetMatrix>, Bits, Bits, Bits, Option<Bits>) = match operation {
        SbwtOperation::Merge { keep_redundant_dummies } => {
            log::info!("Merging SBWTs");
            let (result_sbwt, removed) = if keep_redundant_dummies {
                (sbwt::merge_without_cleanup(sbwt1.clone(), sbwt2.clone(), interleaving.clone(), precalc_len, n_threads), None)
            } else {
                let (result_sbwt, removed) = sbwt::merge_with_removed_positions(sbwt1.clone(), sbwt2.clone(), interleaving.clone(), precalc_len, n_threads);
                (result_sbwt, Some(removed))
            };
            let sbwt::MergeInterleaving { s1, s2, is_leader, .. } = Arc::try_unwrap(interleaving).unwrap();
            (result_sbwt, s1, s2, is_leader, removed)
        },
        SbwtOperation::Intersect => {
            log::info!("Intersecting SBWTs");
            // The interleaving is passed on without keeping a reference, so that sbwt can free it
            // or take it over.
            let (result_sbwt, result_interleaving, removed) = sbwt::intersect_with_removed_positions(sbwt1.clone(), sbwt2.clone(), interleaving, precalc_len, optimize_peak_ram, n_threads)
                .unwrap_or_else(|_| panic!("The SBWT of the first index has no select support"));
            let sbwt::ThreeWayInterleaving { s1, s2, is_leader, .. } = result_interleaving;
            (result_sbwt, s1, s2, is_leader, Some(removed))
        },
        SbwtOperation::Difference => {
            log::info!("Subtracting SBWTs");
            // As for the intersection
            let (result_sbwt, result_interleaving, removed) = sbwt::difference_with_removed_positions(sbwt1.clone(), sbwt2.clone(), interleaving, precalc_len, optimize_peak_ram, n_threads)
                .unwrap_or_else(|_| panic!("The SBWT of the first index has no select support"));
            let sbwt::ThreeWayInterleaving { s1, s2, is_leader, .. } = result_interleaving;
            (result_sbwt, s1, s2, is_leader, Some(removed))
        },
    };

    // Put the coloring structs back together
    let sbwt1 = Arc::try_unwrap(sbwt1).unwrap();
    let sbwt2 = Arc::try_unwrap(sbwt2).unwrap();
    let coloring1 = CompactColexKmers::<CSS>::new(sbwt1, lcs1, map1, sets1, Some(&color_names_1));
    let coloring2 = match colors2 {
        Some((lcs2, map2, sets2, color_names_2)) => Some(CompactColexKmers::<CSS>::new(sbwt2, lcs2, map2, sets2, Some(&color_names_2))),
        None => {
            drop(sbwt2);
            None
        },
    };
    // Every k-mer of a merge is in one of the inputs, so it is visited when marking for the inputs
    // (and the check searches the result, which needs its LCS array)
    let check_visited = matches!(operation, SbwtOperation::Merge { .. });
    colored_result(&coloring1, None, coloring2.as_ref(), result_sbwt, &s1, &s2, is_leader, removed.as_ref(), operation.result_name(), colors, result_color_names, sample_distance, check_visited, n_threads)
}

/// Builds the colored index of the result of a set operation from its SBWT, which consists of the
/// positions of a merge plan of the inputs that were not removed, in order: `s1` and `s2` are the
/// positions of coloring1 and coloring2, and `removed` the removed positions (None if nothing was
/// removed). The color sets of the k-mers of the result combine their color sets in the inputs as
/// given by `colors`, and the result has the colors `result_color_names`. Without coloring2, its
/// k-mers are not used for the key k-mers, and `colors` must be ColorCombination::FromFirst.
/// `is_leader` are the starts of the (k-1)-suffix groups of the merge plan, from which the DBG of
/// the result takes its groups (see sbwt::result_k_minus_1_marks). `boundaries1` are the colored subunitigs of coloring1 if they are already known, so that its
/// unitigs are not walked again for the key k-mers; they are freed once the key k-mers are marked.
/// `check_visited` (see mark_key_kmers) needs that walk, so it can't be combined with `boundaries1`.
/// `result_name` is how the result is called in the log.
#[allow(clippy::too_many_arguments)]
pub(crate) fn colored_result<CSS: ColorSetStorage + Send + Sync>(coloring1: &CompactColexKmers<CSS>, boundaries1: Option<SubunitigBoundaries>, coloring2: Option<&CompactColexKmers<CSS>>, mut result_sbwt: SbwtIndex<SubsetMatrix>, s1: &Bits, s2: &Bits, is_leader: Bits, removed: Option<&Bits>, result_name: &str, colors: ColorCombination, result_color_names: &[String], sample_distance: usize, check_visited: bool, n_threads: usize) -> CompactColexKmers<CSS> {
    // First, so that the merge plan's group starts are freed before anything else is built
    log::info!("Deriving the (k-1)-suffix groups of the {} SBWT from the merge plan", result_name);
    let k_minus_1_marks = sbwt::result_k_minus_1_marks(&is_leader, removed, n_threads);
    drop(is_leader);

    result_sbwt.build_select();

    log::info!("Building the LCS array for the {} SBWT", result_name);
    let result_lcs = LcsArray::from_sbwt(&result_sbwt, n_threads, true);

    log::info!("Initializing DBG for the {} SBWT", result_name);
    // Debug builds only (e.g. the tests): checks the groups from the merge plan against those from
    // the LCS array of the result. Release builds skip this.
    if cfg!(debug_assertions) {
        assert_eq!(&k_minus_1_marks, Dbg::new(&result_sbwt, Some(&result_lcs), n_threads).k_minus_1_marks(), "the groups of the {} from the merge plan", result_name);
    }
    let result_dbg = Dbg::new_with_k_minus_1_marks(&result_sbwt, k_minus_1_marks);

    log::info!("=== Phase 1/3: marking key k-mers ===");
    let key_kmer_marks = mark_key_kmers(coloring1, boundaries1.as_ref(), coloring2, s1, s2, removed, &result_sbwt, &result_lcs, &result_dbg, sample_distance, check_visited, n_threads);
    // Only the key k-mers need these
    drop((result_dbg, boundaries1));
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / result_sbwt.n_kmers() as f64 * 100.0);

    let n_colors = u32::try_from(result_color_names.len()).unwrap_or_else( |_| {
        log::error!("Maximum number of colors 2^32 exceeded");
        panic!();
    });
    let gen = || ElementGeneratorFromInterleaving {
        s1,
        s2,
        removed,
        coloring1,
        coloring2,
        result_key_kmer_marks: &key_kmer_marks,
        filter: None,
        colors,
        n_result_colors: n_colors as usize,
    };

    log::info!("=== PHASE 2/3: Building color set fingerprints for key k-mers ===");
    let random_seed = 123123; // Todo: be more random
    // TODO: do not clone key_kmer_marks
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id, _) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen(), key_kmer_marks.clone(), n_colors, n_threads, random_seed);

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    let css = build_color_set_storage(n_colors as usize, repr_kmer_marks, distinct_set_sizes, gen(), n_threads);

    log::info!("Building rank support for key k-mer marks");
    let mut key_kmer_marks = crate::util::bitvec_to_simple_sds_bitvec(key_kmer_marks);
    key_kmer_marks.enable_rank();
    assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));

    let colex_map = ColexToColorSetMap {
        sampling: key_kmer_marks,
        color_set_ids: key_kmer_idx_to_set_id,
    };

    CompactColexKmers::<CSS>::new(result_sbwt, result_lcs, colex_map, css, Some(result_color_names))
}

/// Marks the key k-mers of the result: the structural key k-mers of the result, and the k-mers of
/// the result where a colored subunitig of an input ends (see mark_key_kmers_for). In an
/// intersection, a color set changes between consecutive k-mers x and y of a unitig of the result
/// only if the color set of x or y changes in one of the inputs, which has the edge from x to y
/// too. Then x is the last k-mer of a colored subunitig of that input, so it is marked. In a
/// difference, the color sets are those of coloring1, which has the edge from x to y too, so
/// marking for coloring1 suffices, and coloring2 is None.
///
/// With `boundaries1`, the key k-mers for coloring1 are marked from its known colored subunitigs
/// instead of a walk over its unitigs.
///
/// With check_visited, checks in debug builds that every k-mer of the result is visited while
/// marking for the inputs (see mark_key_kmers_for). Skipped in release builds: it costs a full
/// extra streaming-index pass over the result and is not needed for marking key k-mers, only for
/// catching bugs during development.
#[allow(clippy::too_many_arguments)]
fn mark_key_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: &CompactColexKmers<CSS>, boundaries1: Option<&SubunitigBoundaries>, coloring2: Option<&CompactColexKmers<CSS>>, s1: &Bits, s2: &Bits, removed: Option<&Bits>, result_sbwt: &SbwtIndex<SubsetMatrix>, result_lcs: &LcsArray, result_dbg: &Dbg<'_, SubsetMatrix>, sample_distance: usize, check_visited: bool, n_threads: usize) -> bitvec::vec::BitVec {
    let k = result_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert!(coloring2.is_none_or(|coloring2| coloring2.get_k() == k));

    let key_kmer_marks = AtomicBitmap::new(result_sbwt.n_sets());
    mark_structural_key_kmers(result_sbwt, result_dbg, &key_kmer_marks, sample_distance, n_threads);

    let visited_marks = (check_visited && cfg!(debug_assertions)).then(|| AtomicBitmap::new(result_sbwt.n_sets()));
    // TODO: This is a second copy of the removed positions (n bits), next to the original that
    // colors 2 and 3 use. Converting them to simple-sds once, with rank support, and sharing that
    // with colors 2 and 3 would save the copy.
    let removed = removed.map(|r| {
        let mut r = u64_bitvec_to_simple_sds(r);
        r.enable_rank();
        r
    });
    if let Some(boundaries1) = boundaries1 {
        let to_result = InputToResult::new(s1, removed.as_ref());
        mark_key_kmers_from_boundaries(boundaries1, |colex| to_result.result_colex(colex), result_dbg, &key_kmer_marks, n_threads);
    }
    let walked_coloring1 = boundaries1.is_none().then_some((coloring1, s1));
    for (coloring, in_input) in walked_coloring1.into_iter().chain(coloring2.map(|coloring2| (coloring2, s2))) {
        let to_result = InputToResult::new(in_input, removed.as_ref());
        mark_key_kmers_for(coloring, |colex| to_result.result_colex(colex), result_sbwt, result_lcs, result_dbg, &key_kmer_marks, visited_marks.as_ref(), n_threads);
    }
    drop(removed);

    if let Some(visited_marks) = visited_marks {
        assert_eq!(visited_marks.into_bitvec().count_ones(), result_sbwt.n_kmers());
    }

    key_kmer_marks.into_bitvec()
}
