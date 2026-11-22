//! The steps that all set operations on colored indexes share: interleaving the two SBWTs,
//! computing the SBWT of the result from the interleaving, and building the colors of the result
//! from the colors of the inputs.

use std::{cmp::max, sync::Arc};

use sbwt::{dbg::Dbg, LcsArray, SbwtIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, Coloring, DbgMarks}, coloring_interface::ColorSetStorage, util::u64_bitvec_to_simple_sds};
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

/// The second input of a set operation: the SBWT, and for all operations other than the
/// difference, the marks of its DBG and its coloring too. No k-mer of a difference is in the second
/// input, and its colors are those of the first input, so the difference needs only the SBWT of the
/// second input.
#[allow(clippy::large_enum_variant)] // Constructed once per set operation and moved right away
pub(crate) enum SecondInput<'c, CSS: ColorSetStorage> {
    Index(SbwtIndex<SubsetMatrix>, DbgMarks, &'c Coloring<CSS>),
    SbwtOnly(SbwtIndex<SubsetMatrix>),
}

/// Where the key k-mers that the colored subunitigs of the first input of a set operation require
/// come from (see mark_key_kmers).
#[allow(clippy::large_enum_variant)] // Constructed once per set operation and moved right away
pub(crate) enum FirstInputSubunitigs<'a> {
    /// A walk over the unitigs of its DBG.
    Walk(Dbg<'a, SubsetMatrix>),
    /// Its colored subunitigs, already known; with its SBWT.
    Known(&'a SbwtIndex<SubsetMatrix>, SubunitigBoundaries),
}

impl<'a> FirstInputSubunitigs<'a> {
    fn sbwt(&self) -> &'a SbwtIndex<SubsetMatrix> {
        match self {
            FirstInputSubunitigs::Walk(dbg) => dbg.sbwt(),
            FirstInputSubunitigs::Known(sbwt, _) => sbwt,
        }
    }
}

/// Computes the result of `operation` on the k-mers of the index with the SBWT `sbwt1`, the DBG
/// marks `marks1` and the coloring `coloring1`, and the index `input2`. The color set of each k-mer
/// of the result combines its color sets in the inputs as given by `colors`, and the result has the
/// colors `result_color_names`. Returns the SBWT, the DBG marks and the coloring of the result, and
/// the LCS array of its SBWT. The SBWTs are taken by value because the set operations of sbwt take
/// them in Arcs, and the marks because the DBGs of the inputs take them over.
///
/// For the intersection and the difference, sbwt1 needs select support. input2 is
/// [SecondInput::SbwtOnly] exactly for the difference.
#[allow(clippy::too_many_arguments)]
pub(crate) fn set_operation<CSS: ColorSetStorage + Send + Sync>(sbwt1: SbwtIndex<SubsetMatrix>, marks1: DbgMarks, coloring1: &Coloring<CSS>, input2: SecondInput<'_, CSS>, operation: SbwtOperation, colors: ColorCombination, result_color_names: &[String], optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> (SbwtIndex<SubsetMatrix>, DbgMarks, Coloring<CSS>, LcsArray) {
    let is_difference = matches!(operation, SbwtOperation::Difference);
    // The SBWTs go into Arcs for the set operations of sbwt, which take Arcs by design (see the
    // comment at sbwt::merge)
    let (sbwt2, rest2) = match input2 {
        SecondInput::Index(sbwt2, marks2, coloring2) => {
            assert!(!is_difference, "The difference needs only the SBWT of the second input");
            (sbwt2, Some((marks2, coloring2)))
        },
        SecondInput::SbwtOnly(sbwt2) => {
            assert!(is_difference, "Only the difference can do without the colors of the second input");
            (sbwt2, None)
        },
    };

    log::info!("Computing the sbwt merge plan");
    let interleaving = sbwt::MergeInterleaving::new(&sbwt1, &sbwt2, optimize_peak_ram, n_threads);
    assert_eq!(interleaving.s1.len(), interleaving.s2.len());

    let precalc_len = max(sbwt1.get_lookup_table().prefix_length, sbwt2.get_lookup_table().prefix_length);
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
            // The clone passed to sbwt is gone by now: as for the SBWTs below.
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

    // The unwraps are ok because sbwt's set operations currently drop all their clones of the Arcs before returning
    let sbwt1 = Arc::try_unwrap(sbwt1).unwrap();
    // Without the rest of the second input (the difference), its SBWT is not needed any more
    let sbwt2 = rest2.is_some().then_some(Arc::try_unwrap(sbwt2).unwrap());
    // The DBGs of the inputs take over their marks, so the marks are freed with the DBGs once the
    // key k-mers are marked
    let input1 = FirstInputSubunitigs::Walk(Dbg::from_marks(&sbwt1, marks1.0, marks1.1));
    let input2 = sbwt2.as_ref().zip(rest2).map(|(sbwt2, ((sufgroups2, dummy_marks2), coloring2))| (Dbg::from_marks(sbwt2, sufgroups2, dummy_marks2), coloring2));
    // Every k-mer of a merge is in one of the inputs, so it is visited when marking for the inputs
    // (and the check searches the result, which needs its LCS array)
    let check_visited = matches!(operation, SbwtOperation::Merge { .. });
    colored_result(input1, coloring1, input2, result_sbwt, &s1, &s2, is_leader, removed.as_ref(), operation.result_name(), colors, result_color_names, sample_distance, check_visited, n_threads)
}

/// Builds the colored index of the result of a set operation from its SBWT, which consists of the
/// positions of a merge plan of the inputs that were not removed, in order: `s1` and `s2` are the
/// positions of the first and the second input, and `removed` the removed positions (None if
/// nothing was removed). The first input has the coloring `coloring1`, and `input1` gives its SBWT
/// and its colored subunitigs: its DBG, walked for them, or the subunitigs if they are already
/// known. `input2` is the DBG and the coloring of the second input. The DBGs and the subunitigs are
/// needed only for the key k-mers, so they are freed once these are marked. The color sets of the
/// k-mers of the result combine their color sets in the inputs as given by `colors`, and the result
/// has the colors `result_color_names`. Without input2, its k-mers are not used for the key k-mers,
/// and `colors` must be ColorCombination::FromFirst. `is_leader` are the starts of the (k-1)-suffix
/// groups of the merge plan, from which the DBG of the result takes its groups (see
/// sbwt::result_k_minus_1_marks). `check_visited` (see mark_key_kmers) needs the walk over the
/// unitigs of the first input, so it can't be combined with known subunitigs. `result_name` is how
/// the result is called in the log. Returns the SBWT, the DBG marks and the coloring of the result,
/// and the LCS array of its SBWT.
#[allow(clippy::too_many_arguments)]
pub(crate) fn colored_result<'a, CSS: ColorSetStorage + Send + Sync>(input1: FirstInputSubunitigs<'a>, coloring1: &Coloring<CSS>, input2: Option<(Dbg<'a, SubsetMatrix>, &Coloring<CSS>)>, mut result_sbwt: SbwtIndex<SubsetMatrix>, s1: &Bits, s2: &Bits, is_leader: Bits, removed: Option<&Bits>, result_name: &str, colors: ColorCombination, result_color_names: &[String], sample_distance: usize, check_visited: bool, n_threads: usize) -> (SbwtIndex<SubsetMatrix>, DbgMarks, Coloring<CSS>, LcsArray) {
    // The color sets of the inputs need their SBWTs after their DBGs are gone
    let sbwt1 = input1.sbwt();
    let input2_colors = input2.as_ref().map(|(dbg2, coloring2)| (dbg2.sbwt(), *coloring2));

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
    let key_kmer_marks = mark_key_kmers(&input1, coloring1, input2.as_ref().map(|(dbg2, coloring2)| (dbg2, *coloring2)), s1, s2, removed, &result_sbwt, &result_lcs, &result_dbg, sample_distance, check_visited, n_threads);
    // Only the key k-mers need the DBGs of the inputs (their marks) and their known subunitigs.
    // The marks of the DBG of the result are kept for the result.
    drop((input1, input2));
    let result_marks = result_dbg.into_marks();
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / result_sbwt.n_kmers() as f64 * 100.0);

    let n_colors = u32::try_from(result_color_names.len()).unwrap_or_else( |_| {
        log::error!("Maximum number of colors 2^32 exceeded");
        panic!();
    });
    let gen = || ElementGeneratorFromInterleaving {
        s1,
        s2,
        removed,
        sbwt1,
        coloring1,
        input2: input2_colors,
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

    let coloring = Coloring::new(colex_map, css, Some(result_color_names));
    (result_sbwt, result_marks, coloring, result_lcs)
}

/// Marks the key k-mers of the result: the structural key k-mers of the result, and the k-mers of
/// the result where a colored subunitig of an input ends (see mark_key_kmers_for). In an
/// intersection, a color set changes between consecutive k-mers x and y of a unitig of the result
/// only if the color set of x or y changes in one of the inputs, which has the edge from x to y
/// too. Then x is the last k-mer of a colored subunitig of that input, so it is marked. In a
/// difference, the color sets are those of the first input, which has the edge from x to y too, so
/// marking for the first input suffices, and input2 is None.
///
/// With [FirstInputSubunitigs::Known], the key k-mers for the first input are marked from its known
/// colored subunitigs instead of a walk over its unitigs.
///
/// With check_visited, checks in debug builds that every k-mer of the result is visited while
/// marking for the inputs (see mark_key_kmers_for). Skipped in release builds: it costs a full
/// extra streaming-index pass over the result and is not needed for marking key k-mers, only for
/// catching bugs during development.
#[allow(clippy::too_many_arguments)]
fn mark_key_kmers<CSS: ColorSetStorage + Send + Sync>(input1: &FirstInputSubunitigs<'_>, coloring1: &Coloring<CSS>, input2: Option<(&Dbg<'_, SubsetMatrix>, &Coloring<CSS>)>, s1: &Bits, s2: &Bits, removed: Option<&Bits>, result_sbwt: &SbwtIndex<SubsetMatrix>, result_lcs: &LcsArray, result_dbg: &Dbg<'_, SubsetMatrix>, sample_distance: usize, check_visited: bool, n_threads: usize) -> bitvec::vec::BitVec {
    let k = result_sbwt.k();
    assert_eq!(k, input1.sbwt().k());
    assert!(input2.is_none_or(|(dbg2, _)| dbg2.sbwt().k() == k));

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
    let walked_input1 = match input1 {
        FirstInputSubunitigs::Walk(dbg1) => Some((dbg1, coloring1, s1)),
        FirstInputSubunitigs::Known(_, boundaries1) => {
            let to_result = InputToResult::new(s1, removed.as_ref());
            mark_key_kmers_from_boundaries(boundaries1, |colex| to_result.result_colex(colex), result_dbg, &key_kmer_marks, n_threads);
            None
        },
    };
    for (dbg, coloring, in_input) in walked_input1.into_iter().chain(input2.map(|(dbg2, coloring2)| (dbg2, coloring2, s2))) {
        let to_result = InputToResult::new(in_input, removed.as_ref());
        mark_key_kmers_for(dbg, coloring, |colex| to_result.result_colex(colex), result_sbwt, result_lcs, result_dbg, &key_kmer_marks, visited_marks.as_ref(), n_threads);
    }
    drop(removed);

    if let Some(visited_marks) = visited_marks {
        assert_eq!(visited_marks.into_bitvec().count_ones(), result_sbwt.n_kmers());
    }

    key_kmer_marks.into_bitvec()
}
