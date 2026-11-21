//! Marking the key k-mers of the result of a set operation: the k-mers whose color set is stored.

use sbwt::{dbg::{Dbg, Node}, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::{ops::{BitVec, Select}, sparse_vector::{SparseBuilder, SparseVector}};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage, unitig_export::break_to_colored_subunitigs};

pub(crate) fn mark_kmer(colex: usize, marks: &AtomicBitmap) {
    marks.set(colex, true);
}

pub(crate) fn mark_in_neighbors<'a>(colex: usize, dbg: &Dbg<'a, SubsetMatrix>, marks: &AtomicBitmap) {
    let mut in_neighbor_buf = Vec::<(Node, u8)>::new(); // TODO: avoid this allocation
    dbg.push_in_neighbors(Node{id: colex}, &mut in_neighbor_buf);
    for (in_node, _) in in_neighbor_buf.iter() {
        marks.set(in_node.id, true);
    }
}

type Bits = bitvec::vec::BitVec<u64, bitvec::order::Lsb0>;

// A colored subunitig of an input requires two kinds of key k-mers in the result: its last k-mer,
// and the in-neighbors in the DBG of the result of its first k-mer. `to_result` maps a colex
// position of the input to its colex position in the result, or None if the k-mer is not there.

/// Marks the key k-mer of the result that the last k-mer of a colored subunitig of an input, at
/// `last` in the input, requires.
fn mark_subunitig_last(last: usize, to_result: &impl Fn(usize) -> Option<usize>, key_kmer_marks: &AtomicBitmap) {
    if let Some(last) = to_result(last) {
        mark_kmer(last, key_kmer_marks);
    }
}

/// Marks the key k-mers of the result that the first k-mer of a colored subunitig of an input, at
/// `first` in the input, requires.
fn mark_subunitig_first(first: usize, to_result: &impl Fn(usize) -> Option<usize>, result_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap) {
    if let Some(first) = to_result(first) {
        mark_in_neighbors(first, result_dbg, key_kmer_marks);
    }
}

/// The first and the last k-mer of colored subunitigs of an input, as colex positions of the input,
/// as found in a walk over its unitigs (see break_to_colored_subunitigs). With these, the key
/// k-mers that the input requires in a result can be marked without walking its unitigs again (see
/// mark_key_kmers_from_boundaries). Subunitigs whose k-mers are not in the result can be left out,
/// and neighboring subunitigs that have the same color set in the result can be given as one (as
/// kmers_with_colors does for a subset), so that the result gets no key k-mers between them.
///
/// The positions are Elias-Fano coded (2 + log(n / m) bits for each of m positions in an input of
/// n nodes): a subunitig has 20 k-mers or so, so this is several times smaller than one bit per
/// node, and much smaller still when few subunitigs are kept.
pub(crate) struct SubunitigBoundaries {
    pub(crate) firsts: SparseVector,
    pub(crate) lasts: SparseVector,
}

impl SubunitigBoundaries {
    /// The positions of the ones of `bits`, Elias-Fano coded.
    pub(crate) fn sparse(bits: &Bits) -> SparseVector {
        let mut builder = SparseBuilder::new(bits.len(), bits.count_ones()).unwrap();
        builder.extend(bits.iter_ones());
        SparseVector::try_from(builder).unwrap()
    }
}

/// Marks the key k-mers that the colored subunitigs given by `boundaries` require, as
/// mark_key_kmers_for does for the colored subunitigs of an input that it walks.
pub(crate) fn mark_key_kmers_from_boundaries(boundaries: &SubunitigBoundaries, to_result: impl Fn(usize) -> Option<usize> + Sync, result_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, n_threads: usize) {
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    const CHUNK: usize = 1 << 12; // Boundaries per parallel task
    assert_eq!(boundaries.firsts.len(), boundaries.lasts.len());
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| {
        (0..boundaries.lasts.count_ones().div_ceil(CHUNK)).into_par_iter().for_each(|chunk| {
            for (_, colex) in boundaries.lasts.select_iter(chunk * CHUNK).take(CHUNK) {
                mark_subunitig_last(colex, &to_result, key_kmer_marks);
            }
        });
        (0..boundaries.firsts.count_ones().div_ceil(CHUNK)).into_par_iter().for_each(|chunk| {
            for (_, colex) in boundaries.firsts.select_iter(chunk * CHUNK).take(CHUNK) {
                mark_subunitig_first(colex, &to_result, result_dbg, key_kmer_marks);
            }
        });
    });
}

/// Marks the key k-mers of the result of a set operation that the colored subunitigs of `coloring`
/// (one of the inputs) require: the last k-mer of each colored subunitig, and the in-neighbors in
/// the DBG of the result of the first k-mer of each colored subunitig. `to_result` maps a colex
/// position of `coloring` to its colex position in the result, or None if the k-mer is not there.
///
/// With visited_marks, also marks there every k-mer of the result that is a k-mer of `coloring`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn mark_key_kmers_for<CSS: ColorSetStorage + Send + Sync>(coloring: &CompactColexKmers<CSS>, to_result: impl Fn(usize) -> Option<usize> + Sync, result_sbwt: &SbwtIndex<SubsetMatrix>, result_lcs: &LcsArray, result_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, visited_marks: Option<&AtomicBitmap>, n_threads: usize) {

    // Only the coverage check searches the result
    let result_si = visited_marks.map(|_| StreamingIndex::new(result_sbwt, result_lcs));

    let k = result_sbwt.k();
    assert_eq!(k, coloring.get_k());

    log::info!("Initializing DBG");
    let dbg = crate::set_operations::input_dbg(coloring, n_threads);

    log::info!("Iterating unitigs");
    let bar = indicatif::ProgressBar::new(coloring.sbwt().n_kmers() as u64);
    dbg.iter_unitigs_with_callback(|nodes|{
        let unitig_colex_ranks = nodes.iter().map(|v| v.id).collect::<Vec<usize>>(); // TODO: avoid this allocation
        let (_, subunitig_ranges) = break_to_colored_subunitigs(&unitig_colex_ranks, &[], coloring.get_map(), coloring.sbwt());

        for subunitig_range in subunitig_ranges {
            // (s,e) = (start of first k-mer, start of the k-mer after the last k-mer)
            let (s,e) = (subunitig_range.start, subunitig_range.end); 
            assert!(s < e);

            // The k-mer at position j of the unitig is node j. Its position in the result comes
            // from the interleaving of the inputs, so neither the unitig string nor a search in
            // the result is needed.
            mark_subunitig_last(nodes[e-1].id, &to_result, key_kmer_marks);
            mark_subunitig_first(nodes[s].id, &to_result, result_dbg, key_kmer_marks);
        }

        // Coverage check: record every k-mer of this unitig as visited, so that the caller can
        // assert afterwards that every k-mer of the result was visited while processing both
        // inputs. pipeline::mark_key_kmers does this for merges in debug builds only: it costs a
        // full streaming-index pass over the unitig. In a merge, it is otherwise unneeded for
        // marking: every k-mer that would be marked here is already marked by the loop above or
        // by mark_structural_key_kmers on the merged DBG (proof: only the second coloring can
        // have visited k-mers, and visited = contained in the first coloring. Let x be the first
        // k-mer of a visited run.
        // - If x is the first k-mer of this unitig, its in-neighbors were marked by the loop above.
        // - Otherwise its predecessor w in this unitig is not in the first coloring. If x has no
        //   in-neighbors in the first coloring, x starts a unitig of the first coloring, so the
        //   loop above marked its in-neighbors when processing the first coloring. Otherwise x
        //   has an in-neighbor v != w there, so x has in-degree >= 2 in the merged DBG and
        //   mark_structural_key_kmers marked its in-neighbors.
        // Let x be the last k-mer of a visited run, followed by y in this unitig. y is not in the
        // first coloring, so either x has out-degree 0 in the first coloring and ends a unitig
        // there (marked by the loop above), or x has out-degree >= 2 in the merged DBG and ends
        // a merged unitig (marked by mark_structural_key_kmers). Checked with an assertion on ~3700
        // random adversarial merges.)
        if let Some(visited_marks) = visited_marks {
            let mut unitig = Vec::<u8>::with_capacity(nodes.len() + k - 1);
            dbg.push_unitig_string(nodes, &mut unitig);
            for (match_len, colex_range) in result_si.as_ref().unwrap().matching_statistics_iter(&unitig).skip(k-1) {
                assert!(match_len == k);
                assert!(colex_range.len() == 1);

                // Modifying the same bitmap we are accessing is alright because we are iterating
                // unitigs and hence this k-mer will not be encountered a second time while
                // processing this coloring. So in the future if we find a visited-bit that is set
                // to 1, then it must have been set during the processing of some previous coloring.
                visited_marks.set(colex_range.start, true);
            }
        }
        bar.inc(nodes.len() as u64);
    }, n_threads);
    bar.finish();
}

/// Marks the key k-mers around branches of the DBG (the last k-mer of each unitig and the
/// in-neighbors of its first k-mer), and samples every sample_distance-th k-mer of each unitig.
/// This part is independent of coloring.
pub(crate) fn mark_structural_key_kmers(sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, sample_distance: usize, n_threads: usize) {
    let bar = indicatif::ProgressBar::new(sbwt.n_kmers() as u64);
    dbg.iter_unitigs_with_callback(|nodes|{
        mark_in_neighbors(nodes.first().unwrap().id, dbg, key_kmer_marks);
        mark_kmer(nodes.last().unwrap().id, key_kmer_marks);

        for v in nodes.iter().rev().step_by(sample_distance) {
            key_kmer_marks.set(v.id,  true);
        }
        bar.inc(nodes.len() as u64);
    }, n_threads);
    bar.finish();
}
