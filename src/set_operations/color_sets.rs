//! Generating the color set elements of the key k-mers of the result of a set operation, from the
//! color sets of the two input indexes, for building the color set storage of the result.

use std::{ops::Range, sync::Arc};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, set_of_sets_construction::{ParallelElementGenerator, SetElement}};

/// How the color sets of a k-mer in the two input indexes are combined into its color set in the
/// result of a set operation.
#[derive(Clone, Copy)]
pub enum ColorCombination<'a> {
    /// The union of the two sets. Maps each color of coloring2 to its color id in the result: the
    /// colors of coloring1 keep their ids, and ids smaller than the number of colors of coloring1
    /// are colors shared with coloring1. The k-mer may be in only one of the inputs.
    Union { color2_to_result: &'a [usize] },
    /// The intersection of the two sets. Maps each color of coloring1 and coloring2 to its color
    /// id in the result, or to None if the color is not in the result. The k-mer must be in both
    /// inputs.
    Intersection { color1_to_result: &'a [Option<usize>], color2_to_result: &'a [Option<usize>] },
    /// The set of coloring1. Maps each color of coloring1 to its color id in the result, or to None
    /// if the color is not in the result; without a mapping, the colors keep their ids. The k-mer
    /// must be in coloring1.
    FromFirst { color1_to_result: Option<&'a [Option<usize>]> },
}

/// Generates the color set elements of the key k-mers of the result of a set operation on the
/// SBWTs of coloring1 and coloring2. The result consists of the positions of an interleaving of the
/// two SBWTs that were not removed, in order (see sbwt::merge_with_removed_positions and
/// sbwt::intersect_with_removed_positions): the j-th position that was not removed is result colex
/// position j.
pub struct ElementGeneratorFromInterleaving<'a, CSS: ColorSetStorage + Sync + Send> {
    pub s1: &'a bitvec::vec::BitVec<u64, bitvec::order::Lsb0>, // The positions of the interleaving in coloring1
    pub s2: &'a bitvec::vec::BitVec<u64, bitvec::order::Lsb0>, // The positions of the interleaving in coloring2
    // The positions of the interleaving that are not in the result. None if nothing was removed,
    // in which case interleaving positions are result colex positions.
    pub removed: Option<&'a bitvec::vec::BitVec<u64, bitvec::order::Lsb0>>,
    pub coloring1: &'a CompactColexKmers<CSS>,
    pub coloring2: Option<&'a CompactColexKmers<CSS>>, // Needed unless the colors are ColorCombination::FromFirst
    pub result_key_kmer_marks: &'a bitvec::vec::BitVec, // Only reporting set elements for these
    pub filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>, // With rank support
    pub colors: ColorCombination<'a>,
    pub n_result_colors: usize,
}

struct ThreadInput {
    interleaving_range: Range<usize>,
    s1_start_rank: usize,
    s2_start_rank: usize,
    result_start_colex: usize, // Result colex position of the first position of interleaving_range that was not removed
}

impl<'a, CSS: ColorSetStorage + Sync + Send> ParallelElementGenerator for ElementGeneratorFromInterleaving<'a, CSS> {

    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        let (s1, s2, removed) = (self.s1, self.s2, self.removed);
        let n = s1.len();
        assert!(s2.len() == n && removed.is_none_or(|r| r.len() == n));

        // The number of bits of s1 and s2 and of result nodes in each piece, counted in parallel
        let ranges = crate::util::segment_range(0..n, n_threads);
        let counts: Vec<(usize, usize, usize)> = ranges.clone().into_par_iter().map(|range| {
            let n_result_nodes = removed.map_or(range.len(), |r| r[range.clone()].count_zeros());
            (s1[range.clone()].count_ones(), s2[range.clone()].count_ones(), n_result_nodes)
        }).collect();
        let mut thread_inputs = Vec::<ThreadInput>::with_capacity(n_threads);
        let (mut n_bits_s1, mut n_bits_s2, mut n_result_nodes) = (0_usize, 0_usize, 0_usize);
        for (range, (c1, c2, c3)) in ranges.into_iter().zip(counts) {
            thread_inputs.push(ThreadInput {
                interleaving_range: range,
                s1_start_rank: n_bits_s1,
                s2_start_rank: n_bits_s2,
                result_start_colex: n_result_nodes,
            });
            n_bits_s1 += c1;
            n_bits_s2 += c2;
            n_result_nodes += c3;
        }
        assert_eq!(n_result_nodes, self.result_key_kmer_marks.len());

        let n_colors1 = self.coloring1.get_set_storage().n_colors();
        let n_colors2 = self.coloring2.map(|coloring2| coloring2.get_set_storage().n_colors());
        // The colors of the result that come from coloring1 and may also come from coloring2, so
        // that they must be reported only once (union) or only if they come from both (intersection).
        let mut is_shared = bitvec::bitvec![0; self.n_result_colors];
        match self.colors {
            ColorCombination::Union { color2_to_result } => {
                assert_eq!(Some(color2_to_result.len()), n_colors2);
                for &color in color2_to_result.iter().filter(|&&c| c < n_colors1) {
                    is_shared.set(color, true);
                }
            },
            ColorCombination::Intersection { color1_to_result, color2_to_result } => {
                assert_eq!(color1_to_result.len(), n_colors1);
                assert_eq!(Some(color2_to_result.len()), n_colors2);
                is_shared.fill(true);
            },
            ColorCombination::FromFirst { color1_to_result: Some(color1_to_result) } => assert_eq!(color1_to_result.len(), n_colors1),
            ColorCombination::FromFirst { color1_to_result: None } => assert_eq!(self.n_result_colors, n_colors1),
        }
        let has_shared_colors = is_shared.any();

        let bar = indicatif::ProgressBar::new(n as u64);
        thread_inputs.into_par_iter().for_each(|input| {
            let mut s1_colex = input.s1_start_rank;
            let mut s2_colex = input.s2_start_rank;
            let mut result_colex = input.result_start_colex;
            // The shared colors of coloring1 at the current k-mer, as a bitmap for O(1) lookups
            // and as a list for clearing the bitmap afterwards. Not allocated without shared colors.
            let mut in_set1 = bitvec::bitvec![0; if has_shared_colors { self.n_result_colors } else { 0 }];
            let mut in_set1_list = Vec::<usize>::new();
            for pos in input.interleaving_range {
                if pos > 0 && pos % 10000 == 0 {
                    bar.inc(10000);
                }
                if !removed.is_some_and(|r| r[pos]) {
                    if self.result_key_kmer_marks[result_colex] {
                        let (in_1, in_2) = (s1[pos], s2[pos]);
                        if let Some(new_set_id) = maybe_apply_filter(self.filter.as_deref(), result_colex) {
                            match self.colors {
                                ColorCombination::Union { color2_to_result } => {
                                    // A shared color may be in both sets. It must be reported only once.
                                    let dedup = has_shared_colors && in_1 && in_2;
                                    if in_1 {
                                        for color in self.coloring1.colex_to_set(s1_colex).iter() {
                                            callback(SetElement{set_id: new_set_id, color});
                                            if dedup && is_shared[color] {
                                                in_set1.set(color, true);
                                                in_set1_list.push(color);
                                            }
                                        }
                                    }
                                    if in_2 {
                                        for color in self.coloring2.unwrap().colex_to_set(s2_colex).iter() {
                                            let color = color2_to_result[color];
                                            if !(dedup && in_set1[color]) {
                                                callback(SetElement{set_id: new_set_id, color});
                                            } // Else already reported from coloring1
                                        }
                                    }
                                },
                                ColorCombination::Intersection { color1_to_result, color2_to_result } => {
                                    assert!(in_1 && in_2, "a key k-mer of an intersection is not in both inputs");
                                    for color in self.coloring1.colex_to_set(s1_colex).iter().filter_map(|c| color1_to_result[c]) {
                                        in_set1.set(color, true);
                                        in_set1_list.push(color);
                                    }
                                    for color in self.coloring2.unwrap().colex_to_set(s2_colex).iter().filter_map(|c| color2_to_result[c]) {
                                        if in_set1[color] {
                                            callback(SetElement{set_id: new_set_id, color});
                                        }
                                    }
                                },
                                ColorCombination::FromFirst { color1_to_result } => {
                                    assert!(in_1, "a key k-mer of the result is not in the first input");
                                    let set1 = self.coloring1.colex_to_set(s1_colex);
                                    match color1_to_result {
                                        None => for color in set1.iter() {
                                            callback(SetElement{set_id: new_set_id, color});
                                        },
                                        Some(color1_to_result) => for color in set1.iter().filter_map(|c| color1_to_result[c]) {
                                            callback(SetElement{set_id: new_set_id, color});
                                        },
                                    }
                                },
                            }
                            for color in in_set1_list.drain(..) {
                                in_set1.set(color, false);
                            }
                        }
                    }
                    result_colex += 1;
                }
                s1_colex += s1[pos] as usize;
                s2_colex += s2[pos] as usize;
            }
        });
        bar.finish();
    }

    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter.clone());
    }

    fn rewind(&mut self) {
        // Nothing needs to done, calling run() again already works
    }
}

/// Returns the id of the set of the given colex position after filtering: its rank among the
/// positions marked in the filter, or None if it is not marked. Without a filter, the colex
/// position itself.
fn maybe_apply_filter(filter: Option<&simple_sds_sbwt::bit_vector::BitVector>, colex: usize) -> Option<usize> {
    match filter {
        Some(filter) if !filter.get(colex) => None, // Do not report this
        Some(filter) => Some(filter.rank(colex)), // Assign new id
        None => Some(colex), // No filter
    }
}
