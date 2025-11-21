use std::collections::HashSet;

use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::io::ChainedInputStream;

trait ColoredSeqStreams {

    /// Returns pairs (color, sequence)
    fn next(&mut self) -> Option<(usize, &[u8])>;
}

// The callback is called on pairs (color, colex)
// The input streams must first iterate all sets of color 0, then all sets of color 1, and so on.
pub fn generate_all_color_set_elements(streaming_index: &StreamingIndex<SbwtIndex<SubsetMatrix>, LcsArray>, mut input: impl ColoredSeqStreams, n_threads: usize, mut callback: impl FnMut(usize, usize)) {

    let k = streaming_index.k();
    let mut colex_ranks = HashSet::<usize>::new(); // For the current color
    let mut cur_color = 0_usize;
    while let Some((color, seq)) = input.next() {
        if color != cur_color {
            for colex in colex_ranks.iter() {
                callback(cur_color, *colex);
            }
            colex_ranks = HashSet::<usize>::new();
            cur_color = color;
        }
        let ms_iter = streaming_index.matching_statistics_iter(seq);
        for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
            assert!(colex.len() == 1);
            colex_ranks.insert(colex.start);
        }
    }

    // Final color
    for colex in colex_ranks.iter() {
        callback(cur_color, *colex);
    }

    /*
    elements.sort_by(|a, b| (a.color, a.set_id).cmp(&(b.color, b.set_id)));
    elements.dedup();
    Self { elements, pos: 0}
    */
}