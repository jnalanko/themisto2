use std::{collections::HashSet, path::PathBuf};

use rayon::iter::IntoParallelIterator;
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::set_of_sets_construction::SetElement;

struct MsElementGenerator<'a> {
    input_files: Vec<PathBuf>,
    streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for MsElementGenerator<'a> {
    fn run(&mut self, callback: impl FnMut(crate::set_of_sets_construction::SetElement), n_threads: usize) {
        let k = self.streaming_index.k();
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        thread_pool.install(|| {
            self.input_files.into_iter().enumerate().into_par_iter().for_each(|(color, file_path)| {
                log::info!("Processing color {}", color);
                let reader = jseqio::reader::DynamicFastXReader::from_file(file_path).unwrap();
                while let Some(seq) = reader.read_next_mut() {
                    let ms_iter = self.streaming_index.matching_statistics_iter(seq);
                    for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
                        assert!(colex.len() == 1);
                        callback(SetElement{
                            set_id: colex.start,
                            color,
                        });
                    }
                }
            })
        });
    }
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