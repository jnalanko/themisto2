use std::{collections::HashSet, path::PathBuf};

use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix, reverse_complement_in_place};

use crate::set_of_sets_construction::SetElement;

pub struct MsElementGenerator<'a> {
    input_files: Vec<PathBuf>,
    streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
}

impl<'a> MsElementGenerator<'a> {
    pub fn new(
        input_files: Vec<PathBuf>,
        streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    ) -> Self {
        Self {
            input_files,
            streaming_index,
        }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for MsElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        let k = self.streaming_index.k();
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(|| {
            self.input_files.par_iter().enumerate().for_each(|(color, file_path)| {
                log::info!("Processing color {}", color);
                let mut reader = jseqio::reader::DynamicFastXReader::from_file(&file_path).unwrap();
                while let Some(rec) = reader.read_next_mut().unwrap() {
                    let ms_iter = self.streaming_index.matching_statistics_iter(rec.seq);
                    for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
                        assert!(colex.len() == 1);
                        callback(SetElement{
                            set_id: colex.start,
                            color,
                        });
                    }

                    reverse_complement_in_place(rec.seq);

                    let ms_iter = self.streaming_index.matching_statistics_iter(rec.seq);
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
