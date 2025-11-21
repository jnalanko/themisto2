use std::{collections::HashSet, path::PathBuf};

use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelBridge as _, ParallelIterator};
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix, reverse_complement_in_place};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::set_of_sets_construction::SetElement;

pub struct MsElementGenerator<'a> {
    input_files: Vec<PathBuf>,
    streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    filter: Option<simple_sds_sbwt::bit_vector::BitVector>,
}

impl<'a> MsElementGenerator<'a> {
    pub fn new(
        input_files: Vec<PathBuf>,
        streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    ) -> Self {
        Self {
            input_files,
            streaming_index,
            filter: None,
        }
    }
}

impl<'a> MsElementGenerator<'a> {
    fn run_seq(&self, seq: &[u8], color: usize, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync) {
        let k = self.streaming_index.k();
        let ms_iter = self.streaming_index.matching_statistics_iter(seq);
        let kmer_iter = ms_iter.skip(k-1).filter(|(len, _colex)| *len == k);
        let filtered_iter = kmer_iter.filter_map(|(_, colex)| {
            assert!(colex.len() == 1);
            let set_id = colex.start;
            if let Some(filter) = &self.filter {
                if !filter.get(set_id) {
                    None // Do not report this
                } else {
                    // Assign new id
                    let new_id = filter.rank(set_id);
                    Some(new_id)
                }
            } else {
                Some(set_id) // No filter
            }
        });

        for id in filtered_iter {
            callback(SetElement{
                set_id: id,
                color,
            });
        }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for MsElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(|| {
            self.input_files.iter().enumerate().par_bridge().for_each(|(color, file_path)| {
                log::info!("Processing color {}", color);
                let mut reader = jseqio::reader::DynamicFastXReader::from_file(&file_path).unwrap();
                while let Some(rec) = reader.read_next_mut().unwrap() {
                    self.run_seq(rec.seq, color, &callback);
                    reverse_complement_in_place(rec.seq);
                    self.run_seq(rec.seq, color, &callback);
                }
            })
        });
    }
    
    fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector) {
        self.filter = Some(filter);
    }

    
}
