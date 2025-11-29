use std::{collections::{hash_set::IntoIter, HashSet}, path::PathBuf};

use bitvec::{order::Lsb0, slice::IterOnes};
use rayon::iter::{ParallelBridge as _, ParallelIterator};
use sbwt::{reverse_complement_in_place, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{io::ChainedInputStream, set_of_sets_construction::SetElement};

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

struct DeduplicatingBuffer {
    universe_size: usize,
    hashset: Option<HashSet<usize>>,
    bitmap: Option<bitvec::vec::BitVec>, // Switch to this when the set is large enough
    empty: bool,
} 

impl DeduplicatingBuffer {
    fn new(universe_size: usize) -> Self {
        Self {
            universe_size,
            hashset: Some(HashSet::new()),
            bitmap: None,
            empty: true,
        }
    }

    fn insert(&mut self, id: usize) {
        self.empty = false;
        if let Some(hs) = &mut self.hashset {
            hs.insert(id);
            if hs.capacity() > self.universe_size / 64 {
                // Switch to bitmap
                let mut bv = bitvec::bitvec![0; self.universe_size];
                for &elem in hs.iter() {
                    bv.set(elem, true);
                }
                self.bitmap = Some(bv);
                self.hashset = None;
            }
        } else if let Some(bv) = &mut self.bitmap {
            bv.set(id, true);
        } else {
            panic!("Both hashset and bitmap are None");
        }
    }

    fn into_iter(self) -> DeduplicatingBufferIter {
        if let Some(hs) = self.hashset {
            DeduplicatingBufferIter {
                hashset_iter: Some(hs.into_iter()),
                bitmap_iter: None,
            }
        } else if let Some(bv) = self.bitmap {
            DeduplicatingBufferIter {
                hashset_iter: None,
                bitmap_iter: Some(OwningBitVecOnesIterator::new(bv)),
            }
        } else {
            panic!("Both hashset and bitmap are None");
        }
    }
}

struct OwningBitVecOnesIterator {
    bv: bitvec::vec::BitVec,
    pos: usize,
}

impl OwningBitVecOnesIterator {
    fn new(bv: bitvec::vec::BitVec) -> Self {
        Self { bv, pos: 0 }
    }
}

impl Iterator for OwningBitVecOnesIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(offset) = self.bv[self.pos..].first_one() {
            let ret = self.pos + offset;
            self.pos += offset + 1; // Starting point of the next iteration
            Some(ret)
        } else {
            None
        }
    }
}   

struct DeduplicatingBufferIter {
    hashset_iter: Option<std::collections::hash_set::IntoIter<usize>>,
    bitmap_iter: Option<OwningBitVecOnesIterator>,
}

impl<'a> Iterator for DeduplicatingBufferIter {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(hs_iter) = &mut self.hashset_iter {
            hs_iter.next()
        } else if let Some(bv_iter) = &mut self.bitmap_iter {
            bv_iter.next()
        } else {
            panic!("Both hashset_iter and bitmap_iter are None");
        }
    }
}

pub struct DeduplicatingColorElementGenerator<'a> {
    streaming_index: sbwt::StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    input: ChainedInputStream,
    cur_color: usize,
    cur_set_ids: DeduplicatingBuffer,
    output_buf: (usize, DeduplicatingBuffer), // (Color, set ids)
    filter: Option<simple_sds_sbwt::bit_vector::BitVector> // Bit vector with rank support
}

impl<'a> DeduplicatingColorElementGenerator<'a> {
    pub fn new(sbwt: &'a sbwt::SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, input: ChainedInputStream) -> Self {
        let streaming_index = StreamingIndex::new(sbwt, lcs); 
        Self {
            streaming_index,
            input,
            cur_color: 0,
            cur_set_ids: DeduplicatingBuffer::new(sbwt.n_sets()),
            output_buf: (0, DeduplicatingBuffer::new(sbwt.n_sets())),
            filter: None,
        }
    }

    fn process_current_seq_in_input(&mut self) {
        let seq = self.input.get_seq_buf_mut();
        let k = self.streaming_index.k();

        let ms_iter = self.streaming_index.matching_statistics_iter(seq);
        for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
            assert!(colex.len() == 1);
            self.cur_set_ids.insert(colex.start);
        }

        reverse_complement_in_place(seq);

        let ms_iter = self.streaming_index.matching_statistics_iter(seq);
        for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
            assert!(colex.len() == 1);
            self.cur_set_ids.insert(colex.start);
        }
    }

    fn hash_set_to_output_buf(&mut self) {
        // Move the set ids to the out buffer
        std::mem::swap(&mut self.cur_set_ids, &mut self.output_buf.1);
        self.output_buf.0 = self.cur_color;

        self.cur_set_ids = DeduplicatingBuffer::new(self.streaming_index.sbwt_len()); // Clear the buffer
        log::info!("Searched color {}", self.cur_color);
        self.cur_color += 1;
    }
}

impl<'a> Iterator for DeduplicatingColorElementGenerator<'a> {
    type Item = SetElement;

    fn next(&mut self) -> Option<Self::Item> {

        if !self.output_buf.1.empty {
            aaa
        }
        if let Some(id) = self.output_buf.1.pop() {
            return Some(SetElement { set_id: id, color: self.output_buf.0 });
        }
        if self.input.done() { return None }

        // Read and process all sequences of the current color
        loop {
            if self.input.stream_next().is_some() {
                let color = self.input.cur_file_idx();
                if color == self.cur_color {
                    self.process_current_seq_in_input();
                } else {
                    // Push to output buffer and start returning
                    self.hash_set_to_output_buf();
                    self.process_current_seq_in_input();
                    return self.next();
                }
            } else {
                // End of input. Push the set ids of the last color to the output buffer
                self.hash_set_to_output_buf();
                return self.next();
            }
        }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for DeduplicatingColorElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        // TODO: make this multithreaded
        while let Some(mut elem) = self.next() {
            if let Some(filter) = &self.filter {
                if !filter.get(elem.set_id) {
                    continue; // This is filtered away
                } else {
                    // Keep and assign new id
                    let new_id = filter.rank(elem.set_id);
                    elem.set_id = new_id;
                }
            }
            callback(elem);
        }
    }

    fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector) {
        self.filter = Some(filter)
    }
}