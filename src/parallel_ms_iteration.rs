use std::{collections::{hash_set::IntoIter, HashSet}, ops::Sub, path::PathBuf};

use rayon::iter::{IntoParallelIterator, ParallelBridge as _, ParallelIterator};
use sbwt::{dbg::Dbg, reverse_complement_in_place, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, io::ChainedInputStream, set_of_sets_construction::{ParallelElementGenerator, SetElement}};

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

pub struct DistinctColexComputation<'a> {
    streaming_index: sbwt::StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    input: jseqio::reader::DynamicFastXReader,
    color: usize,
    set_ids: DeduplicatingBuffer,
}

impl<'a> DistinctColexComputation<'a> {
    pub fn new(sbwt: &'a sbwt::SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, input: jseqio::reader::DynamicFastXReader, color: usize) -> Self {
        let streaming_index = StreamingIndex::new(sbwt, lcs); 
        Self {
            streaming_index,
            input,
            color,
            set_ids: DeduplicatingBuffer::new(sbwt.n_sets()),
        }
    }

    fn process_seq(&mut self, cur_seq: &[u8]) {
        let k = self.streaming_index.k();

        let ms_iter = self.streaming_index.matching_statistics_iter(cur_seq);
        for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
            assert!(colex.len() == 1);
            self.set_ids.insert(colex.start);
        }
    }

    fn run(mut self) -> DeduplicatingBuffer {
        let mut buf = Vec::<u8>::new();
        while let Some(rec) = self.input.read_next_mut().unwrap() {
            buf.clear();
            buf.extend_from_slice(rec.seq);
            self.process_seq(&buf);
            reverse_complement_in_place(&mut buf);
            self.process_seq(&buf);
        }
        self.set_ids
    }
}

pub struct DeduplicatingColorElementGenerator<'a> {
    input_files: Vec<PathBuf>,
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    lcs: &'a LcsArray,
    filter: Option<simple_sds_sbwt::bit_vector::BitVector>,
}

impl<'a> DeduplicatingColorElementGenerator<'a> {
    pub fn new( sbwt: &'a SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, input_files: Vec<PathBuf>) -> Self {
        Self { input_files, sbwt, lcs, filter: None }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for DeduplicatingColorElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(|| {
            self.input_files.iter().enumerate().par_bridge().for_each(|(color, file_path)| {
                log::info!("Processing color {}", color);
                let reader = jseqio::reader::DynamicFastXReader::from_file(&file_path).unwrap();
                let gen = DistinctColexComputation::new(&self.sbwt, &self.lcs, reader, color);
                let distinct_colex_positions = gen.run();
                for colex in distinct_colex_positions.into_iter() {
                    let set_id = if let Some(filter) = &self.filter {
                        if !filter.get(colex) {
                            continue; // Do not report this
                        } else {
                            // Assign new id
                            filter.rank(colex)
                        }
                    } else {
                        colex // No filter
                    };

                    callback(SetElement{
                        set_id,
                        color,
                    });
                }
            })
        });
    }

    fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector) {
        self.filter = Some(filter)
    }
}

pub struct MergedElementGenerator<'a, CSS: ColorSetStorage + Sync + Send> {
    pub merged_sbwt: &'a SbwtIndex<SubsetMatrix>,
    pub merged_lcs: &'a LcsArray,
    pub coloring1: &'a CompactColexKmers<CSS>,
    pub coloring2: &'a CompactColexKmers<CSS>,
    pub dbg1: &'a Dbg<'a, SubsetMatrix>,
    pub dbg2: &'a Dbg<'a, SubsetMatrix>,
    pub filter: Option<simple_sds_sbwt::bit_vector::BitVector>, // With rank support
}

impl<'a, CSS: ColorSetStorage + Sync + Send> MergedElementGenerator<'a, CSS> {

    fn process_dbg<'b>(&'b self, dbg: &Dbg<'a, SubsetMatrix>, coloring: &CompactColexKmers<CSS>, color_offset: usize, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        let si = StreamingIndex::new(&self.merged_sbwt, &self.merged_lcs);
        let k = self.merged_sbwt.k();
        assert!(self.merged_sbwt.k() == k);
        let bar = indicatif::ProgressBar::new(coloring.sbwt().n_kmers() as u64);
        dbg.iter_unitigs_with_callback(|nodes, unitig| {
            let mut merged_colexes_iter = si.matching_statistics_iter(&unitig)
            .skip(k-1)
            .map(|(match_len, range)| {
                assert!(match_len == k, "k-mer not found in merged sbwt"); // TODO: print which k-mer it is
                assert!(range.len() == 1);
                range.start
            });

            for node in nodes {
                let colex = node.id;
                let merged_colex = merged_colexes_iter.next().expect("Programming mistake: merged colex iter has fewer elements than DBG nodes");

                let set_id = if let Some(filter) = &self.filter {
                    if !filter.get(merged_colex) {
                        continue // Do not report this
                    } else {
                        // Assign new id
                        let new_id = filter.rank(merged_colex);
                        new_id
                    }
                } else {
                    merged_colex // No filter
                };

                let cs = coloring.colex_to_set(colex);
                for color in cs.iter() {
                    callback(SetElement{set_id, color: color + color_offset}) 
                }
            }
            bar.inc(nodes.len() as u64);
        }, n_threads);
        bar.finish();

    }

}

impl<'a, CSS: ColorSetStorage + Sync + Send> ParallelElementGenerator for MergedElementGenerator<'a, CSS> {

    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        self.process_dbg(self.dbg1, self.coloring1, 0, &callback, n_threads);
        self.process_dbg(self.dbg2, self.coloring2, self.coloring1.get_set_storage().n_colors(), &callback, n_threads);
    }

    fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector) {
        self.filter = Some(filter);
    }
}