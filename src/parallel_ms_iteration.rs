use std::{collections::HashSet, ops::Range, sync::Arc};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sbwt::{LcsArray, MergeInterleaving, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix, reverse_complement_in_place};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, io::{self, RewindableSeqStreamGenerator}, set_of_sets_construction::{ParallelElementGenerator, SetElement}};

pub struct MsElementGenerator<'a> {
    color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>,
    streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>,
    include_rev_comp: bool,
}

impl<'a> MsElementGenerator<'a> {
    pub fn new(color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>, streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>, include_rev_comp: bool) -> Self {
        Self {
            color_stream_generator,
            streaming_index,
            filter: None,
            include_rev_comp,
        }
    }
}

pub struct WorkBatch<'a> {
    seq_concat: Vec<u8>,
    seq_ends: Vec<usize>,
    seq_colors: Vec<usize>,
    streaming_index: &'a StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    filter: &'a Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>,
    include_rev_comp: bool,
}

impl<'a> WorkBatch<'a> {
    fn run(mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync) {
        let n_seqs = self.seq_ends.len();
        assert_eq!(n_seqs, self.seq_colors.len());

        let mut seq_start = 0_usize;
        for seq_idx in 0..n_seqs {
            let seq_end = self.seq_ends[seq_idx];
            let seq = &self.seq_concat[seq_start..seq_end];
            let color = self.seq_colors[seq_idx];

            self.process_seq(seq, color, &callback);
            if self.include_rev_comp {
                // Reverse-complemet in-place. This is already since we now own `self`, so
                // nobody else can see the sequences anymore.
                let seq_mut = &mut self.seq_concat[seq_start..seq_end];
                reverse_complement_in_place(seq_mut);
                let seq = &self.seq_concat[seq_start..seq_end];
                self.process_seq(seq, color, &callback);
            }

            seq_start = seq_end;
        }
    }

    fn push_seq(&mut self, seq: &[u8], color: usize) {
        self.seq_concat.extend_from_slice(seq);
        self.seq_ends.push(self.seq_concat.len());
        self.seq_colors.push(color);
    }

    fn size_in_bytes(&self) -> usize {
        self.seq_concat.len() + self.seq_ends.len()*size_of::<usize>() + self.seq_colors.len()*size_of::<usize>()
    }

    // Don't call this from outside. Call self.run instead.
    fn process_seq(&self, seq: &[u8], color: usize, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync) {
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

    fn new(streaming_index: &'a StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>, filter: &'a Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>, include_rev_comp: bool) -> WorkBatch<'a> {
        WorkBatch { seq_concat: vec![], seq_ends: vec![], seq_colors: vec![], streaming_index, filter, include_rev_comp }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for MsElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        let (sender, receiver) = crossbeam::channel::bounded::<WorkBatch>(2*n_threads);
        let receiver_ref = &receiver; // To capture a reference

        // Here we need to get a bit tricky to avoid mutable aliasing of self. The issue is that
        // the producer thread needs a mutable reference to the ReWindableSeqStreamGenerator at self,
        // while the consumers need non-mutable access to the rest of self. This is not possible
        // at the same time because borrowing self borrows everything. The workaround is that we
        // swap in a dummy generator into self, so that we can get separate ownership of the generator
        // and pass it into the producer. In the end we swap it back in.
        let mut dummy_color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(io::EmptyRewindableSeqStreamGenerator{});
        std::mem::swap(&mut self.color_stream_generator, &mut dummy_color_stream_generator);
        let mut color_stream_generator = dummy_color_stream_generator; // Now this function owns this

        std::thread::scope(|scope| {
            // Channel of pairs (color id, seq stream)
            let producer_handle = scope.spawn(|| {
                let batch_size: usize = 1 << 23; // 8 MiB
                let mut color = 0_usize;
                let mut cur_batch = WorkBatch::new(&self.streaming_index, &self.filter, self.include_rev_comp);
                while let Some((mut color_stream, _stream_idx)) = color_stream_generator.next() {
                    log::info!("Processing color {}", color);
                    while let Some(seq) = color_stream.stream_next() {
                        cur_batch.push_seq(seq, color);
                        if cur_batch.size_in_bytes() >= batch_size {
                            sender.send(cur_batch).unwrap();
                            cur_batch = WorkBatch::new(&self.streaming_index, &self.filter, self.include_rev_comp);
                        }
                    }
                    color += 1;
                }
                sender.send(cur_batch).unwrap(); // Last batch. May be empty, but that is ok.
                drop(sender); // Finished
            });

            let consumer_handles: Vec<_> = (0..n_threads).map(|_| {
                scope.spawn(|| {
                    while let Ok(batch) = receiver_ref.recv() {
                        batch.run(&callback);
                    }
                })
            }).collect();

            // Wait for threads to finish
            producer_handle.join().unwrap();
            for h in consumer_handles { h.join().unwrap(); }
        });

        self.color_stream_generator = color_stream_generator; // Put this back in (see comment at the start of the function)
    }
    
    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter);
    }

    fn rewind(&mut self) {
        self.color_stream_generator.rewind();
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

impl Iterator for DeduplicatingBufferIter {
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
    input: Box<dyn SeqStream + Sync + Send>,
    set_ids: DeduplicatingBuffer,
}

impl<'a> DistinctColexComputation<'a> {
    pub fn new(sbwt: &'a sbwt::SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, input: Box<dyn SeqStream + Sync + Send>) -> Self {
        let streaming_index = StreamingIndex::new(sbwt, lcs); 
        Self {
            streaming_index,
            input,
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

    fn run(mut self, include_rev_comp: bool) -> DeduplicatingBuffer {
        let mut buf = Vec::<u8>::new();
        while let Some(seq) = self.input.stream_next() {
            buf.clear();
            buf.extend_from_slice(seq);
            self.process_seq(&buf);
            if include_rev_comp {
                reverse_complement_in_place(&mut buf);
                self.process_seq(&buf);
            }
        }
        self.set_ids
    }
}

pub struct DeduplicatingColorElementGenerator<'a> {
    color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>,
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    lcs: &'a LcsArray,
    filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>,
    include_rev_comp: bool,
}

impl<'a> DeduplicatingColorElementGenerator<'a> {
    pub fn new( sbwt: &'a SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>, include_rev_comp: bool) -> Self {
        Self { color_stream_generator, sbwt, lcs, filter: None, include_rev_comp }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for DeduplicatingColorElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        // TODO: a lot of this code is duplicated with  MsElementGenerator.

        let (sender, receiver) = crossbeam::channel::bounded::<(usize, Box<dyn SeqStream + Send + Sync>)>(2*n_threads);
        let receiver_ref = &receiver; // To capture a reference

        // Here we need to get a bit tricky to avoid mutable aliasing of self. The issue is that
        // the producer thread needs a mutable reference to the ReWindableSeqStreamGenerator at self,
        // while the consumers need non-mutable access to the rest of self. This is not possible
        // at the same time because borrowing self borrows everything. The workaround is that we
        // swap in a dummy generator into self, so that we can get separate ownership of the generator
        // and pass it into the producer. In the end we swap it back in.
        let mut dummy_color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(io::EmptyRewindableSeqStreamGenerator{});
        std::mem::swap(&mut self.color_stream_generator, &mut dummy_color_stream_generator);
        let mut color_stream_generator = dummy_color_stream_generator; // Now this function owns this

        std::thread::scope(|scope| {
            // Channel of pairs (color id, seq stream)
            let producer_handle = scope.spawn(|| {
                let mut color = 0_usize;
                while let Some((color_stream, _stream_idx)) = color_stream_generator.next() {
                    sender.send((color, color_stream)).unwrap();
                    color += 1;
                }
                drop(sender); // Finished
            });

            let consumer_handles: Vec<_> = (0..n_threads).map(|_| {
                scope.spawn(|| {
                    while let Ok((color, color_stream)) = receiver_ref.recv() {
                        log::info!("Processing color {}", color);
                        let gen = DistinctColexComputation::new(self.sbwt, self.lcs, color_stream);

                        let distinct_colex_positions = gen.run(self.include_rev_comp);
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
                    }
                })
            }).collect();

            // Wait for threads to finish
            producer_handle.join().unwrap();
            for h in consumer_handles { h.join().unwrap(); }
        });

        self.color_stream_generator = color_stream_generator; // Put this back in (see comment at the start of the function)
    }

    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter.clone())
    }

    fn rewind(&mut self) {
        self.color_stream_generator.rewind();
    }
}

pub struct ElementGeneratorFromMergeInterleaving<'a, CSS: ColorSetStorage + Sync + Send> {
    pub interleaving: &'a MergeInterleaving,
    pub coloring1: &'a CompactColexKmers<CSS>,
    pub coloring2: &'a CompactColexKmers<CSS>,
    pub merged_key_kmer_marks: &'a bitvec::vec::BitVec, // Only reporting set elements for these
    pub filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>, // With rank support
}

struct ThreadInput {
    merged_range: Range<usize>,
    s1_start_rank: usize,
    s2_start_rank: usize,
} 

impl<'a, CSS: ColorSetStorage + Sync + Send> ElementGeneratorFromMergeInterleaving<'a, CSS> {
    fn maybe_apply_filter(&self, merged_colex: usize) -> Option<usize> {
        if let Some(filter) = &self.filter {
            if !filter.get(merged_colex) {
                None // Do not report this
            } else {
                // Assign new id
                let new_id = filter.rank(merged_colex);
                Some(new_id)
            }
        } else {
            Some(merged_colex) // No filter
        }
    }
}

impl<'a, CSS: ColorSetStorage + Sync + Send> ParallelElementGenerator for ElementGeneratorFromMergeInterleaving<'a, CSS> {


    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        assert!(self.interleaving.s1.len() == self.interleaving.s2.len());
        let n = self.interleaving.s1.len();

        let s1 = &self.interleaving.s1;
        let s2 = &self.interleaving.s2;

        let thread_ranges = crate::util::segment_range(0..n, n_threads);
        let mut thread_inputs = Vec::<ThreadInput>::with_capacity(n_threads);

        let mut n_bits_s1 = 0_usize;
        let mut n_bits_s2 = 0_usize;
        for range in thread_ranges.iter() {
            let input = ThreadInput {
                merged_range: range.clone(),
                s1_start_rank: n_bits_s1,
                s2_start_rank: n_bits_s2,
            };
            thread_inputs.push(input);
            n_bits_s1 += s1[range.clone()].count_ones();
            n_bits_s2 += s2[range.clone()].count_ones();
        }

        let bar = indicatif::ProgressBar::new(n as u64);
        thread_inputs.into_par_iter().for_each(|input| {
            let mut s1_colex = input.s1_start_rank;
            let mut s2_colex = input.s2_start_rank;
            let offset_for_colors_from_2 = self.coloring1.get_set_storage().n_colors();
            for merged_colex in input.merged_range {
                if merged_colex > 0 && merged_colex % 10000 == 0 {
                    bar.inc(10000);
                }
                if self.merged_key_kmer_marks[merged_colex] {
                    if self.interleaving.s1[merged_colex] {
                        for color in self.coloring1.colex_to_set(s1_colex).iter() {
                            if let Some(new_set_id) = self.maybe_apply_filter(merged_colex) {
                                callback(SetElement{set_id: new_set_id, color});
                            }
                        }
                    }

                    if self.interleaving.s2[merged_colex] {
                        for color in self.coloring2.colex_to_set(s2_colex).iter() {
                            if let Some(new_set_id) = self.maybe_apply_filter(merged_colex) {
                                callback(SetElement{set_id: new_set_id, color: color + offset_for_colors_from_2});
                            }
                        }
                    }
                }
                s1_colex += self.interleaving.s1[merged_colex] as usize;
                s2_colex += self.interleaving.s2[merged_colex] as usize;
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

#[cfg(test)]
mod tests {
    use sbwt::{BitPackedKmerSortingMem, SeqStream, StreamingIndex, reverse_complement_in_place};
    use crate::{io::RewindableSeqStreamGenerator, set_of_sets_construction::{ParallelElementGenerator, SetElement}, util::VecVecSeqStream};
    use super::MsElementGenerator;

    struct VecColorStream {
        colors: Vec<Vec<Vec<u8>>>,
        color_idx: usize,
    }

    impl VecColorStream {
        fn new(colors: Vec<Vec<Vec<u8>>>) -> Self {
            Self { colors, color_idx: 0 }
        }
    }

    impl RewindableSeqStreamGenerator for VecColorStream {
        fn next(&mut self) -> Option<(Box<dyn SeqStream + Send + Sync>, usize)> {
            if self.color_idx == self.colors.len() {
                return None;
            }
            let seqs = self.colors[self.color_idx].clone();
            let color_idx = self.color_idx;
            self.color_idx += 1;
            Some((Box::new(VecVecSeqStream::new(seqs)), color_idx))
        }
        fn rewind(&mut self) {
            self.color_idx = 0;
        }
    }

    // Compute expected SetElements by running matching statistics directly on forward
    // and reverse-complement sequences, mirroring the WorkBatch logic.
    fn expected_elements(
        si: &StreamingIndex<'_, sbwt::SbwtIndex<sbwt::SubsetMatrix>, sbwt::LcsArray>,
        k: usize,
        color_seqs: &[Vec<Vec<u8>>],
    ) -> Vec<SetElement> {
        let mut out = Vec::new();
        for (color, seqs) in color_seqs.iter().enumerate() {
            for seq in seqs {
                let mut buf = seq.clone();
                // forward
                for (len, range) in si.matching_statistics_iter(&buf).skip(k - 1) {
                    if len == k {
                        assert_eq!(range.len(), 1);
                        out.push(SetElement { set_id: range.start, color });
                    }
                }
                // reverse complement
                reverse_complement_in_place(&mut buf);
                for (len, range) in si.matching_statistics_iter(&buf).skip(k - 1) {
                    if len == k {
                        assert_eq!(range.len(), 1);
                        out.push(SetElement { set_id: range.start, color });
                    }
                }
            }
        }
        out.sort();
        out
    }

    #[test]
    fn ms_element_generator_emits_correct_elements_including_duplicates() {
        let k = 3_usize;

        // seq0: "ACGCG"
        //   forward k-mers:  ACG, CGC, GCG
        //   rev-comp = CGCGT, k-mers: CGC, GCG, CGT
        //   → CGC and GCG each appear twice for color 0 (non-deduplication is observable)
        let seq0: &[u8] = b"ACGCG";
        let seq1: &[u8] = b"TTTGGG";

        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .algorithm(BitPackedKmerSortingMem::new())
            .k(k)
            .add_rev_comp(true)
            .build_lcs(true)
            .run_from_slices(&[seq0, seq1]);
        let lcs = lcs.unwrap();

        let color_seqs: Vec<Vec<Vec<u8>>> = vec![
            vec![seq0.to_vec()], // color 0
            vec![seq1.to_vec()], // color 1
        ];

        let si_ref = StreamingIndex::new(&sbwt, &lcs);
        let expected = expected_elements(&si_ref, k, &color_seqs);

        // Run MsElementGenerator with include_rev_comp = true.
        let gen: Box<dyn RewindableSeqStreamGenerator + Sync + Send> =
            Box::new(VecColorStream::new(color_seqs));
        let si = StreamingIndex::new(&sbwt, &lcs);
        let mut ms_gen = MsElementGenerator::new(gen, si, true);
        let got_mutex = std::sync::Mutex::new(Vec::<SetElement>::new());
        ms_gen.run(|e| got_mutex.lock().unwrap().push(e), 1);
        let mut got = got_mutex.into_inner().unwrap();
        got.sort();

        assert_eq!(got, expected);

        // CGC appears in both the forward and the reverse complement of seq0,
        // so it must be reported twice for color 0 (the non-deduplication property).
        let colex_cgc = si_ref
            .matching_statistics_iter(b"CGC")
            .skip(k - 1)
            .next()
            .map(|(_, r)| r.start)
            .expect("CGC should be in the SBWT");
        let cgc_color0 = got.iter().filter(|e| e.set_id == colex_cgc && e.color == 0).count();
        assert_eq!(cgc_color0, 2, "CGC should be reported twice for color 0 (once per strand)");
    }
}
