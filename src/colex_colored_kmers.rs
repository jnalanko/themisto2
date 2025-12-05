use bitvec::order::Lsb0;
use crossbeam::channel::{Sender, bounded};
use indicatif::ProgressStyle;
use jseqio::reverse_complement;
use jseqio::seq_db::SeqDB;
use sbwt::dbg::Dbg;
use sbwt::reverse_complement_in_place;
use sbwt::LcsArray;
use sbwt::{dbg::Node, SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::serialize::Serialize;
use simple_sds_sbwt::{ops::{BitVec, Rank}, raw_vector::AccessRaw};
use std::io::Write;
use std::ops::Range;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Acquire;
use std::sync::atomic::Ordering::Release;

use crate::atomic_bitmap::AtomicBitmap;
use crate::int_vec::CompactIntVec;
use crate::coloring_interface::{self, ColorSetStorage, ColorSetView};
use crate::index_import;
use crate::iterators::VecVecUsizeIteratorGenerator;

/// This is the main data structure in this file: a set of compressed color sets, and a mapping
/// from SBWT colex ranks to color sets such that we can look up the color set of a k-mer by its
/// colex rank in the SBWT. 
pub struct CompactColexKmers<CSS: coloring_interface::ColorSetStorage> {
    sbwt: SbwtIndex<SubsetMatrix>, 

    lcs: LcsArray,
    sets: CSS, // Distinct color sets
    map: ColexToColorSetMap, // A mapping from the colex rank of a k-mer in the SBWT into a color set id in `sets`
    color_names: Vec<String>, // User-provided names for the colors (e.g. accession numbers)
}

/// A data structure that stores the color set ids for a subset of sampled k-mers in the SBWT such that
/// the color sets of the rest can be obtained by walking forward in the de Bruijn graph to the
/// closest sampled node.
pub struct ColexToColorSetMap {
    pub sampling: simple_sds_sbwt::bit_vector::BitVector, // Marks colex ranks that have a color set stored. Has rank support.
    pub color_set_ids: CompactIntVec, // One color set id for every 1-bit in the sampling
}


struct SeqBatch {
    seqs: SeqDB,
}

impl SeqBatch {
    fn process(self, sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<SubsetMatrix>, marks: &AtomicBitmap) {
        let k = sbwt.k();
        let mut in_neighbor_buf = Vec::<(Node, u8)>::new();
        for rec in self.seqs.iter() {
            let seq = rec.seq;
            for ACGT_run in seq.split(|&c| !crate::util::IS_DNA[c as usize]) {
                if ACGT_run.len() < k { continue }

                let first = sbwt.search(&ACGT_run[0..k]).unwrap_or_else(|| {
                    panic!("k-mer {} not found in SBWT", String::from_utf8_lossy(&ACGT_run[0..k]));
                });

                assert!(first.len() == 1);

                in_neighbor_buf.clear();
                dbg.push_in_neighbors(Node{id: first.start}, &mut in_neighbor_buf);
                for (in_node, _) in in_neighbor_buf.iter() {
                    marks.set(in_node.id, true);
                }

                let last = sbwt.search(&ACGT_run[ACGT_run.len()-k..]).unwrap_or_else(|| {
                    panic!("k-mer {} not found in SBWT", String::from_utf8_lossy(&ACGT_run[ACGT_run.len()-k..]));
                });

                assert!(last.len() == 1);
                marks.set(last.start, true);
            }
        }
    }
}

// key k-mers as defined in the Themisto Bioinformatics paper:
// - Last k-mer of unitig or input sequence
// - In-neighbors of first k-mer of unitig or input sequence
// - Evenly space samples within unitigs
// IMPORTANT: currently assumes that the input `seqs` are all found in the SBWT.
// If not, we would need to search all of them and first the first and last k-mer of
// each run of matches to the index. TODO.
// Also does not mark reverse complements, so you need to provide a SeqStream that
// produces both strands.
pub fn mark_key_kmers(sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, sample_distance: usize, mut seqs: impl sbwt::SeqStream + Send, n_threads: usize) -> bitvec::vec::BitVec {

    log::info!("Initializing DBG");
    let dbg = Dbg::new(sbwt, Some(lcs), n_threads);
    let marks = AtomicBitmap::new(sbwt.n_sets());
    let dbg_ref = &dbg; // To borrow for worker threads
    let marks_ref = &marks; // To borrow for worker threads

    log::info!("Searching first and last k-mer of every input sequence");
    std::thread::scope(|scope| {

        let reader_buf_size = 1_000_000;
        let (batch_send, batch_recv) = crossbeam::channel::bounded::<SeqBatch>(4);
        let reader_handle = scope.spawn(move || {
            let progress = indicatif::ProgressBar::new_spinner();
            progress.set_style(ProgressStyle::with_template("{pos} {msg}").unwrap());
            progress.set_message(" Sequences read");
            let mut buf = SeqDB::new();
            let mut buf_total_len = 0_usize;
            while let Some(seq) = seqs.stream_next(){
                buf.push_seq(seq);
                buf_total_len += seq.len();
                if buf_total_len > reader_buf_size {
                    batch_send.send(SeqBatch{seqs: buf}).unwrap(); 
                    buf = SeqDB::new();
                    buf_total_len = 0;
                }
                progress.inc(1);
            }
            if buf_total_len > 0 { // Last batch
                batch_send.send(SeqBatch{seqs: buf}).unwrap(); 
            }
            progress.finish();
            drop(batch_send);
        });

        let mut worker_handles = Vec::new();
        for _ in 0..n_threads {
            let recv_clone = batch_recv.clone();
            let worker_handle = scope.spawn(move || {
                while let Ok(batch) = recv_clone.recv() {
                    batch.process(sbwt, dbg_ref, marks_ref);
                }
            });
            worker_handles.push(worker_handle);
        }

        reader_handle.join().unwrap();
        for h in worker_handles {
            h.join().unwrap();
        }
    });

    log::info!("Sampling along unitigs");
    dbg.iter_unitigs_with_callback(|nodes| {
        for (dist_from_end, node) in nodes.iter().rev().enumerate() {
            if dist_from_end % sample_distance == 0 {
                marks.set(node.id, true);
            }
        }
    }, n_threads);

    marks.into_bitvec()
}

impl ColexToColorSetMap {

    // sets maps from color set to its index in the distinct color sets
    #[allow(dead_code)]
    fn new(sbwt: Arc<SbwtIndex<SubsetMatrix>>, lcs: Option<&LcsArray>, sample_distance: usize, colex_to_color_set_id: Vec<usize>, n_distinct_color_sets: usize, n_threads: usize) -> Self {

        let get_colorset_fn = |colex| colex_to_color_set_id[colex]; // TODO: this actually returns a color set id. Rename here and later.
        let mut sampling_marks = Self::pick_sampled_kmers(sample_distance, &sbwt, lcs, get_colorset_fn, n_threads);

        let color_set_id_bit_width = n_distinct_color_sets.next_power_of_two().trailing_zeros() as usize;
        let mut sampled_color_set_ids = CompactIntVec::new(sampling_marks.count_ones(), color_set_id_bit_width); // In colex order
        let mut n_ids_stored = 0_usize;
        for colex in 0..sbwt.n_sets() {
            if sampling_marks.get(colex) {
                sampled_color_set_ids.set(n_ids_stored, colex_to_color_set_id[colex]);
                n_ids_stored += 1;
            }
        }

        log::info!("Building rank support for sampling marks");
        sampling_marks.enable_rank();

        Self{sampling: sampling_marks, color_set_ids: sampled_color_set_ids}
    }

    // We don't have the suffix group leader marks or the LCS array stored, so
    // we can't detect if there is an out-neighbor, but assuming there is, we
    // can find it. If the assumption does not hold, this function will still return
    // something but it will be wrong.
    // TODO: if we had a borrow to the LCS array, we could do this faster, and also
    // detect whether the node is a sink. That would mean wrapping the LCS array in
    // an Arc to avoid a self-reference in ColexColoredKmers.
    fn dbg_outneighbor_assuming_there_is_one(&self, mut colex: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> usize {
        loop {
            for char_idx in 0..sbwt.alphabet().len() {
                if sbwt.sbwt().set_contains(colex, char_idx as u8) {
                    // Found the outedge label
                    let new_colex = sbwt.lf_step(colex, char_idx);
                    return new_colex;
                }
            }
            // No outedges found -> colex is not a suffix group leader position
            assert!(colex > 0);
            colex -= 1;
        }
    }

    // Returns the colex rank of the next sampled node from here, and
    // the number of nodes walked (0 means the starting point was already marked).
    fn walk_to_next_sample(&self, mut colex: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> (usize, usize) {
        let mut depth = 0;
        loop {
            if self.sampling.get(colex) {
                return (colex, depth)
            } else {
                // This set is not stored -> walk forward in the de Bruijn graph
                // Since this is not sampled,  this can not be a sink.
                colex = self.dbg_outneighbor_assuming_there_is_one(colex, sbwt);
                depth += 1;
            }
        }
    }

    fn colex_to_color_set_id(&self, colex: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> usize {
        let pos = self.walk_to_next_sample(colex, sbwt).0;
        self.color_set_ids.get(self.sampling.rank(pos))
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        self.sampling.serialize(out).unwrap();
        self.color_set_ids.serialize(out);
    }

    pub fn load(input: &mut impl std::io::Read) -> Self {
        let sampling = simple_sds_sbwt::bit_vector::BitVector::load(input).unwrap();
        let color_set_ids = CompactIntVec::load(input);

        assert_eq!(color_set_ids.len(), sampling.count_ones());

        Self{sampling, color_set_ids}
    }

    /// Utility function used in construction
    fn pick_sampled_kmers<F: Fn(usize) -> usize + Sync + Send>(sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>, lcs: Option<&LcsArray>, get_colorset_fn: F, n_threads: usize) -> simple_sds_sbwt::bit_vector::BitVector {
        // Find starts of unitigs. Walk forward to the end of the unitig. Segment by color sets.
        
        let marks = simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
        let marks_mutex = Mutex::new(marks); // Need thread-safe modifications
        let marks_mutex_borrow = &marks_mutex; // Passed into the callback

        let callback = |nodes: &[Node]| {
            let mut marks = marks_mutex_borrow.lock().unwrap();

            let mut prev_set: Option<usize> = None;
            let mut prev_sample_pos = usize::MAX;
            for (v_pos, v) in nodes.iter().enumerate().rev() {
                let colex = v.id; 
                let cur_set = get_colorset_fn(colex);

                // Sample this node if any of the following are true:
                // - v is the last node of the unitig
                // - v has a different color set than the previous node in iteration order 
                // - v is far enough from the previous sampled node 
                if prev_set.is_none() || cur_set != prev_set.unwrap() || prev_sample_pos - v_pos >= sample_distance {
                    marks.set_bit(colex, true);
                    prev_sample_pos = v_pos;
                }
                prev_set = Some(cur_set);
            }
        };

        log::info!("Initializing the de Bruijn graph");
        let dbg = sbwt::dbg::Dbg::new(sbwt, lcs, n_threads);

        log::info!("Iterating unitigs");
        dbg.iter_unitigs_with_callback(callback, n_threads);

        let marks = marks_mutex.into_inner().unwrap();
        let marks = simple_sds_sbwt::bit_vector::BitVector::from(marks);

        let n_sampled = marks.count_ones();
        log::info!("Sampled {} out of {} k-mers ({:.2}%)", n_sampled, sbwt.n_kmers(), n_sampled as f64 / sbwt.n_kmers() as f64 * 100.0);

        log::info!("Unitig sampling finished");

        marks
    }
}

struct UnitigImportSeqBatch {
    pub concat: Vec<u8>,
    pub starts: Vec<usize>, // Has concat.len() at the end 
    pub color_set_ids: Vec<usize>, // The color set id for each sequence in the concatenation
}

impl UnitigImportSeqBatch {
    fn get_seq_mut(&mut self, idx: usize) -> &mut [u8] {
        &mut self.concat[self.starts[idx]..self.starts[idx+1]]
    }

    fn n_seqs(&self) -> usize {
        self.starts.len() - 1 // Has concat.len() at the end
    }

    fn process(mut self, results_out: &mut Vec<(usize, usize)>, index: &sbwt::StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>, sample_distance: usize) { // Todo this should consume that batch since it's edited
        let k = index.k();

        let mut set_ids_fn = |seq: &[u8], color_set_id| {
            if seq.len() < k { return }
            let mut distance_from_end = seq.len()-k+1;
            for (start, (match_len, colex_range)) in index.matching_statistics_iter(seq).skip(k-1).enumerate() {
                distance_from_end -= 1;
                if match_len == k {
                    assert_eq!(colex_range.len(), 1);
                    if distance_from_end % (sample_distance-1) == 0 {
                        results_out.push((colex_range.start, color_set_id));
                    }
                } else {
                    panic!( // TODO: return error instead
                        "Error reading unitigs from dump: k-mer {} not found in SBWT", 
                        String::from_utf8_lossy(&seq[start..start+k])
                    );
                }
            }
            assert!(distance_from_end == 0);
        };

        for seq_idx in 0..self.n_seqs() {
            let color_set_id = self.color_set_ids[seq_idx];
            let seq = self.get_seq_mut(seq_idx);

            // Process both forward and reverse complement directions
            set_ids_fn(seq, color_set_id);
            reverse_complement_in_place(seq); // TODO could revcomp all at once
            set_ids_fn(seq, color_set_id);
        }

    }
}

// Todo: the logic associated with this is all over the place. Refactor into one place.
fn unitig_import_parser_thread(unitig_dump: impl std::io::BufRead + Send + 'static, buf_cap: usize, out: Sender<UnitigImportSeqBatch>){
        
    let mut seqs = jseqio::reader::DynamicFastXReader::new(unitig_dump).unwrap();

    let mut cur_concat = Vec::<u8>::with_capacity(buf_cap);
    let mut cur_starts = Vec::<usize>::new();
    let mut cur_color_set_ids = Vec::<usize>::new();
    
    while let Some(rec) = seqs.read_next().unwrap() {
        let color_set_id = index_import::get_color_set_id_from_fasta_header(rec.head);
        cur_color_set_ids.push(color_set_id);

        // Add to concatenation
        cur_starts.push(cur_concat.len());
        cur_concat.extend(rec.seq);

        if cur_concat.len() >= buf_cap {
            cur_starts.push(cur_concat.len()); // End sentinel, as required
            let batch = UnitigImportSeqBatch{concat: cur_concat, starts: cur_starts, color_set_ids: cur_color_set_ids};
            out.send(batch).unwrap();

            // Start a new batch
            cur_concat = Vec::<u8>::with_capacity(buf_cap);
            cur_starts = Vec::<usize>::new();
            cur_color_set_ids = Vec::<usize>::new();
        }
    }

    if !cur_concat.is_empty() {
        // Send remaining batch
        cur_starts.push(cur_concat.len()); // End sentinel, as required
        let batch = UnitigImportSeqBatch{concat: cur_concat, starts: cur_starts, color_set_ids: cur_color_set_ids};
        out.send(batch).unwrap();
    }

    log::info!("Producer thread: all work pushed to work queue");
    drop(out);
}

impl<CSS: ColorSetStorage> CompactColexKmers<CSS> {

    pub fn new(sbwt: SbwtIndex<SubsetMatrix>, lcs: LcsArray, colex_map: ColexToColorSetMap, color_sets: CSS, color_names: Option<&[String]>)
    -> CompactColexKmers<CSS> {
        let color_names = if let Some(names) = color_names {
            assert!(names.len() == color_sets.n_colors());
            names.to_vec()
        } else {
            // Assign default color names
            (0..color_sets.n_colors()).map(|x| format!("color_{}", x)).collect::<Vec<String>>()
        };
        Self {sbwt, lcs, sets: color_sets, map: colex_map, color_names}
    }

    pub fn into_parts(self) -> (SbwtIndex<SubsetMatrix>, LcsArray, ColexToColorSetMap, CSS, Vec<String>) {
        (self.sbwt, self.lcs, self.map, self.sets, self.color_names)
    }

    pub fn sbwt(&self) -> &SbwtIndex<SubsetMatrix> {
        &self.sbwt
    }

    pub fn lcs(&self) -> &LcsArray {
        &self.lcs
    }

    #[allow(dead_code)]
    pub fn get_map(&self) -> &ColexToColorSetMap {
        &self.map
    }

    pub fn new_from_colored_unitig_dump(
        sbwt: SbwtIndex<SubsetMatrix>, 
        lcs: LcsArray, 
        sample_distance: usize,
        n_threads: usize,
        metadata_dump: impl std::io::BufRead, 
        unitig_dump: impl std::io::BufRead + Send + 'static, 
        color_dump: impl std::io::BufRead) 
        -> Self {

        assert!(sample_distance > 0);

        log::info!("Reading metadata");
        let metadata = index_import::read_index_dump_metadata(metadata_dump);

        log::info!("Building (colex, color set id) pairs");
        let mut colex_to_color_set_id = std::thread::scope(|scope| {
            let (parser_out, worker_in) = bounded(n_threads);

            // Create producer
            let producer_handle = scope.spawn(move || {
                unitig_import_parser_thread(unitig_dump, 1 << 20, parser_out);
            });

            // Create workers
            let mut worker_handles = Vec::<std::thread::ScopedJoinHandle::<Vec::<_>>>::new();

            for _ in 0..n_threads {
                let worker_in_clone = worker_in.clone();
                let index = sbwt::StreamingIndex::new(&sbwt, &lcs);
                worker_handles.push(scope.spawn(move || {
                    let mut our_colex_to_color_set_id: Vec<(usize, usize)> = vec![]; // (colex, coled_set_id)
                    while let Ok(batch) = worker_in_clone.recv(){
                        batch.process(&mut our_colex_to_color_set_id, &index, sample_distance);
                    }
                    our_colex_to_color_set_id
                }));
            }

            producer_handle.join().unwrap(); // Wait for the producer to finish

            let mut colex_to_color_set_id = Vec::<(usize, usize)>::new(); // Collect thread outputs here
            for h in worker_handles { // Wait for the workers to finish
                colex_to_color_set_id.extend(h.join().unwrap());
            }
            colex_to_color_set_id
        });

        let n_colors = metadata.num_colors;

        log::info!("Sorting (colex, color set id) pairs");
        colex_to_color_set_id.sort(); // Sorts by colex
        let bit_width = metadata.num_color_sets.next_power_of_two().trailing_zeros() as usize;
        log::info!("Building compressed representation for color set ids");
        let mut stored_color_set_ids = CompactIntVec::new(colex_to_color_set_id.len(), bit_width);
        let mut sample_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
        for (rank, (colex, id)) in colex_to_color_set_id.into_iter().enumerate() {
            stored_color_set_ids.set(rank, id);
            sample_marks.set_bit(colex, true);
        }
        let mut sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(sample_marks);
        sample_marks.enable_rank();
        log::info!("Marked {:.2} % of all k-mers", sample_marks.count_ones() as f64 / sbwt.n_kmers() as f64 * 100.0);

        log::info!("Reading distinct color sets");
        let color_set_stream = index_import::ColorSetDumpIterGenerator::new(color_dump);
        let distinct_css = CSS::new(color_set_stream, n_colors);
        let distinct_css = *distinct_css; // Unbox

        let colex_map = ColexToColorSetMap {
            sampling: sample_marks,
            color_set_ids: stored_color_set_ids,
        };

        let color_names: Vec<String> = (0..distinct_css.n_colors()).map(|x| x.to_string()).collect();
        Self {sbwt, lcs, sets: distinct_css, map: colex_map, color_names}
    }

    #[allow(dead_code)] // Might still be useful in the future
    pub fn new_single_colored(sbwt: SbwtIndex<SubsetMatrix>, lcs: LcsArray, sample_distance: usize, n_threads: usize, color_name: String) -> Self {
        let n_colors = 1;
        let int_bitwidth = 1;

        // The only color set is {0}
        let vv = VecVecUsizeIteratorGenerator{sets: vec![vec![0]], pos: 0};
        let sets = CSS::new(vv, n_colors);

        log::info!("Sampling nodes");
        let mut unitig_samples = ColexToColorSetMap::pick_sampled_kmers(sample_distance, &sbwt, Some(&lcs), |_colex| 0, n_threads);
        unitig_samples.enable_rank();
        log::info!("Storing color set ids for sampled nodes");
        //let color_set_ids = IntVector::with_len(unitig_samples.count_ones(), int_bitwidth, 0).unwrap();
        let color_set_ids = CompactIntVec::new(unitig_samples.count_ones(), int_bitwidth);
        let colex_map = ColexToColorSetMap{
            sampling: unitig_samples,
            color_set_ids,
        };
        Self {sbwt, lcs, sets: *sets, map: colex_map, color_names: vec![color_name]}
    }

    pub fn colex_to_set_id(&self, colex: usize) -> usize {
        self.map.colex_to_color_set_id(colex, &self.sbwt)
    }

    pub fn set_id_to_set<'a>(&'a self, id: usize) -> CSS::SetView<'a> {
        self.sets.get_set_view(id)
    }

    pub fn colex_to_set<'a>(&'a self, colex: usize) -> CSS::SetView<'a> {
        self.sets.get_set_view(self.colex_to_set_id(colex))
    }

    fn serialization_magic_string() -> &'static [u8; 4] {
        b"CCKM" // Compact Colex Kmers
    }

    fn serialization_version() -> u64 {
        1_u64
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        out.write_all(Self::serialization_magic_string()).unwrap();
        out.write_all(&Self::serialization_version().to_le_bytes()).unwrap();

        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();
        self.sets.serialize(out);
        self.map.serialize(out);
        
        // Serialize color names: first the number of names, then each name length and name
        let n_names = self.color_names.len() as u64;
        out.write_all(&n_names.to_le_bytes()).unwrap();
        for name in self.color_names.iter() {
            let name_bytes = name.as_bytes();
            let name_len = name_bytes.len() as u64;
            out.write_all(&name_len.to_le_bytes()).unwrap();
            out.write_all(name_bytes).unwrap();
        }
    }

    /// If this struct is going to be merged with [crate::coloring::merge_colorings], it will
    /// need select support on the sbwt. We need to build it already during loading because
    /// once the sbwt is put on to the heap into an Arc, it cannot be modified anymore.
    /// Unless we make it an Arc<Refcell<...>>, but that might have overhead because then
    /// it will do run-time borrow checking on every access if I understand correctly.
    pub fn load(input: &mut impl std::io::Read, enable_select: bool) -> Self {
        // Read and verify magic string
        let mut magic_string = [0_u8; 4];
        input.read_exact(&mut magic_string).unwrap();
        if magic_string != *Self::serialization_magic_string() {
            panic!("Error loading CompactColexKmers: expected bytes {:?} but found {:?}", Self::serialization_magic_string(), magic_string);
        }
        let mut version_bytes = [0_u8; 8];
        input.read_exact(&mut version_bytes).unwrap();
        let version = u64::from_le_bytes(version_bytes);
        if version != Self::serialization_version() {
            panic!("Error loading CompactColexKmers: expected version {} but found {}", Self::serialization_version(), version);
        }

        let mut sbwt = SbwtIndex::<SubsetMatrix>::load(input).unwrap();
        if enable_select {
            log::info!("Building select support");
            sbwt.build_select();
        }
        let lcs = LcsArray::load(input).unwrap();
        let sets = CSS::load(input);
        let map = ColexToColorSetMap::load(input);
        assert_eq!(map.sampling.len(), sbwt.n_sets());

        // Load color names
        let mut n_names_bytes = [0_u8; 8];
        input.read_exact(&mut n_names_bytes).unwrap();
        let n_names = u64::from_le_bytes(n_names_bytes) as usize;
        assert_eq!(n_names, sets.n_colors());
        let mut color_names = Vec::<String>::with_capacity(n_names);
        for _ in 0..n_names {
            let mut name_len_bytes = [0_u8; 8];
            input.read_exact(&mut name_len_bytes).unwrap();
            let name_len = u64::from_le_bytes(name_len_bytes) as usize;
            let mut name_bytes = vec![0_u8; name_len];
            input.read_exact(&mut name_bytes).unwrap();
            let name = String::from_utf8(name_bytes).unwrap();
            color_names.push(name);
        }

        CompactColexKmers{sbwt, lcs, sets, map, color_names}
    }

    pub fn lookup_kmer_color_sets(&self, seq: &[u8]) -> Vec<Option<CSS::SetView<'_>>> {
        let mut buffer = Vec::<Option<usize>>::new();
        self.push_color_set_ids_to_buffer(seq, &mut buffer);
        buffer.into_iter().map(|opt| opt.map(|x| self.set_id_to_set(x))).collect()
    }

    // Does not clear the buffer
    pub fn push_color_set_ids_to_buffer(&self, seq: &[u8], buffer: &mut Vec<Option<usize>>) {
        let k = self.sbwt.k();
        if seq.len() < k {
            return;
        }

        let si = sbwt::StreamingIndex::new(&self.sbwt, &self.lcs);

        #[derive(Eq, PartialEq, Debug)]
        enum ColexPos {
            Sampled(usize),
            SameAsNext(usize),
            AbsentKmer
        }

        let mut colex_positions = Vec::<ColexPos>::with_capacity(seq.len()-k+1);

        // Pass 1: Compute colex ranks and whether the k-mers are sampled
        for (len, range) in si.matching_statistics_iter(seq).skip(k-1) {
            if len == k {
                assert!(range.len() == 1);
                let colex = range.start;
                if self.map.sampling.get(colex) {
                    colex_positions.push(ColexPos::Sampled(colex));
                } else {
                    colex_positions.push(ColexPos::SameAsNext(colex));
                }
            } else {
                colex_positions.push(ColexPos::AbsentKmer);
            }
        }

        // Pass 2: look up set ids right-to-left, reusing set ids for the next k-mer whenever possible
        let old_buf_len = buffer.len();
        buffer.resize(old_buf_len + colex_positions.len(), None); // Make space for results
        let set_ids_output = &mut buffer[old_buf_len..];
        for i in (0..colex_positions.len()).rev() {
            match colex_positions[i] {
                ColexPos::Sampled(colex) => {
                    let set_id = self.map.colex_to_color_set_id(colex, &self.sbwt);
                    set_ids_output[i] = Some(set_id);
                },
                ColexPos::SameAsNext(colex) => {
                    if i+1 == set_ids_output.len() || set_ids_output[i+1].is_none() {
                        // Can not copy color set id from position i+1
                        let set_id = self.map.colex_to_color_set_id(colex, &self.sbwt);
                        set_ids_output[i] = Some(set_id);
                    } else {
                        set_ids_output[i] = set_ids_output[i+1];
                    }
                },
                ColexPos::AbsentKmer => {
                    set_ids_output[i] = None;
                },
            }
        }

        // Debug verification:
        /*
        for (i, kmer) in seq.windows(k).enumerate() {
            match self.sbwt.search(kmer) {
                Some(colex) => {
                    assert!(set_ids_output[i] == Some(self.map.colex_to_color_set_id(colex.start, &self.sbwt))); 
                },
                None => {
                    assert!(set_ids_output[i].is_none());    
                },
            }
        }
        */
    }

    pub fn get_k(&self) -> usize {
        self.sbwt.k()
    }

    pub fn get_set_storage(&self) -> &CSS {
        &self.sets
    }

    pub fn get_color_names(&self) -> &Vec<String> {
        &self.color_names
    }

    pub fn break_to_colored_subunitigs(&self, unitig_colex_ranks: &[usize], _unitig_string: &[u8]) -> (Vec<usize>, Vec<Range<usize>>){
        if unitig_colex_ranks.len() == 0 {
            // Make this a special case to ensure that there is always at least
            // one run to avoid a special case at the end.
            return (vec![], vec![]);
        }
        let mut subunitig_color_set_ids: Vec<usize> = vec![];
        let mut subunitigs: Vec<Range<usize>> = vec![]; // Ranges of k-mers (= starts of k-mers)
        let mut current_run_set_id = usize::MAX; // Will be set at the start of the first iteration
        let mut current_run_end = unitig_colex_ranks.len(); 

        // Iterate from end to start, updating the color set when the current
        // node is marked.
        for (pos, &colex) in unitig_colex_ranks.iter().enumerate().rev() {
            if pos == unitig_colex_ranks.len()-1 {
                assert!(self.map.sampling.get(colex)); // Last position of a unitig should always be marked
            }

            if self.map.sampling.get(colex) {
                // Update the set id
                let new_set_id = self.colex_to_set_id(colex);

                if new_set_id != current_run_set_id {
                    // Close the active run (if exists)
                    let start = pos + 1;
                    if current_run_end > start { // Active run exists
                        subunitigs.push(start..current_run_end);
                        subunitig_color_set_ids.push(current_run_set_id);
                        current_run_end = pos + 1;
                    }
                }
                current_run_set_id = new_set_id;
            }
        }

        // Close the active run (exists because of the assert at the start)
        assert!(current_run_set_id != usize::MAX);
        subunitigs.push(0..current_run_end);
        subunitig_color_set_ids.push(current_run_set_id);

        subunitigs.reverse();
        subunitig_color_set_ids.reverse();

        (subunitig_color_set_ids, subunitigs)
    }

    #[allow(clippy::type_complexity)] // Yeah yeah I know
    fn search_unitig_from(&self, v: Node, dbg: &Dbg<'_, SubsetMatrix>) -> (Vec<usize>, Vec<usize>, Vec<u8>, Vec<Range<usize>>, Vec<usize>) {
        // Walk the unitig in forward orientation, and then backwards
        let k = self.sbwt.k();
        let mut workspace = Vec::<u8>::new();
        let nodes = dbg.walk_unitig_from(v, &mut workspace);
        workspace.clear();
        let mut unitig_string = Vec::<u8>::new();
        dbg.push_unitig_string(&nodes, &mut unitig_string);

        let string_len = unitig_string.len();
        assert!(string_len >= k);
        let last_kmer = &unitig_string[string_len-k..];
        let last_kmer_rc = reverse_complement(last_kmer);
        let last_kmer_rc_colex = self.sbwt.search(&last_kmer_rc).unwrap_or_else(|| panic!(
            "Reverse complement of k-mer {} not found in index", 
            String::from_utf8_lossy(last_kmer))
        ).start;
        let rc_nodes = dbg.walk_unitig_from(sbwt::dbg::Node{id: last_kmer_rc_colex}, &mut workspace);

        let fw_colex: Vec<usize> = nodes.into_iter().map(|v| v.id).collect();
        let rc_colex: Vec<usize> = rc_nodes.into_iter().rev().map(|v| v.id).collect();
        assert_eq!(fw_colex.len(), rc_colex.len());

        // Figure out color set id runs in the forward strand 
        let (subuniting_color_set_ids, subunitig_kmer_ranges) = self.break_to_colored_subunitigs(&fw_colex, &unitig_string);

        (fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subuniting_color_set_ids)
    }

    #[allow(clippy::too_many_arguments)] // Yeah yeah I know
    fn visit_and_output_kmers(&self, unitig_string: &[u8], subunitig_kmer_ranges: &[Range<usize>], subunitig_color_set_ids: &[usize], fw_colex: &[usize], rc_colex: &[usize], visited: &mut bitvec::vec::BitVec, unitigs_out: &mut impl Write, unitig_id: &mut usize) {

        let k = self.sbwt.k();

        for (subunitig_idx, r) in subunitig_kmer_ranges.iter().enumerate() {
            // All k-mers in this subunitig have the same color set id.
            // It would be nice if we could just figure out the unvisited
            // runs of k-mers and visit and output those, but there is a subtle problem:
            // A subunitig may loop back to itself in reverse complement orientation.
            // Printing the subunitig would print the same k-mer in both orientations.
            // So, we need to keep track of the visited bit vector also while processing
            // a subunitig, and end the subunitig when we encounter a visited k-mer.
            let subunitig = &unitig_string[r.start..r.end+k-1];
            let color_set_id = subunitig_color_set_ids[subunitig_idx];
            let fw_colex_slice = &fw_colex[r.start..r.end];
            let rc_colex_slice = &rc_colex[r.start..r.end];

            let mut subsubunitig_start: Option<usize> = None;
            for kmer_idx in 0..fw_colex_slice.len() {
                if !visited[fw_colex_slice[kmer_idx]] {
                    // Extend the current subunitig and visit this k-mer
                    if subsubunitig_start.is_none() {
                        subsubunitig_start = Some(kmer_idx);
                    }
                    visited.set(fw_colex_slice[kmer_idx], true);
                    visited.set(rc_colex_slice[kmer_idx], true);
                } else {
                    // Already visited! Output the current subunitig
                    if let Some(s) = subsubunitig_start {
                        let e = kmer_idx + k - 1;
                        writeln!(unitigs_out, "> unitig_id={} color_set_id={}", unitig_id, color_set_id).unwrap();
                        unitigs_out.write_all(&subunitig[s..e]).unwrap();
                        unitigs_out.write_all(b"\n").unwrap();
                        *unitig_id += 1;
                    }
                    subsubunitig_start = None;
                }
            }

            // Write the last subunitig if it's still open
            if let Some(s) = subsubunitig_start {
                let e = fw_colex_slice.len() + k - 1;
                writeln!(unitigs_out, "> unitig_id={} color_set_id={}", unitig_id, color_set_id).unwrap();
                unitigs_out.write_all(&subunitig[s..e]).unwrap();
                unitigs_out.write_all(b"\n").unwrap();
                *unitig_id += 1;
            }
        }
    }

    /// Canonical here means whichever strand is visited first.
    /// This assumes that the color set of a forward k-mer and a reverse k-mer is the same.
    /// Returns the number of unitigs written
    fn export_canonical_unitigs_with_shared_color_set(&self, mut unitigs_out: impl Write + Sync + Send, n_threads: usize) -> usize where CSS : Sync {

        log::info!("Initializing the de Bruijn graph");
        let dbg = Dbg::new(&self.sbwt, Some(&self.lcs), n_threads);
        let dbg_ref = &dbg;

        let k = self.get_k();

        log::info!("Computing unitigs");
        let n_unitig_searches = std::sync::atomic::AtomicUsize::new(0);
        let n_unitig_searches_ref = &n_unitig_searches;

        let bar = indicatif::ProgressBar::new(self.sbwt.n_sets() as u64);
        let n_unitigs = std::thread::scope(|scope| {

            // Channels of tuples of with these fields: 
            //   * forward colex ranks 
            //   * reverse complement colex ranks
            //   * unitig string
            //   * colored subunitig k-mer ranges
            //   * color set ids of the colored subunitig ranges
            // TODO: less heap allocation
            let (worker_out, collector_in) = bounded::<(Vec<usize>, Vec<usize>, Vec<u8>, Vec<Range<usize>>, Vec<usize>)>(n_threads);

            // Create unitig search threads 
            let mut worker_handles = Vec::<_>::new();
            let bar_ref = &bar;
            for thread_id in 0..n_threads { 
                let worker_out_clone = worker_out.clone();
                let handle = scope.spawn(move || {
                    // Iterating all colex positions that have remainder thread_id modulo number of threads
                    let mut colex = thread_id;
                    while colex < self.sbwt.n_sets() {
                        let v = Node { id: colex };
                        if !dbg_ref.is_dummy_colex_position(colex) && dbg_ref.is_first_kmer_of_unitig(v) {
                            n_unitig_searches_ref.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            worker_out_clone.send(self.search_unitig_from(v, dbg_ref)).unwrap();
                        }
                        colex += n_threads;
                        if ((colex - thread_id)/n_threads) % 10000 == 0 {
                            bar_ref.inc(10000);
                            //eprintln!("number of unitig searches: {}", n_unitig_searches_ref.load(std::sync::atomic::Ordering::Relaxed));
                        }
                    }
                    log::info!("Thread {} finished", thread_id);
                });
                worker_handles.push(handle);
            }

            let collector_handle = scope.spawn(move || {
                // We maintain the visited bit vector so that when we mark a k-mer, we also mark its
                // reverse complement.
                let mut unitig_id = 0_usize;

                // Bitvector marking visited colex ranks 
                let mut visited = bitvec::bitvec![usize, Lsb0; 0; self.sbwt.n_sets()];

                // Process all non-cyclic unitigs
                while let Ok((fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subunitig_color_set_ids)) = collector_in.recv() {
                    self.visit_and_output_kmers(&unitig_string, &subunitig_kmer_ranges, &subunitig_color_set_ids, &fw_colex, &rc_colex, &mut visited, &mut unitigs_out, &mut unitig_id); 
                }

                // Process remaining cyclic unitigs
                log::info!("Processing remaining cyclic unitigs");
                let n_acyclic = unitig_id; // This many unitigs have been written so far
                let mut colex = 0_usize;
                while colex < visited.len() {
                    colex = match visited[colex..].first_zero() {
                        Some(i) => colex + i,
                        None => break,
                    };
                    if !dbg_ref.is_dummy_colex_position(colex) {
                        let (fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subunitig_color_set_ids)
                        = self.search_unitig_from(Node { id: colex }, dbg_ref);

                        // Make sure it's really cyclic
                        assert!(unitig_string[..k-1] == unitig_string[unitig_string.len()-(k-1)..]);

                        self.visit_and_output_kmers(&unitig_string, &subunitig_kmer_ranges, &subunitig_color_set_ids, &fw_colex, &rc_colex, &mut visited, &mut unitigs_out, &mut unitig_id);
                    }
                    colex += 1;
                }
                unitigs_out.flush().unwrap();
                log::info!("Found {} cyclic unitigs", unitig_id - n_acyclic);
                unitig_id
            });

            for h in worker_handles { // Wait for the workers to finish
                h.join().unwrap();
            }

            drop(worker_out);

            // Wait for the collector to finish
            let n_unitigs = collector_handle.join().unwrap();

            #[allow(clippy::let_and_return)] // It's renaming of the variable. Clearer this way.
            n_unitigs
        });
        bar.finish();

        log::info!("Wrote {} unitigs", n_unitigs);
        n_unitigs
    }


    /// Same format as [crate::index_import].
    /// Select support must be built before calling this!
    pub fn export_colored_unitigs(&self, mut metadata_out: impl Write + Sync + Send, unitigs_out: impl Write + Sync + Send, mut colors_out: impl Write + Sync + Send, n_threads: usize)
        where CSS: Sync {
        // The metadata should look like this:
        // num_colors=3682
        // num_unitigs=9314735
        // num_color_sets=5591009
        // k=31

        log::info!("Exporting to colored unitigs");

        let n_unitigs = self.export_canonical_unitigs_with_shared_color_set(unitigs_out, n_threads);

        // Write color sets
        // Lines should look like this:
        // color_set_id=9 size=7 3 4 9 12 14 15 16
        for set_id in 0..self.sets.n_sets() {
            let set_view = self.sets.get_set_view(set_id);
            write!(colors_out, "color_set_id={} size={}", set_id, set_view.len()).unwrap();
            for color in set_view.iter() { 
                write!(colors_out, " {}", color).unwrap(); // TODO: faster IO
            }
            writeln!(colors_out).unwrap();
        }

        metadata_out.write_all(format!("num_colors={}\n", self.sets.n_colors()).as_bytes()).unwrap();
        metadata_out.write_all(format!("num_unitigs={}\n", n_unitigs).as_bytes()).unwrap();
        metadata_out.write_all(format!("num_color_sets={}\n", self.sets.n_sets()).as_bytes()).unwrap();
        metadata_out.write_all(format!("k={}\n", self.sbwt.k()).as_bytes()).unwrap();
    }

    pub fn compute_index_stats(&self, n_threads: usize) -> IndexStats where CSS: Sync {
        log::info!("Computing size of distinct color sets");
        let total_size_of_distinct_color_sets: usize = (0..self.sets.n_sets()).map(|i| self.sets.get_set_view(i).len()).sum();

        log::info!("Initializing de Bruijn graph");
        let dbg = Dbg::new(&self.sbwt, Some(&self.lcs), n_threads);
        let total_kmer_color_set_size = AtomicUsize::new(0);
        let n_unitigs = AtomicUsize::new(0);
        let total_unitig_length = AtomicUsize::new(0);
        let min_unitig_length = AtomicUsize::new(usize::MAX);
        let max_unitig_length = AtomicUsize::new(0);

        let bar = indicatif::ProgressBar::new(self.sbwt.n_kmers() as u64);
        dbg.iter_unitigs_with_callback(|nodes|{
            let mut cur_color_set_len: Option<usize> = None;
            for v in nodes.iter().rev() {
                if self.map.sampling.get(v.id) { 
                    let id = self.colex_to_set_id(v.id);
                    cur_color_set_len = Some(self.set_id_to_set(id).len());
                } else {
                    // Same set as previous
                    assert!(cur_color_set_len.is_some()); // End of unitig should always be sampled
                }
                total_kmer_color_set_size.fetch_add(cur_color_set_len.unwrap(), Release);
            }

            let unitig_length = nodes.len()+self.sbwt().k()-1;
            total_unitig_length.fetch_add(unitig_length, Release);
            min_unitig_length.fetch_min(unitig_length, Release);
            max_unitig_length.fetch_max(unitig_length, Release);
            n_unitigs.fetch_add(1, Release);

            bar.inc(nodes.len() as u64);
        }, n_threads);
        bar.finish();

        let total_kmer_color_set_size = total_kmer_color_set_size.load(Acquire);
        let n_sampled_kmers = self.map.sampling.count_ones();

        IndexStats {
            n_colors: self.sets.n_colors(),
            n_kmers: self.sbwt.n_kmers(),
            n_sbwt_sets: self.sbwt.n_sets(),
            n_distinct_color_sets: self.sets.n_sets(),
            n_sampled_kmers,
            total_size_of_distinct_color_sets,
            total_kmer_color_set_size,
            n_unitigs: n_unitigs.load(Acquire),
            total_unitig_length: total_unitig_length.load(Acquire),
            max_unitig_length: max_unitig_length.load(Acquire),
            min_unitig_length: min_unitig_length.load(Acquire),
        }
    }
}

pub struct IndexStats {
    pub n_colors: usize,
    pub n_kmers: usize,
    pub n_sbwt_sets: usize,
    pub n_distinct_color_sets: usize,
    pub n_sampled_kmers: usize,
    pub total_size_of_distinct_color_sets: usize,
    pub total_kmer_color_set_size: usize,
    pub n_unitigs: usize,
    pub total_unitig_length: usize,
    pub max_unitig_length: usize,
    pub min_unitig_length: usize,
}

impl IndexStats {
    pub fn mean_size_of_distinct_sets(&self) -> f64 {
        self.total_size_of_distinct_color_sets as f64 / self.n_distinct_color_sets as f64
    }

    pub fn mean_kmer_color_set_size(&self) -> f64 {
        self.total_kmer_color_set_size as f64 / self.n_kmers as f64
    }

    pub fn sample_fraction(&self) -> f64 {
        self.n_sampled_kmers as f64 / self.n_kmers as f64
    }

    pub fn mean_unitig_length(&self) -> f64 {
        self.total_unitig_length as f64 / self.n_unitigs as f64
    }
}
