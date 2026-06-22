use bitvec::order::Lsb0;
use crossbeam::channel::{Sender, bounded};
use indicatif::ProgressStyle;
use jseqio::reverse_complement;
use jseqio::seq_db::SeqDB;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;
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
use crate::io::RewindableSeqStreamGenerator;
use crate::iterators::VecVecUsizeIteratorGenerator;
use crate::work_dispatcher;

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

struct FirstLastKmerWorker<'a> {
    in_neighbor_buf: Vec::<(Node, u8)>,
    rev_comp_buf: Vec<u8>,
    k: usize,
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    dbg: &'a Dbg<'a, SubsetMatrix>,
    marks: &'a AtomicBitmap,
    add_rev_comp: bool,
}

impl<'a> crate::work_dispatcher::Worker for FirstLastKmerWorker<'a> {
    fn process(&mut self, seq: &[u8], _color: usize) {
        self.process_internal(seq, false);
        if self.add_rev_comp {
            self.process_internal(seq, true);
        }
    }
}

impl<'a> FirstLastKmerWorker<'a> {
    // This is a bit ugly. We process seq, unless rev_comp_instead is true,
    // in which case we compute from the reverse complement. It's done this
    // way because it was hard to make the borrow checker happy otherwise
    // because we need to mutate the reverse complement buffer here, so
    // this function needs a mutable &self, but it can't simultaneously
    // take an immutable reference into its own reverse complement buffer.
    fn process_internal(&mut self, seq: &[u8], rev_comp_instead: bool) {

        let seq = if rev_comp_instead {
            self.rev_comp_buf.clear();
            self.rev_comp_buf.extend_from_slice(seq);
            reverse_complement_in_place(&mut self.rev_comp_buf);
            &self.rev_comp_buf
        } else { 
            seq 
        };

        for ACGT_run in seq.split(|&c| !crate::util::IS_DNA[c as usize]) {
            if ACGT_run.len() < self.k { continue }

            let first = self.sbwt.search(&ACGT_run[0..self.k]).unwrap_or_else(|| {
                panic!("k-mer {} not found in SBWT", String::from_utf8_lossy(&ACGT_run[0..self.k]));
            });

            assert!(first.len() == 1);

            self.in_neighbor_buf.clear();
            self.dbg.push_in_neighbors(Node{id: first.start}, &mut self.in_neighbor_buf);
            for (in_node, _) in self.in_neighbor_buf.iter() {
                self.marks.set(in_node.id, true);
            }

            let last = self.sbwt.search(&ACGT_run[ACGT_run.len()-self.k..]).unwrap_or_else(|| {
                panic!("k-mer {} not found in SBWT", String::from_utf8_lossy(&ACGT_run[ACGT_run.len()-self.k..]));
            });

            assert!(last.len() == 1);
            self.marks.set(last.start, true);
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
pub fn mark_key_kmers(sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, sample_distance: usize, seq_stream_gen: &mut Box<dyn RewindableSeqStreamGenerator + Sync + Send>, n_threads: usize, n_parser_threads: usize, add_rev_comp: bool) -> bitvec::vec::BitVec {

    log::info!("Initializing DBG");
    let dbg = Dbg::new(sbwt, Some(lcs), n_threads);
    let marks = AtomicBitmap::new(sbwt.n_sets());

    log::info!("Searching first and last k-mer of every input sequence");
    let workers: Vec<_> = (0..n_threads).map(|_| FirstLastKmerWorker {
        in_neighbor_buf: vec![],
        rev_comp_buf: vec![],
        k: sbwt.k(),
        sbwt,
        dbg: &dbg,
        marks: &marks,
        add_rev_comp
    }).collect();

    crate::work_dispatcher::dispatch_work(seq_stream_gen, workers, n_parser_threads, 1 << 23);

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
    // Requires select support on the sbwt
    #[cfg(test)]
    fn new(sbwt: &SbwtIndex<SubsetMatrix>, lcs: Option<&LcsArray>, sample_distance: usize, colex_to_color_set_id: Vec<usize>, n_distinct_color_sets: usize, n_threads: usize) -> Self {

        let get_colorset_fn = |colex| colex_to_color_set_id[colex]; // TODO: this actually returns a color set id. Rename here and later.
        let mut sampling_marks = Self::pick_sampled_kmers(sample_distance, sbwt, lcs, get_colorset_fn, n_threads);

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

    pub fn colex_to_color_set_id(&self, colex: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> usize {
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

pub(crate) struct UnitigImportSeqBatch {
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

    pub(crate) fn process(mut self, results_out: &mut Vec<(usize, usize)>, index: &sbwt::StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>, sample_distance: usize) { // Todo this should consume that batch since it's edited
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
fn unitig_import_parser_thread(unitig_dump: impl std::io::BufRead + Send + Sync + 'static, buf_cap: usize, out: Sender<UnitigImportSeqBatch>){
        
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

    /// Easy but inefficient constructor for tests
    /// colored_seqs has pairs (seq, color)
    #[cfg(test)]
    pub fn new_from_small_input(colored_seqs: &[(&[u8], usize)], k: usize, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {
        use std::collections::HashMap;
        use sbwt::StreamingIndex;
        
        // n_colors = max color id + 1
        let n_colors = colored_seqs.iter().map(|(_seq, color)| *color).max().unwrap() + 1;

        let just_seqs: Vec<&[u8]> = colored_seqs.iter().map(|x| x.0).collect();

        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .algorithm(sbwt::BitPackedKmerSortingMem::new())
            .k(k)
            .build_select_support(true) // Required for colex map creation from sbwt
            .build_lcs(true)
            .run_from_slices(&just_seqs);

        let lcs = lcs.unwrap(); // .build_lcs(true)

        // Build the color sets for each colex position
        log::info!("Building color sets");
        let mut color_sets = Vec::<Vec::<usize>>::new();
        color_sets.resize(sbwt.n_sets(), vec![]);

        let si = StreamingIndex::new(&sbwt, &lcs);
        for (seq, color) in colored_seqs.iter() {
            si.matching_statistics_iter(seq)
            .skip(sbwt.k()-1)
            .filter(|(len,_)| *len == sbwt.k())
            .for_each(|(_, range)|{
                assert_eq!(range.len(), 1);
                color_sets[range.start].push(*color);
            })
        }

        for cset in color_sets.iter_mut() {
            cset.sort();
            cset.dedup();
        }

        let mut distinct_sets_to_id = HashMap::<Vec<usize>, usize>::new(); // Color set -> id
        let mut distinct_sets_in_order = Vec::<Vec::<usize>>::new();
        for cset in color_sets.iter() {
            if !distinct_sets_to_id.contains_key(cset) {
                distinct_sets_to_id.insert(cset.clone(), distinct_sets_to_id.len());
                distinct_sets_in_order.push(cset.clone());
            }
        }

        let n_distinct_sets = distinct_sets_to_id.len();
        log::info!("{} distinct color sets found", n_distinct_sets);

        let mut colex_to_color_set_id = vec![0_usize; sbwt.n_sets()];
        for colex in 0..sbwt.n_sets() {
            colex_to_color_set_id[colex] = distinct_sets_to_id[&color_sets[colex]];
        }

        let map = ColexToColorSetMap::new(&sbwt, Some(&lcs), sample_distance, colex_to_color_set_id, n_distinct_sets, n_threads);

        let color_set_iter = VecVecUsizeIteratorGenerator::new(distinct_sets_in_order);
        let css = CSS::new(color_set_iter, n_colors);
        Self::new(sbwt, lcs, map, *css, None)

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
        unitig_dump: impl std::io::BufRead + Send + Sync + 'static, 
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
        rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap().install(|| {
            colex_to_color_set_id.par_sort_unstable(); // Sorts by colex
        });

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
        2_u64
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        out.write_all(Self::serialization_magic_string()).unwrap();
        out.write_all(&Self::serialization_version().to_le_bytes()).unwrap();

        // Color names are written first (right after magic + version) so that
        // they can be loaded without having to read the rest of the index.
        let n_names = self.color_names.len() as u64;
        out.write_all(&n_names.to_le_bytes()).unwrap();
        for name in self.color_names.iter() {
            let name_bytes = name.as_bytes();
            let name_len = name_bytes.len() as u64;
            out.write_all(&name_len.to_le_bytes()).unwrap();
            out.write_all(name_bytes).unwrap();
        }

        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();
        self.sets.serialize(out);
        self.map.serialize(out);
    }

    /// Reads and validates the magic string and serialization version from `input`.
    fn check_magic_and_version(input: &mut impl std::io::Read) {
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
    }

    fn load_color_names_internal(input: &mut impl std::io::Read) -> Vec<String> {
        let mut n_names_bytes = [0_u8; 8];
        input.read_exact(&mut n_names_bytes).unwrap();
        let n_names = u64::from_le_bytes(n_names_bytes) as usize;
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
        color_names
    }

    /// Loads only the color names from a serialized index, skipping the rest of
    /// the on-disk structure. The remaining bytes of `input` are not consumed.
    pub fn load_color_names_only(input: &mut impl std::io::Read) -> Vec<String> {
        Self::check_magic_and_version(input);
        Self::load_color_names_internal(input)
    }

    /// Load the index from the serialization format.
    /// If this struct is going to be merged with [crate::coloring::merge_colorings], it will
    /// need select support on the sbwt. We need to build it already during loading because
    /// once the sbwt is put on to the heap into an Arc, it cannot be modified anymore.
    /// Unless we make it an Arc<Refcell<...>>, but that might have overhead because then
    /// it will do run-time borrow checking on every access if I understand correctly.
    pub fn load(input: &mut impl std::io::Read, enable_select: bool) -> Self {
        Self::check_magic_and_version(input);
        let color_names = Self::load_color_names_internal(input);

        log::info!("Loading the SBWT");
        let mut sbwt = SbwtIndex::<SubsetMatrix>::load(input).unwrap();
        if enable_select {
            log::info!("Building select support");
            sbwt.build_select();
        }
        log::info!("Loading the LCS array");
        let lcs = LcsArray::load(input).unwrap();
        log::info!("Loading color sets");
        let sets = CSS::load(input);
        log::info!("Loading color set mapping");
        let map = ColexToColorSetMap::load(input);
        assert_eq!(map.sampling.len(), sbwt.n_sets());
        assert_eq!(color_names.len(), sets.n_colors());

        log::info!("Index loaded");
        CompactColexKmers{sbwt, lcs, sets, map, color_names}
    }

    pub fn lookup_kmer_color_sets(&self, seq: &[u8]) -> Vec<Option<CSS::SetView<'_>>> {
        let mut set_ids = Vec::<Option<usize>>::new();
        self.push_color_set_ids_to_buffer(seq, &mut set_ids);

        let mut sets = Vec::<Option::<CSS::SetView<'_>>>::with_capacity(set_ids.len());
        for (idx, id) in set_ids.iter().enumerate() {
            if idx == 0 || (set_ids[idx] != set_ids[idx-1]) {
                // If the kmer existed, its color set different than the previous
                match id {
                    Some(id) => {
                        // k-mer exists -> fetch color set
                        sets.push(Some(self.set_id_to_set(*id)))
                    },
                    None => {
                        // Absent k-mer
                        sets.push(None)
                    }
                }
            } else {
                // Same set as previous (possibly None)
                sets.push(sets.last().unwrap().clone())
            }
        }
        sets
    }

    // Does not clear the buffer
    pub fn push_color_set_ids_to_buffer(&self, seq: &[u8], buffer: &mut Vec<Option<usize>>) {
        let k = self.sbwt.k();
        if seq.len() < k {
            return;
        }

        let si = sbwt::StreamingIndex::new(&self.sbwt, &self.lcs);

        #[derive(Eq, PartialEq, Debug)]
        enum ColorId { // Contains colex rank in the value, if exists
            Sampled(usize),
            SameAsNext(usize),
            AbsentKmer
        }

        let mut colex_positions = Vec::<ColorId>::with_capacity(seq.len()-k+1);

        // Pass 1: Compute colex ranks and whether the k-mers are sampled
        for (len, range) in si.matching_statistics_iter(seq).skip(k-1) {
            if len == k {
                assert!(range.len() == 1);
                let colex = range.start;
                if self.map.sampling.get(colex) {
                    colex_positions.push(ColorId::Sampled(colex));
                } else {
                    colex_positions.push(ColorId::SameAsNext(colex));
                }
            } else {
                colex_positions.push(ColorId::AbsentKmer);
            }
        }

        // Pass 2: look up set ids right-to-left, reusing set ids for the next k-mer whenever possible
        let old_buf_len = buffer.len();
        buffer.resize(old_buf_len + colex_positions.len(), None); // Make space for results
        let set_ids_output = &mut buffer[old_buf_len..];
        for i in (0..colex_positions.len()).rev() {
            match colex_positions[i] {
                ColorId::Sampled(colex) => {
                    let set_id = self.map.colex_to_color_set_id(colex, &self.sbwt);
                    set_ids_output[i] = Some(set_id);
                },
                ColorId::SameAsNext(colex) => {
                    if i+1 == set_ids_output.len() || set_ids_output[i+1].is_none() {
                        // Can not copy color set id from position i+1
                        let set_id = self.map.colex_to_color_set_id(colex, &self.sbwt);
                        set_ids_output[i] = Some(set_id);
                    } else {
                        set_ids_output[i] = set_ids_output[i+1];
                    }
                },
                ColorId::AbsentKmer => {
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
