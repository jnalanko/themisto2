use std::{cmp::Reverse, io::Write, ops::Range, path::Path, sync::atomic::{AtomicUsize, Ordering::Relaxed}};
use crossbeam::channel::{RecvTimeoutError, Sender};
use jseqio::seq_db::SeqDB;
use sbwt::SeqStream;
use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView}};

pub trait Pseudoaligner<CSS: ColorSetStorage> {

    // The &mut self is to access and modify thread-local buffers
    // owned by the algorithm.
    fn push_compatible_colors(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, output: &mut Vec<usize>);
}

#[derive(Clone, Copy, Debug)]
pub enum Denominator {
    All,
    Relevant,
    MaxHits,
}

#[derive(Clone)]
pub struct ThresholdPseudoaligner {
    counts: Vec<usize>,
    nonzero_count_indices: Vec<usize>,
    threshold: f64,
    denominator: Denominator,
    min_hits: usize,
}

impl ThresholdPseudoaligner {
    pub fn new(n_colors: usize, threshold: f64, min_hits: usize, denominator: Denominator) -> Self {
        Self {
            counts: vec![0; n_colors],
            nonzero_count_indices: vec![],
            threshold,
            min_hits,
            denominator,
        }
    }
}

impl<CSS: ColorSetStorage> Pseudoaligner<CSS> for ThresholdPseudoaligner {
    fn push_compatible_colors(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>) {
        let mut n_relevant = 0_usize;
        let mut n_all = 0_usize;
        for set in index.lookup_kmer_color_sets(seq) {
            if let Some(set) = set {
                for color in set.iter() {
                    self.counts[color] += 1;
                    if self.counts[color] == 1 {
                        self.nonzero_count_indices.push(color);
                    }
                }
                n_relevant += 1;
            }
            n_all += 1;
        }
        self.nonzero_count_indices.sort_unstable(); // Sort to output in sorted order by colors 

        // Add to output all colors that pass the threshold
        if n_relevant >= self.min_hits && self.nonzero_count_indices.len() > 0 {
            let den = match self.denominator {
                Denominator::All => n_all as f64,
                Denominator::Relevant => n_relevant as f64,
                Denominator::MaxHits => {
                    let maxhits = self.nonzero_count_indices.iter().map(|color| self.counts[color]).max();
                    maxhits.unwrap() as f64 // Safe because here nonzero_count_indices.len() > 0
                },
            };
            for color in self.nonzero_count_indices.iter() {
                if self.counts[color] as f64 / den >= self.threshold {
                    out.push(color);
                }
            }
        }

        // Clean up
        for &color in &self.nonzero_count_indices {
            self.counts[color] = 0;
        }
        self.nonzero_count_indices.clear();
    }
}

#[derive(Clone)]
pub struct IntersectionPseudoaligner {
    min_hits: usize,
}

impl IntersectionPseudoaligner{
    pub fn new(min_hits: usize) -> Self {
        Self {
            min_hits,
        }
    }
}

impl<CSS: ColorSetStorage> Pseudoaligner<CSS> for IntersectionPseudoaligner {
    fn push_compatible_colors(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, output: &mut Vec<usize>) {
        let mut intersection = index.get_set_storage().get_full_set();
        let mut n_hits = 0_usize;
        #[allow(clippy::manual_flatten)] // Clearer this way
        for set in index.lookup_kmer_color_sets(seq) {
            if let Some(set) = set {
                index.get_set_storage().intersect(&mut intersection, &set);
                n_hits += 1;
            }
        }

        if n_hits >= self.min_hits {
            for color in intersection.iter() {
                output.push(color);
            }
        }
    }
}


struct PseudoalignmentBatch {
    seqs: SeqDB,
    seq_ranks: Range<usize>
}

struct PseudoalignmentBatchResult {
    concat: Vec<usize>,
    starts: Vec<usize>, // Has an end sentinel one past the end
    seq_ranks: Range<usize>, // Original sequence ranks for these queries
}

impl PartialEq for PseudoalignmentBatchResult {
    fn eq(&self, other: &Self) -> bool {
        self.seq_ranks == other.seq_ranks
    }
}
impl Eq for PseudoalignmentBatchResult {}

impl PartialOrd for PseudoalignmentBatchResult {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other)) // Using the total order
    }
}

impl Ord for PseudoalignmentBatchResult {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.seq_ranks.start.cmp(&other.seq_ranks.start)
    }
}

impl PseudoalignmentBatch {
    const SEND_TO_OUTPUT_THRESHOLD: usize = 1_000_000;

    fn new(first_seq_rank: usize) -> Self {
        Self { seqs: SeqDB::new(), seq_ranks: first_seq_rank..first_seq_rank }
    }

    fn push(&mut self, seq: &[u8]) {
        self.seqs.push_seq(seq);
        self.seq_ranks.end += 1;
    }

    fn process<CSS: ColorSetStorage>(self, aligner: &mut Box<dyn Pseudoaligner<CSS> + Send>, index: &CompactColexKmers<CSS>, n_bases_processed: &AtomicUsize, output_channel: &Sender<PseudoalignmentBatchResult>) {
        let mut result = PseudoalignmentBatchResult::new(self.seq_ranks.start);
        for rec in self.seqs.iter() {
            aligner.push_compatible_colors(rec.seq, index, &mut result.concat);
            result.starts.push(result.concat.len());
            result.seq_ranks.end += 1;

            if result.concat.len() > Self::SEND_TO_OUTPUT_THRESHOLD {
                // Send current results for output
                let next_start = result.seq_ranks.end;
                output_channel.send(result).unwrap();
                result = PseudoalignmentBatchResult::new(next_start);
            }

            n_bases_processed.fetch_add(rec.seq.len(), Relaxed);
        }

        // Output remaining results
        if result.seq_ranks.len() > 0 {
            output_channel.send(result).unwrap();
        }
    }
}

impl PseudoalignmentBatchResult {

    fn new(first_seq_rank: usize) -> Self {
        Self {
            concat: vec![],
            starts: vec![0],
            seq_ranks: first_seq_rank..first_seq_rank,
        }
    }

    fn get_result_set(&self, idx: usize) -> &[usize] {
        &self.concat[self.starts[idx]..self.starts[idx+1]]
    } 

    fn write(&self, out: &mut impl std::io::Write) {
        for seq_rank in self.seq_ranks.clone() {
            write!(out, "{}", seq_rank).unwrap();
            for cid in self.get_result_set(seq_rank - self.seq_ranks.start){
                write!(out, " {}", cid).unwrap();
            }
            writeln!(out).unwrap();
        }
    }
}

// This uses a factory pattern for creating new pseudoaligners. I'm so sorry.
// But it actually makes sense here: I want that the pseudoalignment function
// can create a separate aligner for each worker thread, but so that
// it does not have to care how they are constructed.
// The output callback takes pairs (read rank in input, pseudoaligned color ids)
fn run_all_queries<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, mut queries: impl SeqStream + Send + 'static, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send>, mut output_callback: impl FnMut(PseudoalignmentBatchResult) + Send, n_workers: usize) {

    let batch_size = 10_000_usize;
    let (work_send, work_recv) = crossbeam::channel::bounded::<PseudoalignmentBatch>(n_workers);
    let (results_send, results_recv) = crossbeam::channel::bounded::<PseudoalignmentBatchResult>(n_workers);

    let (progress_printer_quit_signal_send, progress_printer_quit_signal_recv) = crossbeam::channel::bounded::<()>(1);

    let n_bases_processed = AtomicUsize::new(0); 

    std::thread::scope(|scope| {
        let parser_handle = scope.spawn(move || {
            let mut n_seqs_read = 0_usize;
            let mut cur_batch = PseudoalignmentBatch::new(0);
            while let Some(q) = queries.stream_next() {
                n_seqs_read += 1;
                cur_batch.push(q);
                if cur_batch.seqs.total_seq_len() >= batch_size {
                    work_send.send(cur_batch).unwrap();
                    cur_batch = PseudoalignmentBatch::new(n_seqs_read);
                }
            }
            if cur_batch.seqs.total_seq_len() > 0 { // Last batch
                work_send.send(cur_batch).unwrap();
            }
            drop(work_send); // Signal that no more work is going to be pushed
        });

        let mut worker_handles = vec![];
        for _worker_id in 0..n_workers {
            let mut aligner = create_new_aligner();
            let work_recv_clone = work_recv.clone();
            let results_send_clone = results_send.clone();
            let index_ref = index;
            let n_bases_processed_ref = &n_bases_processed;
            let handle = scope.spawn(move || {
                while let Ok(batch) = work_recv_clone.recv() {
                    batch.process(&mut aligner, index_ref, n_bases_processed_ref, &results_send_clone);
                }
            });
            worker_handles.push(handle);
        }

        let outputter_handle = scope.spawn(|| {
            while let Ok(result) = results_recv.recv() {
                output_callback(result)
            }
        });

        let progress_printer_handle = scope.spawn(|| {
            let mut last_wakeup_time = std::time::Instant::now();
            let mut last_n_bases_processed = n_bases_processed.load(Relaxed);
            let print_interval = std::time::Duration::from_secs(10);
            let start_time = std::time::Instant::now();
            loop {
                match progress_printer_quit_signal_recv.recv_timeout(print_interval) {
                    Ok(_) => break, // Received the quit signal
                    Err(RecvTimeoutError::Timeout) => { // Time to print
                        let n = n_bases_processed.load(Relaxed);
                        let t = last_wakeup_time.elapsed().as_secs_f64();
                        let throughput = (n - last_n_bases_processed) as f64 / t / (1 << 20) as f64;
                        log::info!("Current throughput {:.3} Mbases/s ({} bases processed total)", throughput, n);
                        last_n_bases_processed = n;
                        last_wakeup_time = std::time::Instant::now();
                    },
                    Err(RecvTimeoutError::Disconnected) => {
                        // I'm not sure when this would happen, but let's just quit
                       break
                    }
                }
            }

            // Print total statistics
            let total_n = n_bases_processed.load(Relaxed);
            let total_t = start_time.elapsed().as_secs_f64();
            let total_throughput = total_n as f64 / total_t / (1 << 20) as f64;
            log::info!("Total bases {} bases processed in {:.3} seconds", total_n, total_t);
            log::info!("Total throughput: {:.3} Mbases/s", total_throughput);
        });
        
        parser_handle.join().unwrap(); // Wait for the parser to finish
        for h in worker_handles { h.join().unwrap() } // Wait for the workers to finish
        drop(results_send); // Signal that no more results will be pushed
        outputter_handle.join().unwrap(); // Wait for the outputter to finish
        progress_printer_quit_signal_send.send(()).unwrap(); // Interrupt the progress printer from sleep
        progress_printer_handle.join().unwrap();
    }); 
}

struct OutputWriter {
    buffer: std::collections::BinaryHeap::<Reverse<PseudoalignmentBatchResult>>, // Reverse makes this a min-heap
    next_seq_rank: usize,
}

impl OutputWriter {
    fn new() -> Self {
        Self {
            buffer: std::collections::BinaryHeap::<Reverse<PseudoalignmentBatchResult>>::new(),
            next_seq_rank: 0,
        }
    }

    // Does not immediately write if we're missing some batch result that should
    // be written ealier.
    fn push_batch(&mut self, result_batch: PseudoalignmentBatchResult, mut output: &mut impl Write) {
        self.buffer.push(Reverse(result_batch)); // Reverse makes this a min heap
        loop { // Print all batches that can now be printed
            let min_batch = self.buffer.peek();
            if let Some(min_batch) = min_batch {
                let min_batch = &min_batch.0; // Unwrap from Reverse
                if min_batch.seq_ranks.start == self.next_seq_rank {
                    min_batch.write(&mut output);
                    self.next_seq_rank = min_batch.seq_ranks.end;
                    self.buffer.pop();
                } else {
                    break; // Not ready to print min_batch yet
                }
            } else {
                break; // Batch buffer is empty
            }
        }
    }
}

pub fn run_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, input_file: &Path, mut output: impl Write + Send, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send> + 'static, n_aligners: usize) {
    let reader = crate::io::ChainedInputStream::new(vec![input_file.to_path_buf()]);

    // Create the output callback
    let mut writer = OutputWriter::new();
    let output_callback = |result_batch| {
        writer.push_batch(result_batch, &mut output);
    };

    run_all_queries(index, reader, create_new_aligner, output_callback, n_aligners);
    assert!(writer.buffer.is_empty());
}