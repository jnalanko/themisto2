use std::{io::Write, ops::Range, path::Path, sync::atomic::{AtomicUsize, Ordering::Relaxed}};
use crossbeam::channel::RecvTimeoutError;
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

impl PseudoalignmentBatch {
    fn new(first_seq_rank: usize) -> Self {
        Self { seqs: SeqDB::new(), seq_ranks: first_seq_rank..first_seq_rank }
    }

    fn push(&mut self, seq: &[u8]) {
        self.seqs.push_seq(seq);
        self.seq_ranks.end += 1;
    }

    fn process<CSS: ColorSetStorage>(self, aligner: &mut Box<dyn Pseudoaligner<CSS> + Send>, index: &CompactColexKmers<CSS>, n_bases_processed: &AtomicUsize) -> PseudoalignmentBatchResult {
        let mut result = PseudoalignmentBatchResult {
            concat: vec![],
            starts: vec![0],
            seq_ranks: self.seq_ranks,
        };

        for rec in self.seqs.iter() {
            aligner.push_compatible_colors(rec.seq, index, &mut result.concat);
            result.starts.push(result.concat.len());
            n_bases_processed.fetch_add(rec.seq.len(), Relaxed);
        }
        result
    }
}

impl PseudoalignmentBatchResult {
    fn get_result_set(&self, idx: usize) -> &[usize] {
        &self.concat[self.starts[idx]..self.starts[idx+1]]
    } 
}

// This uses a factory pattern for creating new pseudoaligners. I'm so sorry.
// But it actually makes sense here: I want that the pseudoalignment function
// can create a separate aligner for each worker thread, but so that
// it does not have to care how they are constructed.
// The output callback takes pairs (read rank in input, pseudoaligned color ids)
fn run_all_queries<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, mut queries: impl SeqStream + Send + 'static, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send>, mut output_callback: impl FnMut((usize, &[usize])) + Send, n_workers: usize) {

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
                    let result = batch.process(&mut aligner, index_ref, n_bases_processed_ref);
                    results_send_clone.send(result).unwrap();
                }
            });
            worker_handles.push(handle);
        }

        let outputter_handle = scope.spawn(|| {
            while let Ok(result) = results_recv.recv() {
                for (idx, query_rank) in result.seq_ranks.clone().enumerate() {
                    output_callback((query_rank, result.get_result_set(idx)));
                }
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
                        let n = n_bases_processed.load(Relaxed) - last_n_bases_processed;
                        let t = last_wakeup_time.elapsed().as_secs();
                        let throughput = n as f64 / t as f64 / (1 << 20) as f64;
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
            let total_t = start_time.elapsed().as_secs();
            let total_throughput = total_n as f64 / total_t as f64 / (1 << 20) as f64;
            log::info!("Total bases processed: {}", total_n);
            log::info!("Total throughput: {:.3}", total_throughput);
            //std::thread::sleep(std::time::Duration::from_secs(sleep_interval_seconds));
        });
        
        parser_handle.join().unwrap(); // Wait for the parser to finish
        for h in worker_handles { h.join().unwrap() } // Wait for the workers to finish
        drop(results_send); // Signal that no more results will be pushed
        outputter_handle.join().unwrap(); // Wait for the outputter to finish
        progress_printer_quit_signal_send.send(()).unwrap(); // Interrupt the progress printer from sleep
        progress_printer_handle.join().unwrap();
    }); 
}

pub fn run_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, input_file: &Path, mut output: impl Write + Send, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send> + 'static, n_aligners: usize) {
    let reader = crate::io::ChainedInputStream::new(vec![input_file.to_path_buf()]);

    let output_callback = |(read_rank, color_ids): (usize, &[usize])| {
        write!(output, "{}", read_rank).unwrap();
        for cid in color_ids {
            write!(output, " {}", cid).unwrap();
        }
        writeln!(output).unwrap();
    };
    run_all_queries(index, reader, create_new_aligner, output_callback, n_aligners);
}