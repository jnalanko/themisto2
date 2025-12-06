use std::{io::Write, marker::PhantomData, path::Path, sync::atomic::{AtomicUsize, Ordering::Relaxed}};

use crossbeam::channel::{Receiver, RecvTimeoutError, Sender};
use jseqio::{record::Record, seq_db::SeqDB};
use rand_distr::num_traits::ConstOne;

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView}};

trait Pseudoaligner<CSS: ColorSetStorage> {
    // The &mut self is to allow internal state containing reused buffers
    fn push_compatibility_set(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>);
}

struct IntersectionPseudoalignment<CSS: ColorSetStorage> {
    min_hits: usize,
    css: PhantomData<CSS>, // TODO: does this make sense?
}

#[derive(Clone, Copy, Debug)]
pub enum Denominator {
    All,
    Relevant,
    MaxHits,
}

#[derive(Clone)]
pub struct ThresholdPseudoaligner<CSS: ColorSetStorage> {
    counts: Vec<usize>,
    nonzero_count_indices: Vec<usize>,
    threshold: f64,
    denominator: Denominator,
    min_hits: usize,
    css: PhantomData<CSS>, // TODO: does this make sense?
}

impl<CSS: ColorSetStorage> ThresholdPseudoaligner<CSS> {
    pub fn new(n_colors: usize, threshold: f64, min_hits: usize, denominator: Denominator) -> Self {
        Self {
            counts: vec![0; n_colors],
            nonzero_count_indices: vec![],
            threshold,
            min_hits,
            denominator,
            css: PhantomData,
        }
    }
}

impl<CSS: ColorSetStorage> Pseudoaligner<CSS> for ThresholdPseudoaligner<CSS> {
    fn push_compatibility_set(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>) {
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

impl<CSS: ColorSetStorage> Pseudoaligner<CSS> for IntersectionPseudoalignment<CSS> {
    fn push_compatibility_set(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>) {
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
                out.push(color);
            }
        }
    }
}

#[derive(Copy, Clone)]
enum Metric {
    KmerHits,
    BasesCovered,
    AlignmentLength,
    LongestMatchRun,
    ShortestGap
}

// todo: can be much faster 
#[allow(clippy::manual_flatten)]
fn compute_kmer_hits_to_compatible_colors<CSS: ColorSetStorage>(color_set_ids: &[Option<usize>], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<usize> {
    let mut hits = vec![0; index.get_set_storage().n_colors()];
    for color_set_id_opt in color_set_ids {
        if let Some(color_set_id) = color_set_id_opt {
            let color_set = index.set_id_to_set(*color_set_id);
            for color in color_set.iter() {
                hits[color] += 1;
                // TODO: Faster: If same color id appears multiple times, increment by the multiplicity
            }
        }
    }

    // Return only hits to compatible colors
    compatible_colors.iter().map(|&c| hits[c]).collect()
}

struct QueryBatch {
    seqs: SeqDB,
    metrics: Vec<Metric>, // These metric should be computed
}

impl QueryBatch {
    fn new() -> Self { // TODO: take metrics
        Self {
            seqs: SeqDB::new(),
            metrics: vec![],
        }
    }

    // Returns JSON-formatted bytes
    fn process<CSS: ColorSetStorage>(self, index: &CompactColexKmers<CSS>, cc: &mut Box<dyn Pseudoaligner<CSS> + Send>, n_bases_processed: &AtomicUsize) -> Vec<u8> {
        let mut result = QueryResult::new();
        let mut compat_set_buf = Vec::<usize>::new();
        for rec in self.seqs.iter() {
            compat_set_buf.clear();

            cc.push_compatibility_set(rec.seq, index, &mut compat_set_buf);
            result.push(&compat_set_buf, rec.name());

            let metric_vecs = self.compute_metrics(&rec.seq, &compat_set_buf, index);
            result.metrics.push(metric_vecs);

            n_bases_processed.fetch_add(rec.seq.len(), Relaxed);
        }
        let mut bytes_out = Vec::<u8>::new();
        result.into_json(&mut bytes_out);

        bytes_out
    }

    fn compute_metrics<CSS: ColorSetStorage>(&self, seq: &[u8], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<(Metric, Vec<usize>)>{
        let mut ans = Vec::<(Metric, Vec<usize>)>::new();

        if self.metrics.len() > 0 {
            let mut color_set_ids = Vec::<Option<usize>>::new();
            index.push_color_set_ids_to_buffer(seq, &mut color_set_ids);
            for metric in self.metrics.iter() {
                let values = match metric {
                    Metric::KmerHits => compute_kmer_hits_to_compatible_colors(&color_set_ids, compatible_colors, index),
                    Metric::BasesCovered => todo!(),
                    Metric::AlignmentLength => todo!(),
                    Metric::LongestMatchRun => todo!(),
                    Metric::ShortestGap => todo!(),
                };
                ans.push((*metric, values));
            }
        }
        ans
    }
}

struct QueryResult {
    query_names_concat: Vec<u8>,
    query_names_starts: Vec<usize>,

    compatibility_class_concat: Vec<usize>,
    compatibility_class_starts: Vec<usize>,

    // Optional metrics: For each query sequence, a vector of pairs
    // (Metric, Vec<usize>), where the length of the vector in the pair
    // is equal to the compatibility class size, and has metric values
    // for each color in the compatibility class in the same order as
    // the colors appear in the compatibility class.
    metrics: Vec<Vec<(Metric, Vec<usize>)>> // TODO: concatenation
}

impl QueryResult {

    fn new() -> Self {
        Self { 
            query_names_concat: vec![], 
            query_names_starts: vec![0], 
            compatibility_class_concat: vec![], 
            compatibility_class_starts: vec![0], 
            metrics: vec![],
        }
    }

    fn into_json(self, out: &mut impl Write) {
        todo!();
    }

    fn push(&mut self, compat_set: &[usize], seq_name: &[u8]) {
        self.compatibility_class_concat.extend_from_slice(compat_set);
        self.compatibility_class_starts.push(self.compatibility_class_concat.len());

        self.query_names_concat.extend_from_slice(seq_name);
        self.query_names_starts.push(self.query_names_concat.len());
    }
}

struct Worker<'a, CSS: ColorSetStorage> {
    index: &'a  CompactColexKmers<CSS>,
}

pub fn run_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, input_file: &Path, mut output: impl Write + Send, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send> + 'static, n_aligners: usize) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&input_file).unwrap();

    let batch_size = 10_000_usize;
    let (work_send, work_recv) = crossbeam::channel::bounded::<QueryBatch>(n_aligners);
    let (results_send, results_recv) = crossbeam::channel::bounded::<Vec<u8>>(n_aligners); // Json-formatted blocks of text

    let (progress_printer_quit_signal_send, progress_printer_quit_signal_recv) = crossbeam::channel::bounded::<()>(1);

    let n_bases_processed = AtomicUsize::new(0); 

    std::thread::scope(|scope| {
        let parser_handle = scope.spawn(move || {
            let mut cur_batch = QueryBatch::new();
            while let Some(q) = reader.read_next().unwrap() {
                cur_batch.seqs.push_record(q);
                if cur_batch.seqs.total_seq_len() >= batch_size {
                    work_send.send(cur_batch).unwrap();
                    cur_batch = QueryBatch::new();
                }
            }
            if cur_batch.seqs.total_seq_len() > 0 { // Last batch
                work_send.send(cur_batch).unwrap();
            }
            drop(work_send); // Signal that no more work is going to be pushed
        });

        let mut worker_handles = vec![];
        for _worker_id in 0..n_aligners {
            let mut aligner = create_new_aligner();
            let work_recv_clone = work_recv.clone();
            let results_send_clone = results_send.clone();
            let index_ref = index;
            let n_bases_processed_ref = &n_bases_processed;
            let handle = scope.spawn(move || {
                while let Ok(batch) = work_recv_clone.recv() {
                    let json = batch.process(index_ref, &mut aligner, n_bases_processed_ref);
                    results_send_clone.send(json);
                }
            });
            worker_handles.push(handle);
        }

        let outputter_handle = scope.spawn(|| {
            while let Ok(json) = results_recv.recv() {
                output.write_all(&json).unwrap();
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