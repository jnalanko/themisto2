use std::{io::Write, path::Path, sync::atomic::{AtomicUsize, Ordering::Relaxed}};

use crossbeam::channel::{RecvTimeoutError};
use jseqio::{record::Record, seq_db::SeqDB};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView}, pseudoalignment_metrics::{Metric, PseudoalignmentMetricProcessor, create_metric_processor}};

struct SortedVec {
    v: Vec<usize>
}

impl SortedVec {
    fn new(mut v: Vec<usize>) -> Self {
        v.sort();
        Self{v}
    }
}

pub trait Pseudoaligner<CSS: ColorSetStorage> {
    // The &mut self is to allow internal state containing reused buffers
    fn push_compatibility_set(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>);
}

pub struct IntersectionPseudoaligner {
    min_hits: usize,
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
    fn push_compatibility_set(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>) {
        let mut n_relevant = 0_usize;
        let mut n_all = 0_usize;

        let mut color_set_ids = Vec::<Option::<usize>>::new();
        index.push_color_set_ids_to_buffer(seq, &mut color_set_ids);
        crate::util::for_each_run(&color_set_ids, |run_range| {
            let run_len = run_range.len(); 
            assert!(run_len > 0);

            let first_id = color_set_ids[run_range.start];
            if let Some(set_id) = first_id {
                for color in index.set_id_to_set(set_id).iter() {
                    if self.counts[color] == 0 {
                        self.nonzero_count_indices.push(color);
                    }
                    self.counts[color] += run_len;
                }
                n_relevant += run_len;
            }

            n_all += run_len; // Runs of None count here

        });

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

impl IntersectionPseudoaligner {
    pub fn new(min_hits: usize) -> Self {
        Self {
            min_hits,
        }
    }
}

impl<CSS: ColorSetStorage> Pseudoaligner<CSS> for IntersectionPseudoaligner {
    fn push_compatibility_set(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>) {
        let mut intersection = index.get_set_storage().get_full_set();
        let mut n_hits = 0_usize;

        #[allow(clippy::manual_flatten)] // Clearer this way
        let mut color_set_ids = Vec::<Option::<usize>>::new();
        index.push_color_set_ids_to_buffer(seq, &mut color_set_ids);
        // TODO: do not look everything at once to be able to exit early
        // if the intersection becomes empty.
        // TODO: intersect smallest sets first to speed up the computation.
        crate::util::for_each_run(&color_set_ids, |run_range| {
            // If the current intersection in empty, it will stay empty, so we
            // need to do any work only if the intersection is nonempty.
            if intersection.len() > 0 {
                let run_len = run_range.len(); 
                assert!(run_len > 0);

                let first_id = color_set_ids[run_range.start];
                if let Some(first_id) = first_id {
                    let set = index.set_id_to_set(first_id);
                    index.get_set_storage().intersect(&mut intersection, &set);
                    n_hits += 1;
                }
            }
        });

        if n_hits >= self.min_hits {
            for color in intersection.iter() {
                out.push(color);
            }
        }
    }
}

struct QueryBatch {
    seqs: SeqDB,
    batch_id: usize
}

impl QueryBatch {
    fn new(batch_id: usize) -> Self { // TODO: take metrics
        Self {
            seqs: SeqDB::new(),
            batch_id,
        }
    }

    // Returns JSON-formatted bytes, and the batch id.
    fn process<CSS: ColorSetStorage>(self, index: &CompactColexKmers<CSS>, aligner: &mut Box<dyn Pseudoaligner<CSS> + Send>, n_bases_processed: &AtomicUsize, metrics: &mut [Box<dyn PseudoalignmentMetricProcessor<CSS>>]) -> (Vec<u8>, usize) {
        let mut result = QueryResult::new();
        let mut compat_set_buf = Vec::<usize>::new();
        for rec in self.seqs.iter() {
            compat_set_buf.clear();

            aligner.push_compatibility_set(rec.seq, index, &mut compat_set_buf);

            let sorted_compat_set = SortedVec::new(compat_set_buf);
            let metric_values_concat = self.compute_metrics(rec.seq, &sorted_compat_set, index, metrics);
            compat_set_buf = sorted_compat_set.v; // Move ownership back

            result.push(&compat_set_buf, rec.name(), metric_values_concat);

            n_bases_processed.fetch_add(rec.seq.len(), Relaxed);
        }

        result.computed_metric_names.extend(metrics.iter().map(|proc| proc.metric_id()));

        let mut bytes_out = Vec::<u8>::new();
        result.into_json(&mut bytes_out);

        (bytes_out, self.batch_id)
    }

    fn compute_metrics<CSS: ColorSetStorage>(&self, seq: &[u8], compatible_colors: &SortedVec, index: &CompactColexKmers<CSS>, metrics: &mut [Box<dyn PseudoalignmentMetricProcessor<CSS>>]) -> Vec<usize> {
        let mut ans_concats = Vec::<usize>::new();

        if metrics.len() > 0 {
            let mut color_set_ids = Vec::<Option<usize>>::new();
            index.push_color_set_ids_to_buffer(seq, &mut color_set_ids);
            for processor in metrics.iter_mut() {
                let mut compatible_colors_idx = 0_usize;
                let mut color_value_pairs = processor.process(&color_set_ids, index);
                color_value_pairs.sort();
                for (color, value) in color_value_pairs {
                    // Push this only if this color is compatible. Here we make use of the property that
                    // compatible_colors is sorted.
                    if compatible_colors_idx < compatible_colors.v.len() && compatible_colors.v[compatible_colors_idx] == color {
                        ans_concats.push(value);
                        compatible_colors_idx += 1;
                    }
                }

                // All of the metrics are assumed to be so that if a color is compatible
                // according to pseudoalignment, then we report some value for it.
                // TODO: this is a brittle assumption and might break in the future. Do something about it.
                // For now we just assert that it's true.
                assert_eq!(compatible_colors_idx, compatible_colors.v.len());
            }
        }
        ans_concats
    }
}

struct QueryResult {
    query_names_concat: Vec<u8>,
    query_names_starts: Vec<usize>,

    compatibility_class_concat: Vec<usize>,
    compatibility_class_starts: Vec<usize>,

    // List of metric names that were computed for this batch
    computed_metric_names: Vec<Metric>,

    // For each query, the concatenation like:
    //   [metric 1, color 1] [metric 1, color 2] ... [metric 1, color m]
    //   [metric 2, color 1] [metric 2, color 2] ... [metric 2, color m]
    //   ... 
    //   [metric r, color 1] [metric r, color 2] ... [metric r, color m]
    // The metric values are reported only for the m compatible colors, in the
    // same order as the colors appear in the compatibility class.
    metrics_concat: Vec<Vec<usize>>, 
}

impl QueryResult {

    fn new() -> Self {
        Self { 
            query_names_concat: vec![], 
            query_names_starts: vec![0], 
            compatibility_class_concat: vec![], 
            compatibility_class_starts: vec![0], 
            computed_metric_names: vec![],
            metrics_concat: vec![],
        }
    }

    fn write_slice_as_ascii(v: &[usize], out: &mut Vec<u8>) {
        out.push(b'[');
        for (i, val) in v.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }

            // Push digits to out without intermediate allocation
            write!(out, "{}", val).unwrap(); // Does this allocate?
        }
        out.push(b']');
    }

    fn into_json(self, out: &mut impl Write) {
        // I'm rolling my own JSON serialization because I don't trust
        // that the serde json implementation is fast enough because it
        // must be quite generic.
        let mut bytes = Vec::<u8>::new();
        assert_eq!(self.compatibility_class_starts.len(), self.query_names_starts.len());

        assert!(self.query_names_starts.len() > 0); // There should be at least the 0 at the beginning even if it's empty
        let n_queries = self.query_names_starts.len()-1;

        for seq_idx in 0..n_queries {
            bytes.push(b'{');
            let name_s = self.query_names_starts[seq_idx];
            let name_e = self.query_names_starts[seq_idx+1];
            let compat_s = self.compatibility_class_starts[seq_idx]; 
            let compat_e = self.compatibility_class_starts[seq_idx+1]; 

            bytes.extend(b"\"name\": \"");
            // Escape double quotes. '"' is ASCII (0x22), which never appears as a byte
            // inside a multi-byte UTF-8 sequence, so escaping byte-by-byte is Unicode-safe.
            for &b in &self.query_names_concat[name_s..name_e] {
                if b == b'"' {
                    bytes.push(b'\\');
                }
                bytes.push(b);
            }

            bytes.extend(b"\", \"colors\": ");
            Self::write_slice_as_ascii(&self.compatibility_class_concat[compat_s..compat_e], &mut bytes);

            for (metric_idx, metric) in self.computed_metric_names.iter().enumerate() {
                match metric {
                    Metric::KmerHits => {
                        bytes.extend(b", \"kmer_hits\": ");
                    },
                    Metric::BasesCovered => {
                        bytes.extend(b", \"bases_covered\": ");
                    },
                }
                let slice_start = metric_idx * (compat_e-compat_s);
                let slice_end = (metric_idx+1) * (compat_e-compat_s);
                let values = &self.metrics_concat[seq_idx][slice_start..slice_end];
                Self::write_slice_as_ascii(values, &mut bytes);
            }
            bytes.push(b'}');
            bytes.push(b'\n');
        }
        out.write_all(&bytes).unwrap();
    }

    fn push(&mut self, compat_set: &[usize], seq_name: &[u8], metric_values_concat: Vec<usize>) {
        self.compatibility_class_concat.extend_from_slice(compat_set);
        self.compatibility_class_starts.push(self.compatibility_class_concat.len());

        self.query_names_concat.extend_from_slice(seq_name);
        self.query_names_starts.push(self.query_names_concat.len());

        self.metrics_concat.push(metric_values_concat);
    }
}

fn output_thread(results_recv: crossbeam::channel::Receiver<(Vec<u8>, usize)>, mut output: impl Write + Send, sort_output: bool, n_bytes_output: &AtomicUsize){
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;

    if !sort_output {
        // Write batches in the order they are received.
        while let Ok((json, _batch_id)) = results_recv.recv() {
            output.write_all(&json).unwrap();
            n_bytes_output.fetch_add(json.len(), Relaxed);
        }
        return;
    }

    let mut next_batch_id = 0_usize;
    // Min-heap of batches received out of order, keyed by batch id.
    let mut buffer: BinaryHeap<Reverse<(usize, Vec<u8>)>> = BinaryHeap::new();

    let mut next_buffer_size_warning = 1_usize << 31; // 2 GiB 
    let mut total_buffer_size = 0_usize;

    while let Ok((json, batch_id)) = results_recv.recv() {

        total_buffer_size += json.len();
        buffer.push(Reverse((batch_id, json)));

        if total_buffer_size >= next_buffer_size_warning {
            let human_bytes = human_bytes::human_bytes(total_buffer_size as f64);
            log::warn!("A large number of bytes waiting to be written to disk: {}. Consider using fewer threads, or running without output sorting.", human_bytes);
            next_buffer_size_warning *= 2;
        }

        // Write out all batches that are now consecutive starting from next_batch_id.
        while let Some(Reverse((id, _))) = buffer.peek() {
            if *id != next_batch_id {
                break;
            }
            let Reverse((_, json)) = buffer.pop().unwrap();
            assert!(total_buffer_size >= json.len());
            total_buffer_size -= json.len();
            output.write_all(&json).unwrap();
            n_bytes_output.fetch_add(json.len(), Relaxed);
            next_batch_id += 1;
        }
    }

    if buffer.len() > 0 {
        log::error!("Missing output batch. Please file a bug report to the maintainers.");
        panic!("Missing output batch");
    }
}

pub fn run_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, input_file: &Path, output: impl Write + Send, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send> + 'static, metrics: &[Metric], n_aligners: usize, sort_output: bool) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&input_file).unwrap();

    let batch_size = 10_000_usize;
    let (work_send, work_recv) = crossbeam::channel::bounded::<QueryBatch>(n_aligners);
    let (results_send, results_recv) = crossbeam::channel::bounded::<(Vec<u8>, usize)>(n_aligners); // Json-formatted blocks of text, and the batch id that produced this json.

    let (progress_printer_quit_signal_send, progress_printer_quit_signal_recv) = crossbeam::channel::bounded::<()>(1);

    let n_bases_processed = AtomicUsize::new(0);
    let n_bytes_output = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        let parser_handle = scope.spawn(move || {
            let mut batch_id = 0_usize;
            let mut cur_batch = QueryBatch::new(batch_id);
            while let Some(q) = reader.read_next().unwrap() {
                cur_batch.seqs.push_record(q);
                if cur_batch.seqs.total_seq_len() >= batch_size {
                    work_send.send(cur_batch).unwrap();
                    batch_id += 1;
                    cur_batch = QueryBatch::new(batch_id);
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
                let mut metric_processors: Vec<Box<dyn PseudoalignmentMetricProcessor<CSS>>> = vec![];
                for metric in metrics {
                    metric_processors.push(create_metric_processor(*metric, index.get_set_storage().n_colors()));
                }
                while let Ok(batch) = work_recv_clone.recv() {
                    let (json, batch_id) = batch.process(index_ref, &mut aligner, n_bases_processed_ref, &mut metric_processors);
                    results_send_clone.send((json, batch_id)).unwrap();
                }
            });
            worker_handles.push(handle);
        }

        let outputter_handle = scope.spawn(|| {
            output_thread(results_recv, output, sort_output, &n_bytes_output);
        });

        let progress_printer_handle = scope.spawn(|| {
            let mut last_wakeup_time = std::time::Instant::now();
            let mut last_n_bases_processed = n_bases_processed.load(Relaxed);
            let mut last_n_bytes_output = n_bytes_output.load(Relaxed);
            let print_interval = std::time::Duration::from_secs(10);
            let start_time = std::time::Instant::now();
            loop {
                match progress_printer_quit_signal_recv.recv_timeout(print_interval) {
                    Ok(_) => break, // Received the quit signal
                    Err(RecvTimeoutError::Timeout) => { // Time to print
                        let n_in = n_bases_processed.load(Relaxed);
                        let n_out = n_bytes_output.load(Relaxed);
                        let t = last_wakeup_time.elapsed().as_secs_f64();
                        let in_throughput = (n_in - last_n_bases_processed) as f64 / t / 1000000 as f64;
                        let out_throughput = (n_out - last_n_bytes_output) as f64 / t / (1 << 20) as f64;
                        log::info!("Input {:.3} Mbases/s, output {:.3} MiB/s (total: {} bases, {} out)", in_throughput, out_throughput, n_in, human_bytes::human_bytes(n_out as f64));
                        last_n_bases_processed = n_in;
                        last_n_bytes_output = n_out;
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
            let total_n_out = n_bytes_output.load(Relaxed);
            let total_t = start_time.elapsed().as_secs_f64();
            let total_in_throughput = total_n as f64 / total_t / (1 << 20) as f64;
            let total_out_throughput = total_n_out as f64 / total_t / (1 << 20) as f64;
            log::info!("Total {} bases processed and {} output written in {:.3} seconds", total_n, human_bytes::human_bytes(total_n_out as f64), total_t);
            log::info!("Total throughput: input {:.3} Mbases/s, output {:.3} MiB/s", total_in_throughput, total_out_throughput);
        });
        
        parser_handle.join().unwrap(); // Wait for the parser to finish
        for h in worker_handles { h.join().unwrap() } // Wait for the workers to finish
        drop(results_send); // Signal that no more results will be pushed
        outputter_handle.join().unwrap(); // Wait for the outputter to finish
        progress_printer_quit_signal_send.send(()).unwrap(); // Interrupt the progress printer from sleep
        progress_printer_handle.join().unwrap();
    }); 
}