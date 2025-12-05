use std::{io::{Read, Write}, ops::Range, path::Path};

use jseqio::seq_db::SeqDB;
use sbwt::SeqStream;

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}};

trait Pseudoaligner<CSS: ColorSetStorage> {

    // The &mut self is to access and modify thread-local buffers
    // owned by the algorithm.
    fn process<'a>(&'a mut self, seq: &[u8], index: &CompactColexKmers<CSS>) -> CSS::SetView<'a>;
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Denominator { // Options for the CLI
    All,
    Relevant,
    MaxHits,
}

#[derive(Clone)]
struct ThresholdPseudoaligner<CSS: ColorSetStorage> {
    counts: Vec<usize>,
    nonzero_count_indices: Vec<usize>,
    threshold: f64,
    denominator: Denominator,
    answer: CSS::OwnedSet, // Answer to the current query
}

impl<'a, CSS: ColorSetStorage> ThresholdPseudoaligner<CSS> {
    fn new(index: &'a CompactColexKmers<CSS>, threshold: f64, denominator: Denominator) -> Self {
        Self {
            counts: vec![0; index.get_set_storage().n_colors()],
            nonzero_count_indices: vec![],
            threshold,
            denominator,
            answer: index.get_set_storage().get_empty_set(),
        }
    }
}

impl<CSS: ColorSetStorage> Pseudoaligner<CSS> for ThresholdPseudoaligner<CSS> {
    fn process<'a>(&'a mut self, seq: &[u8], index: &CompactColexKmers<CSS>) -> CSS::SetView<'a> {
        todo!();
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

    fn process<CSS: ColorSetStorage>(self, aligner: &mut Box<dyn Pseudoaligner<CSS> + Send>) -> PseudoalignmentBatchResult {
        let mut result = PseudoalignmentBatchResult {
            concat: vec![],
            starts: vec![0],
            seq_ranks: self.seq_ranks,
        };

        for rec in self.seqs.iter() {
            let set = aligner.process(rec.seq);
            for color in set.iter() {
                result.concat.push(color);
            }
            result.starts.push(result.concat.len());
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
fn run_all_queries<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, mut queries: impl SeqStream + Send + 'static, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send>, mut output_callback: impl FnMut((usize, &[usize])) + Send, n_workers: usize) {

    let batch_size = 10_000_usize;
    let (work_send, work_recv) = crossbeam::channel::bounded::<PseudoalignmentBatch>(n_workers);
    let (results_send, results_recv) = crossbeam::channel::bounded::<PseudoalignmentBatchResult>(n_workers);

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
            let handle = scope.spawn(move || {
                while let Ok(batch) = work_recv_clone.recv() {
                    let result = batch.process(&mut aligner);
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
        
        parser_handle.join().unwrap(); // Wait for the parser to finish
        for h in worker_handles { h.join().unwrap() } // Wait for the workers to finish
        drop(results_send); // Signal that no more results will be pushed
        outputter_handle.join().unwrap(); // Wait for the outputter to finish
    }); 
}

pub fn run_pseudoalignment<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, input_file: &Path, mut output: impl Write + Send, create_new_aligner: impl Fn() -> Box<dyn Pseudoaligner<CSS> + Send>, n_aligners: usize) {
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