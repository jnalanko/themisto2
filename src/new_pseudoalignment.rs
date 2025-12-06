use crossbeam::channel::{Receiver, Sender};
use jseqio::seq_db::SeqDB;

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

trait CompatibilityCriterion{
    // The &mut self is to allow internal state containing reused buffers
    fn push_compatibility_set(&mut self, seq: &[u8], out: &mut Vec<usize>);
}

struct QueryBatch {
    seqs: SeqDB,
    report_hit_counts: bool,
    report_bases_covered: bool,
    report_alignment_length: bool,
    report_longest_match_run: bool,
    report_shortest_gap: bool,
}

impl QueryBatch {
    fn process(self) -> QueryResult {
        todo!();
    }
}

struct QueryResult {
    compatibility_class_concat: Vec<usize>,
    compatibility_class_starts: Vec<usize>,
    hit_counts: Option<Vec<usize>>,
    bases_covered: Option<Vec<usize>>,
    alignment_length: Option<Vec<usize>>,
    longest_match_runs: Option<Vec<usize>>,
    shortest_gaps: Option<Vec<usize>>,
}

impl QueryResult {
    fn into_json(self, out: &mut Vec<u8>) {
        todo!();
    }
}

struct Worker<'a, CSS: ColorSetStorage> {
    index: &'a  CompactColexKmers<CSS>,
}

impl<'a, CSS: ColorSetStorage> Worker<'a, CSS> {
    fn run(&mut self, input: Receiver<QueryBatch>, output: &Sender<Vec<u8>>){
        while let Ok(batch) = input.recv() {
            let result = batch.process();
            let mut json_buf = Vec::<u8>::new();
            result.into_json(&mut json_buf);
            json_buf.push(b'\n');
            output.send(json_buf).unwrap();
        }
    }
}