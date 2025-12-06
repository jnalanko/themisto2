use crossbeam::channel::{Receiver, Sender};
use jseqio::{record::Record, seq_db::SeqDB};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

trait CompatibilityCriterion {
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
    fn process<CSS: ColorSetStorage>(self, index: &CompactColexKmers<CSS>, cc: &mut impl CompatibilityCriterion) -> QueryResult {
        let mut result = QueryResult::new();
        let mut color_set_id_buf = Vec::<Option<usize>>::new();
        let mut compat_set_buf = Vec::<usize>::new();
        for rec in self.seqs.iter() {
            compat_set_buf.clear();
            color_set_id_buf.clear();

            let color_set_ids = index.push_color_set_ids_to_buffer(rec.seq, &mut color_set_id_buf);
            todo!();
            cc.push_compatibility_set(rec.seq, &mut compat_set_buf);
            result.push(&compat_set_buf, rec.name());
        }
        result
    }
}

struct QueryResult {
    query_names_concat: Vec<u8>,
    query_names_starts: Vec<usize>,

    compatibility_class_concat: Vec<usize>,
    compatibility_class_starts: Vec<usize>,

    // Optional metrics
    hit_counts: Option<Vec<usize>>,
    bases_covered: Option<Vec<usize>>,
    alignment_length: Option<Vec<usize>>,
    longest_match_runs: Option<Vec<usize>>,
    shortest_gaps: Option<Vec<usize>>,
}

impl QueryResult {

    fn new() -> Self {
        Self { 
            query_names_concat: vec![], 
            query_names_starts: vec![0], 
            compatibility_class_concat: vec![], 
            compatibility_class_starts: vec![0], 
            hit_counts: None, 
            bases_covered: None, 
            alignment_length: None, 
            longest_match_runs: None, 
            shortest_gaps: None,
        }
    }

    fn into_json(self, out: &mut Vec<u8>) {
        todo!();
    }

    fn push(&mut self, compat_set: &[usize], seq_name: &[u8]) {
        self.compatibility_class_concat.extend_from_slice(compat_set);
        self.compatibility_class_starts.push(self.compatibility_class_concat.len());

        self.query_names_concat.extend_from_slice(seq_name);
        self.query_names_starts.push(self.query_names_concat.len());
    }
}

struct Worker<'a, CSS: ColorSetStorage, CC: CompatibilityCriterion> {
    index: &'a  CompactColexKmers<CSS>,
    compatibility_criterion: CC,
}

impl<'a, CSS: ColorSetStorage, CC: CompatibilityCriterion> Worker<'a, CSS, CC> {
    fn run(&mut self, input: Receiver<QueryBatch>, output: &Sender<Vec<u8>>){
        while let Ok(batch) = input.recv() {
            let result = batch.process(&mut self.compatibility_criterion);
            let mut json_buf = Vec::<u8>::new();
            result.into_json(&mut json_buf);
            json_buf.push(b'\n');
            output.send(json_buf).unwrap();
        }
    }
}