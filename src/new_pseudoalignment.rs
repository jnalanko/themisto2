use crossbeam::channel::{Receiver, Sender};
use jseqio::{record::Record, seq_db::SeqDB};
use rand_distr::num_traits::ConstOne;

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}};

trait CompatibilityCriterion {
    // The &mut self is to allow internal state containing reused buffers
    fn push_compatibility_set<CSS: ColorSetStorage>(&mut self, seq: &[u8], index: &CompactColexKmers<CSS>, out: &mut Vec<usize>);
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
    fn process<CSS: ColorSetStorage>(self, index: &CompactColexKmers<CSS>, cc: &mut impl CompatibilityCriterion) -> QueryResult {
        let mut result = QueryResult::new();
        let mut compat_set_buf = Vec::<usize>::new();
        for rec in self.seqs.iter() {
            compat_set_buf.clear();

            cc.push_compatibility_set(rec.seq, index, &mut compat_set_buf);
            result.push(&compat_set_buf, rec.name());

            let metric_vecs = self.compute_metrics(&rec.seq, &compat_set_buf, index);
            result.metrics.push(metric_vecs);
        }
        result
    }

    fn compute_metrics<CSS: ColorSetStorage>(&self, seq: &[u8], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<(Metric, Vec<usize>)>{
        let mut color_set_ids = Vec::<Option<usize>>::new();
        index.push_color_set_ids_to_buffer(seq, &mut color_set_ids);

        let mut ans = Vec::<(Metric, Vec<usize>)>::new();

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
    metrics: Vec<Vec<(Metric, Vec<usize>)>>
}

impl QueryResult {

    fn new() -> Self {
        Self { 
            query_names_concat: vec![], 
            query_names_starts: vec![0], 
            compatibility_class_concat: vec![], 
            compatibility_class_starts: vec![0], 
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