use std::{cmp::min, iter::Map};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, util::for_each_run};

pub trait PseudoalignmentMetricProcessor<CSS: ColorSetStorage> {
    /// Returns pairs (color, metric value) for those colors that were found
    fn process(&mut self, color_set_ids: &[Option<usize>], index: &CompactColexKmers<CSS>) -> Vec<(usize, usize)>;

    fn metric_id(&self) -> Metric;
}

#[derive(Copy, Clone, Debug)]
pub enum Metric {
    KmerHits,
    BasesCovered,
    AlignmentLength,
    LongestMatchRun,
    ShortestGap
}

pub fn create_metric_processor<CSS: ColorSetStorage>(metric: Metric, n_colors: usize) -> Box<dyn PseudoalignmentMetricProcessor<CSS>> {
    match metric {
        Metric::KmerHits => Box::new(HitCountProcessor::new(n_colors)),
        Metric::BasesCovered => Box::new(BasesCoveredProcessor::new(n_colors)),
        Metric::AlignmentLength => todo!(),
        Metric::LongestMatchRun => todo!(),
        Metric::ShortestGap => todo!(),
    }
}

// An array of integers that tracks which indices have been touched.
pub struct NonzeroTrackingIntArray {
    data: Vec<usize>,
    nonzero_indices: Vec<usize>,
} 

impl NonzeroTrackingIntArray {
    pub fn get(&self, idx: usize) -> usize {
        self.data[idx]
    }

    pub fn add_positive_number(&mut self, idx: usize, value: usize) {
        debug_assert!(value > 0);
        if self.data[idx] == 0 {
            self.nonzero_indices.push(idx)
        }
        self.data[idx] += value;
    }

    pub fn reset(&mut self) {
        for &i in self.nonzero_indices.iter() {
            self.data[i] = 0;
        }
        self.nonzero_indices.clear();
    }

    // Iterates pairs (index, value)
    #[allow(dead_code)]
    pub fn iter_nonzero<'a>(&'a self) -> Map<std::slice::Iter<'a, usize>, impl FnMut(&'a usize) -> (usize, usize)> {
        self.nonzero_indices.iter().map(|&i| (i, self.data[i]))
    }

    pub fn new(len: usize) -> Self {
        Self {
            data: vec![0; len],
            nonzero_indices: vec![],
        }
    }
}

struct HitCountProcessor {
    hits: NonzeroTrackingIntArray, // Reused space between calls
}

impl HitCountProcessor {
    fn new(n_colors: usize) -> Self {
        Self { hits: NonzeroTrackingIntArray::new(n_colors) }
    }
}

impl<CSS: ColorSetStorage> PseudoalignmentMetricProcessor<CSS> for HitCountProcessor {
    fn process(&mut self, color_set_ids: &[Option<usize>], index: &CompactColexKmers<CSS>) -> Vec<(usize, usize)> {
        self.hits.reset();
        for_each_run(color_set_ids, |run_range| {
            let run_len = run_range.len(); 
            assert!(run_len > 0);

            let first_id = color_set_ids[run_range.start];
            if let Some(set_id) = first_id {
                for color in index.set_id_to_set(set_id).iter() {
                    self.hits.add_positive_number(color, run_len);
                }
            } // Runs of None are ignored
        });

        self.hits.iter_nonzero().collect() // Pair color, value
    }
    
    fn metric_id(&self) -> Metric {
        Metric::KmerHits
    }
    
}

struct BasesCoveredProcessor {

    // Reused space between calls

    // For each color.
    bases_covered: NonzeroTrackingIntArray,

    // For each color. Exclusive end.
    end_of_last_covered_kmer: NonzeroTrackingIntArray,
}

impl BasesCoveredProcessor {
    fn new(n_colors: usize) -> Self {
        Self { 
            bases_covered: NonzeroTrackingIntArray::new(n_colors),
            end_of_last_covered_kmer: NonzeroTrackingIntArray::new(n_colors),
        }
    }
}

impl<CSS: ColorSetStorage> PseudoalignmentMetricProcessor<CSS> for BasesCoveredProcessor {
    fn process(&mut self, color_set_ids: &[Option<usize>], index: &CompactColexKmers<CSS>) -> Vec<(usize, usize)> {

        self.bases_covered.reset();
        self.end_of_last_covered_kmer.reset();

        for_each_run(color_set_ids, |run_range| {
            let run_len = run_range.len(); 
            assert!(run_len > 0);
            let first_id = color_set_ids[run_range.start];

            let first_kmer_end = run_range.start + index.get_k();
            let last_kmer_end = run_range.end + index.get_k() - 1;
            if let Some(set_id) = first_id {
                let color_set = index.set_id_to_set(set_id);
                for color in color_set.iter() {

                    // New bases covered by the first k-mer
                    let mut n_new_covered = min(first_kmer_end - self.end_of_last_covered_kmer.get(color), index.get_k());

                    // Add new bases covered by the rest of the k-mers (1 base each)
                    n_new_covered += run_len - 1;

                    self.bases_covered.add_positive_number(color, n_new_covered);

                    // Set end of last covered k-mer. Here we use add_positive_integer because
                    // there is no method to set a value, because if there was, tracking nonzeros
                    // would be more complicated.
                    let old_end = self.end_of_last_covered_kmer.get(color);
                    self.end_of_last_covered_kmer.add_positive_number(color, last_kmer_end - old_end);
                }
            } // Runs of None are ignored
        });
        
        self.bases_covered.iter_nonzero().collect() // Pair color, value
    }
    
    fn metric_id(&self) -> Metric {
        Metric::BasesCovered 
    }

    
}


#[cfg(test)]
mod tests {

    use crate::{colex_colored_kmers, sparse_dense_storage::SparseDenseStorage};

    use super::*;

    #[test]
    fn test_kmer_hits() {
        let k = 8;
        let sample_distance = 3;
        let n_threads = 1;

        let s0 = b"AACTACGTACGTACGACATCGTACGATCGATTATGCTAGCTAGCTGAT".as_slice(); // "Random" sequence
        let s1 =           b"GTACGACATCGTACGATCGATTATGCTAGCTAGCTGAT".as_slice();
        let s2 =           b"GTACGACATCGTACGATCGATT".as_slice();
        let s3 = b"AACTACGTACGTACGACATCGTACGATCGATTAT".as_slice();

        let colored_seqs: Vec<(&[u8], usize)> = vec![(s0,0), (s1,1), (s2,2), (s3,3)];
        
        let index = colex_colored_kmers::CompactColexKmers::<SparseDenseStorage>::new_from_small_input(&colored_seqs, k, sample_distance, n_threads);

        let query = s0;
        let mut cset_ids = Vec::new();
        index.push_color_set_ids_to_buffer(query, &mut cset_ids);

        let mut processor = create_metric_processor(Metric::KmerHits, index.get_set_storage().n_colors());
        let mut hit_counts = processor.process(&cset_ids, &index);
        hit_counts.sort();

        dbg!(&hit_counts);

        assert_eq!(hit_counts[0], (0, s0.len()-k+1));
        assert_eq!(hit_counts[1], (1, s1.len()-k+1));
        assert_eq!(hit_counts[2], (2, s2.len()-k+1));
        assert_eq!(hit_counts[3], (3, s3.len()-k+1));
    }

    #[test]
    fn test_bases_covered() {
        let k = 8;
        let sample_distance = 3;
        let n_threads = 1;

        // s1 and s3 have substitutions with X
        let s0 = b"AACTACGTACGTACGACATCGTACGATCGATTATGCTAGCTAGCTGAT".as_slice(); // "Random" sequence
        let s1 =           b"GTACGACATCGTACGATCGXTTAXGCTAGCTAGCTGAT".as_slice();
        let s2 =           b"GTACGACATCGTACGATCGATT".as_slice();
        let s3 = b"AACTACGTAXXTACGACATCGTACXATCGATTAT".as_slice();

        let colored_seqs: Vec<(&[u8], usize)> = vec![(s0,0), (s1,1), (s2,2), (s3,3)];
        
        let index = colex_colored_kmers::CompactColexKmers::<SparseDenseStorage>::new_from_small_input(&colored_seqs, k, sample_distance, n_threads);

        let query = s0;
        let mut cset_ids = Vec::new();
        index.push_color_set_ids_to_buffer(query, &mut cset_ids);
        let mut processor = create_metric_processor(Metric::BasesCovered, index.get_set_storage().n_colors());
        let mut bases_covered = processor.process(&cset_ids, &index);
        bases_covered.sort();

        dbg!(&bases_covered);

        assert_eq!(bases_covered[0], (0, s0.len()));
        assert_eq!(bases_covered[1], (1, s1.len()-5));
        assert_eq!(bases_covered[2], (2, s2.len()));
        assert_eq!(bases_covered[3], (3, s3.len()-3));

    }
}