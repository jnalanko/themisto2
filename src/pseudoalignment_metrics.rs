use std::{cmp::min, ops::Range};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, util::for_each_run};

#[derive(Copy, Clone, Debug)]
pub enum Metric {
    KmerHits,
    BasesCovered,
    AlignmentLength,
    LongestMatchRun,
    ShortestGap
}

#[allow(clippy::manual_flatten)]
pub fn compute_kmer_hits_to_compatible_colors<CSS: ColorSetStorage>(color_set_ids: &[Option<usize>], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<usize> {
    let mut hits = vec![0; index.get_set_storage().n_colors()];
    for_each_run(color_set_ids, |run_range| {
        let run_len = run_range.len(); 
        assert!(run_len > 0);

        let first_id = color_set_ids[run_range.start];
        if let Some(set_id) = first_id {
            for color in index.set_id_to_set(set_id).iter() {
                hits[color] += run_len;
            }
        } // Runs of None are ignored
    });

    // Return only hits to compatible colors
    compatible_colors.iter().map(|&c| hits[c]).collect()
}

pub fn compute_bases_covered<CSS: ColorSetStorage>(color_set_ids: &[Option<usize>], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<usize> {
    let mut bases_covered = vec![0; index.get_set_storage().n_colors()];
    let mut end_of_last_covered_kmer = vec![0; index.get_set_storage().n_colors()]; // For each color. Exclusive end.

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
                let mut n_new_covered = min(first_kmer_end - end_of_last_covered_kmer[color], index.get_k());

                // Add new bases covered by the rest of the k-mers (1 base each)
                n_new_covered += run_len - 1;

                bases_covered[color] += n_new_covered;
                end_of_last_covered_kmer[color] = last_kmer_end;
            }
        } // Runs of None are ignored
    });
    
    // Return only counts to compatible colors
    compatible_colors.iter().map(|&c| bases_covered[c]).collect()
}

#[cfg(test)]
mod tests {

    use crate::{colex_colored_kmers, pseudoalignment_metrics::{compute_bases_covered, compute_kmer_hits_to_compatible_colors}, sparse_dense_storage::SparseDenseStorage};
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
        let hit_counts = compute_kmer_hits_to_compatible_colors(&cset_ids, &[1,2,3], &index);

        dbg!(&hit_counts);

        assert_eq!(hit_counts[0], s1.len()-k+1);
        assert_eq!(hit_counts[1], s2.len()-k+1);
        assert_eq!(hit_counts[2], s3.len()-k+1);
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
        let bases_covered = compute_bases_covered(&cset_ids, &[1,2,3], &index);

        dbg!(&bases_covered);

        assert_eq!(bases_covered[0], s1.len()-5);
        assert_eq!(bases_covered[1], s2.len());
        assert_eq!(bases_covered[2], s3.len()-3);

    }
}