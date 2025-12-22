use std::cmp::min;

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}};

#[derive(Copy, Clone, Debug)]
pub enum Metric {
    KmerHits,
    BasesCovered,
    AlignmentLength,
    LongestMatchRun,
    ShortestGap
}

fn break_into_runs<T: Eq>(items: &[T]) -> Vec<&[T]> {
    let mut runs: Vec<&[T]> = vec![];
    if items.is_empty() {
        return runs;
    }

    let mut run_start = 0;
    for i in 1..items.len() {
        if items[i] != items[i-1] {
            runs.push(&items[run_start..i]);
            run_start = i;
        }
    }

    // Final run
    runs.push(&items[run_start..items.len()]);

    runs
}

#[allow(clippy::manual_flatten)]
pub fn compute_kmer_hits_to_compatible_colors<CSS: ColorSetStorage>(color_set_ids: &[Option<usize>], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<usize> {
    let mut hits = vec![0; index.get_set_storage().n_colors()];
    for color_set_run in break_into_runs(color_set_ids) {
        let run_len = color_set_run.len(); 
        match color_set_run.first() {
            None => panic!("Empty color set run"), // Should never happen
            Some(first_id) => {
                if let Some(set_id) = first_id {
                    for color in index.set_id_to_set(*set_id).iter() {
                        hits[color] += run_len;
                    }
                }
            }
        }
    }

    // Return only hits to compatible colors
    compatible_colors.iter().map(|&c| hits[c]).collect()
}

// todo: can be much faster 
pub fn compute_bases_covered<CSS: ColorSetStorage>(color_set_ids: &[Option<usize>], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<usize> {
    let mut bases_covered = vec![0; index.get_set_storage().n_colors()];
    let mut end_of_last_covered_kmer = vec![0; index.get_set_storage().n_colors()]; // For each color. Exclusive end.
    for (kmer_start, color_set_id_opt) in color_set_ids.iter().enumerate() {
        let kmer_end = kmer_start + index.get_k();
        if let Some(color_set_id) = color_set_id_opt {
            let color_set = index.set_id_to_set(*color_set_id);
            for color in color_set.iter() {
                let n_new_covered = min(kmer_end - end_of_last_covered_kmer[color], index.get_k());
                bases_covered[color] += n_new_covered;
                end_of_last_covered_kmer[color] = kmer_end;
            }
        }
    }
    
    // Return only counts to compatible colors
    compatible_colors.iter().map(|&c| bases_covered[c]).collect()
}

#[cfg(test)]
mod tests {

    use crate::{colex_colored_kmers, pseudoalignment_metrics::{compute_bases_covered, compute_kmer_hits_to_compatible_colors}, sparse_dense_storage::SparseDenseStorage};

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