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

// todo: can be much faster 
#[allow(clippy::manual_flatten)]
pub fn compute_kmer_hits_to_compatible_colors<CSS: ColorSetStorage>(color_set_ids: &[Option<usize>], compatible_colors: &[usize], index: &CompactColexKmers<CSS>) -> Vec<usize> {
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

/*
#[cfg(test)]
mod tests {
    use sbwt::SbwtConstructionAlgorithm;

    use crate::colex_colored_kmers::{self, ColexToColorSetMap};


    #[test]
    fn test_kmer_hits() {
        let s1 = b"AACTACGTACGTACGACATCGTACGATCGATTATGCTAGCTAGCTGAT".as_slice(); // "Random" sequence
        let s2 =           b"GTACGACATCGTACGATCGATTATGCTAGCTAGCTGAT".as_slice();
        let s3 =                           b"GTACGACATCGTACGATCGATT".as_slice();
        let s4 = b"AACTACGTACGTACGACATCGTACGATCGATTAT".as_slice();

        let inputs = vec![s1,s2,s3,s4];
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .algorithm(sbwt::BitPackedKmerSortingMem::new())
            .run_from_slices(&inputs);

        let index = colex_colored_kmers::CompactColexKmers::new(sbwt, lcs, colex_map, color_sets, color_names);
    }
}
    */