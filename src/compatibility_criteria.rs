use std::cmp::max;

fn indices_with_minimum_support(supports: &[usize], minimum_support: usize) -> Vec<usize> {
    supports.iter().enumerate().filter(|(_, s)| **s >= minimum_support).map(|(i,_)| i).collect()
}

// query_kmers can be e.g. either the number of k-mers in the query, or the number of relevant
// k-mers in the query.
pub fn basic_threshold_method(supports: &[usize], query_kmers: usize, mut min_support: usize, threshold: f64) -> Vec<usize> {
    min_support = max(min_support, (query_kmers as f64 * threshold) as usize);
    indices_with_minimum_support(supports, min_support)
}

/// Shared support counts include the unique support hits (todo: rename?)
pub fn unique_support_method(
    unique_supports: &[usize], 
    shared_supports: &[usize], 
    min_unique_support: usize, 
    min_shared_support: usize, 
    fraction_of_max: f64) -> Vec<usize> {

    // Try to determine based on unique support
    let max_unique_support = *unique_supports.iter().max().expect("Programming error: Empty unique support array");
    let unique_support_threshold = max(min_unique_support, (max_unique_support as f64 * fraction_of_max) as usize);
    if unique_supports.iter().any(|&s| s >= unique_support_threshold) {
        return indices_with_minimum_support(unique_supports, unique_support_threshold);
    }
        
    // Not enough evidence in unique support -> try shared support
    let max_shared_support = *shared_supports.iter().max().expect("Programming error: Empty shared support array");
    let shared_support_threshold = max(min_shared_support, (max_shared_support as f64 * fraction_of_max) as usize);

    if shared_supports.iter().any(|&s| s >= shared_support_threshold){
        return indices_with_minimum_support(shared_supports, shared_support_threshold);
    }

    vec![] // Not enough evidence for anything -> return the empty set
}

pub fn resolve_consensus_compatibility(
    compatibility_sets: &[&[usize]], 
    n_colors: usize, 
    min_unique_segments: usize, 
    min_shared_segments: usize,
    fraction_of_max: f64) -> Vec<usize> {

    let mut shared_supports = vec![0_usize; n_colors];
    let mut unique_supports = vec![0_usize; n_colors];
    for &set in compatibility_sets.iter() {
        for &x in set {
            shared_supports[x] += 1;
            if set.len() == 1 {
                unique_supports[x] += 1; 
            }
        }
    }

    unique_support_method(unique_supports.as_slice(), shared_supports.as_slice(), min_unique_segments, min_shared_segments, fraction_of_max)
}