use std::cmp::min;

use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};


// Returns (len, colex, position in sfs slice)
#[allow(clippy::manual_flatten)] // More readable
fn pick_finimizer(sfs_slice: &[Option<(usize, std::ops::Range<usize>)>]) -> (usize, usize, usize){
    // The finimizer is the shortest unique suffix, with ties broken by colex

    // The full k-mer should have an existing unique match
    assert!(sfs_slice.last().expect("Empty slice").as_ref().expect("Last SFS pos is None").1.len() == 1); 

    let mut best: (usize, usize) = (usize::MAX, usize::MAX); // (length, colex)
    let mut best_i: usize = usize::MAX;
    for (i,x) in sfs_slice.iter().enumerate() {
        if let Some((len, range)) = x { 
            if i + 1 >= *len && (*len, range.start) < best {
                best = min(best, (*len, range.start));
                best_i = i;
            }
        }
    }   

    if best == (usize::MAX, usize::MAX){
        dbg!(sfs_slice);
        panic!("Finimizer not found for kmer");
    }   

    (best.0, best.1, best_i)

}

fn explore(sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, f_colex: usize, f_len: usize, f_pos: usize) -> Vec<usize> {
    let si = StreamingIndex::new(sbwt, lcs);
    let k = sbwt.k();

    let mut kmer_colex_with_same_finimizer = Vec::<usize>::new();
    
    let mut initial_kmer = sbwt.access_kmer(f_colex);
    let mut initial_kmer_colex = f_colex;
    let mut dfs_stack = Vec::<(usize, Vec<u8>, usize)>::new(); // Depth, k-mer, colex
    dfs_stack.push((0, initial_kmer, initial_kmer_colex));

    //for _depth in 0..(k - f_len + 1) {
    while let Some((depth, kmer, colex)) = dfs_stack.pop() {
        if depth == k - f_len + 1 { continue } // Finimizer has falled out of the k-mer
        let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
        if pick_finimizer(&sfs).1 == f_colex { // Same finimizer
            kmer_colex_with_same_finimizer.push(colex);
        }
    }
    todo!();
}


pub fn finimizer_stats<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, n_threads: usize) {
    let sbwt = index.sbwt();
    let lcs = index.lcs();
    let si = StreamingIndex::new(sbwt, lcs);
    for colex in 0..sbwt.n_sets() {
        let kmer = sbwt.access_kmer(colex); // TODO: need to build select support for this
        let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
        let (f_len, f_colex, f_pos) = pick_finimizer(&sfs);
        let finimizer_kmer_set = explore(sbwt, f_colex, f_len, f_pos); 
    }
}