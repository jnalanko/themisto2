use std::{cmp::{max, min}, collections::HashMap};

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
// Explore from a colex position that has a finimizer as a suffix, and return colex ranks of
// all k-mers that have the finimizer as their finimizer.
#[allow(clippy::collapsible_else_if)]
fn finimizer_inverse_function(sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, f_colex: usize, f_len: usize) -> Vec<usize> {
    let si = StreamingIndex::new(sbwt, lcs);
    let k = sbwt.k();

    let mut kmer_colex_with_same_finimizer = Vec::<usize>::new();
    
    let initial_kmer = sbwt.access_kmer(f_colex);
    let initial_kmer_colex = f_colex;

    let mut dfs_stack = Vec::<(usize, Vec<u8>, usize, bool)>::new(); // Depth, k-mer, colex, selected
    dfs_stack.push((0, initial_kmer, initial_kmer_colex, false));

    while let Some((depth, mut kmer, colex, selected_before)) = dfs_stack.pop() {
        if depth == k - f_len + 1 { continue } // Finimizer has falled out of the k-mer
        let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
        let selected_here = pick_finimizer(&sfs).1 == f_colex;
        if selected_here { 
            kmer_colex_with_same_finimizer.push(colex);
        } else { 
            if selected_before {
                // The finimizer with colex rank f_colex was selected in a previous
                // k-mer, but is not selected anymore. This means that there is now
                // a smaller finimizer in the same window, so we must wait for that
                // to fall our of the k-mer window before we can select f_colex again.
                // But if this happens, then f_colex is a suffix of the current k-mer,
                // which means we are back to where we started from because we have
                // unique finimizers.
                continue; 
            }
        }

        // Push out-neighbors to the dfs stack
        for c in [b'A', b'C', b'G', b'T'] {
            kmer.push(c);
            if let Some(r) = sbwt.search(&kmer[1..k+1]) {
                dfs_stack.push((depth+1, kmer[1..k+1].to_vec(), r.start, selected_here));
            }
            kmer.pop().unwrap(); 
        }
    }
    kmer_colex_with_same_finimizer
}


// Requires select support
pub fn finimizer_stats<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, n_threads: usize, verify: bool) {
    let sbwt = index.sbwt();
    let lcs = index.lcs();
    let si = StreamingIndex::new(sbwt, lcs);
    let mut visited = bitvec::bitvec![0; sbwt.n_sets()];
    let bar = indicatif::ProgressBar::new(sbwt.n_sets() as u64);

    let mut finimizer_to_kmers: Option<HashMap<usize, Vec<usize>>> = if verify { 
        Some(HashMap::new()) 
    } else { None };

    for colex in 0..sbwt.n_sets() {
        bar.inc(1);
        if visited[colex] { continue }
        let kmer = sbwt.access_kmer(colex);
        if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
            let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
            let (f_len, f_colex, _f_pos) = pick_finimizer(&sfs);
            let kmer_equivalence_class = finimizer_inverse_function(sbwt, lcs, f_colex, f_len); 
            for &p in kmer_equivalence_class.iter() {
                visited.set(p, true);
            }
            println!("{}", kmer_equivalence_class.len());

            for &kmer_colex in kmer_equivalence_class.iter() {
                if let Some(map) = finimizer_to_kmers.as_mut() {
                    let class = map.entry(f_colex).or_insert_with(Vec::new); // Create new if does not exist yet
                    class.push(kmer_colex);
                }
            }
        }
    }
    bar.finish();

    let bar = indicatif::ProgressBar::new(sbwt.n_sets() as u64);
    if let Some(finimizer_to_kmers) = finimizer_to_kmers {
        log::info!("Verifying that the classes contain every k-mer they are supposed to.");
        let mut n_kmers_checked = 0;
        for colex in 0..sbwt.n_sets() {
            bar.inc(1);
            let kmer = sbwt.access_kmer(colex); // TODO: need to build select support for this
            if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
                let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
                let (_f_len, f_colex, _f_pos) = pick_finimizer(&sfs);
                let our_class = &finimizer_to_kmers[&f_colex];
                assert!(our_class.contains(&colex));
                n_kmers_checked += 1;
            }
        }

        log::info!("Checking that classes are disjoint have total size equal to the number of k-mers in the sbwt");
        let mut seen_colex_ranks = bitvec::bitvec![0; sbwt.n_sets()];
        let mut total_class_size = 0;
        for (_, class) in finimizer_to_kmers.iter() {
            for &r in class.iter() {
                assert!(!seen_colex_ranks[r]);
                seen_colex_ranks.set(r, true);
            }
            total_class_size += class.len();
        }
        assert_eq!(n_kmers_checked, sbwt.n_kmers());
        assert_eq!(total_class_size, sbwt.n_kmers());
    }
    bar.finish();
}