use std::{cmp::{max, min}, collections::HashMap, sync::{Arc, Mutex, atomic::{AtomicUsize, Ordering::Relaxed}}};

use rayon::iter::{IntoParallelIterator, ParallelBridge, ParallelIterator};
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView}};


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
        panic!("Finimizer not found for SFS slice");
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
    
    // Build the finimizer string. Note that this might not be a suffix
    // of any full k-mer, but a suffix of a dummy k-mer.
    let mut initial_suffix_match = sbwt.access_kmer(f_colex).to_vec();
    while *initial_suffix_match.first().unwrap() == b'$' {
        initial_suffix_match.remove(0);
    }

    //eprintln!("finimizer {}", String::from_utf8_lossy(&initial_suffix_match));

    let mut dfs_stack = Vec::<(usize, Vec<u8>, usize, bool)>::new(); // Depth, k-mer, colex, selected
    dfs_stack.push((0, initial_suffix_match, f_colex, false));


    while let Some((depth, suffix_match, colex, selected_before)) = dfs_stack.pop() {
        if depth == k - f_len + 1 { continue } // Finimizer has fallen out of the k-mer
        //eprintln!("suffix match {}", String::from_utf8_lossy(&suffix_match));

        let sfs = si.shortest_freq_bound_suffixes(&suffix_match, 1);
        let selected_here = pick_finimizer(&sfs).1 == f_colex;
        if selected_here { 
            if suffix_match.len() == k {
                kmer_colex_with_same_finimizer.push(colex);
            }
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
            let mut new_suffix_match = suffix_match.clone();
            new_suffix_match.push(c);
            if new_suffix_match.len() > k {
                assert!(new_suffix_match.len() == k+1);
                new_suffix_match.remove(0); // Pop front to get back to length k
            }

            if let Some(r) = sbwt.search(&new_suffix_match) {
                assert!(r.len() <= 1); // Still unique
                if r.len() > 0 { // Extension with c successful
                    dfs_stack.push((depth+1, new_suffix_match, r.start, selected_here));
                }
            }
        }
    }

    // Duplicate k-mers can happen if the DBG loops back to itself before the dfs depth limit
    kmer_colex_with_same_finimizer.sort();
    kmer_colex_with_same_finimizer.dedup();

    kmer_colex_with_same_finimizer
}


// Returns (n_correct, n_wrong)
fn evaluate_equivalence_class<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, class: &[usize]) -> (usize, usize) {
    let storage = index.get_set_storage();
    let mut finimizer_colors = storage.get_empty_set();
    for &colex in class {
        storage.union(&mut finimizer_colors, &index.colex_to_set(colex)); 
    }
    
    let mut n_correct = 0_usize;
    let mut n_wrong = 0_usize;
    for &colex in class {
        let view = index.colex_to_set(colex);
        if view.iter().eq(finimizer_colors.iter()) {
            n_correct += 1; // This k-mer has the same color set as the finimizer
        } else {
            n_wrong += 1;
        }
    }
    (n_correct, n_wrong)
}

// Requires select support
pub fn finimizer_stats<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, n_threads: usize, verify: bool) {
    let sbwt = index.sbwt();
    let lcs = index.lcs();
    let si = StreamingIndex::new(sbwt, lcs);


    // Data for the critical section
    let shared_visited_marks = bitvec::bitvec![0; sbwt.n_sets()];
    let mut n_correct_by_finimizer_len: Vec<usize> = vec![0; sbwt.k()+1]; 
    let mut n_wrong_by_finimizer_len: Vec<usize> = vec![0; sbwt.k()+1]; 
    let mut class_size_by_finimizer_len: Vec<usize> = vec![0; sbwt.k()+1]; 
    let mut n_finimizers_by_len: Vec<usize> = vec![0; sbwt.k()+1];
    let mut critical_data = Arc::new(Mutex::new((shared_visited_marks, n_correct_by_finimizer_len, n_wrong_by_finimizer_len, class_size_by_finimizer_len, n_finimizers_by_len)));

    let bar = indicatif::ProgressBar::new(sbwt.n_sets() as u64);

    /*
    let mut finimizer_to_kmers: Option<HashMap<usize, Vec<usize>>> = if verify { 
        Some(HashMap::new()) 
    } else { None };
     */

    //eprintln!("{}", String::from_utf8_lossy(&sbwt.access_kmer(364223)));

    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    pool.install(|| {
        (0..sbwt.n_sets()).into_par_iter().for_each(|colex| {
            if colex % 10000 == 0 {
                bar.inc(10000);
            }
            if critical_data.lock().unwrap().0[colex] { return } // Already visited
            let kmer = sbwt.access_kmer(colex);
            if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
                let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
                let (f_len, f_colex, _f_pos) = pick_finimizer(&sfs);
                let kmer_equivalence_class = finimizer_inverse_function(sbwt, lcs, f_colex, f_len); 
                assert!(kmer_equivalence_class.len() > 0); // At least the k-mer itself should be here
                let (n_correct, n_wrong) = evaluate_equivalence_class(index, &kmer_equivalence_class);
                assert_eq!(n_correct + n_wrong, kmer_equivalence_class.len());

                // Critical section: we must make sure that each class is counted only once
                let (visited, n_correct_vec, n_wrong_vec, class_size_by_finimizer_len, n_finimizers_by_len) = &mut *critical_data.lock().unwrap();
                // We need to check the visited bit again here because some other thread could have visited
                // this k-mer since we last checked the visited bit. The earlier check is redundant in the
                // sense that is does not change the result of the computation, but it does save unnecessary
                // computation.
                if !visited[colex] {
                    // Count this class
                    n_correct_vec[f_len] += n_correct;
                    n_wrong_vec[f_len] += n_wrong;
                    class_size_by_finimizer_len[f_len] += kmer_equivalence_class.len();
                    n_finimizers_by_len[f_len] += 1;
                    for p in kmer_equivalence_class.iter() {
                        visited.set(p, true);
                    }
                }
                // End of critical section

                /*
                if let Some(map) = finimizer_to_kmers.as_mut() {
                    for kmer_colex in kmer_equivalence_class.iter() {
                        let class = map.entry(f_colex).or_insert_with(Vec::new); // Create new if does not exist yet
                        class.push(kmer_colex);
                    }
                }
                */

            }
        })
    });
    bar.finish();

    let (_visited, n_correct_vec, n_wrong_vec, class_size_vec, n_finimizers_vec) = &*critical_data.lock().unwrap();
    let n_correct_total: usize = n_correct_vec.iter().sum();
    let n_wrong_total: usize = n_wrong_vec.iter().sum();
    assert_eq!(n_correct_total + n_wrong_total, sbwt.n_kmers()); // No double counting
    eprintln!("Fraction correct: {:.2}%", n_correct_total as f64 / (n_correct_total + n_wrong_total) as f64 * 100.0);
    for f_len in 0..=sbwt.k() {
        let n_correct: usize = n_correct_vec[f_len];
        let n_wrong: usize = n_wrong_vec[f_len];
        let mean_class_size = class_size_vec[f_len] as f64 / n_finimizers_vec[f_len] as f64;
        eprintln!("{}\t{:.5}\t{}\t{}\t{}\t{}", 
            f_len, 
            n_correct as f64 / (n_correct + n_wrong) as f64, 
            n_correct, 
            n_correct + n_wrong,
            mean_class_size,
            n_finimizers_vec[f_len]
        );
    }

    /*
    if let Some(finimizer_to_kmers) = finimizer_to_kmers {
        log::info!("Verifying that the classes contain every k-mer they are supposed to.");
        let mut n_kmers_checked = 0;
        let bar = indicatif::ProgressBar::new(sbwt.n_sets() as u64);
        for colex in 0..sbwt.n_sets() {
            bar.inc(1);
            let kmer = sbwt.access_kmer(colex);
            if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
                let sfs = si.shortest_freq_bound_suffixes(&kmer, 1);
                let (_f_len, f_colex, _f_pos) = pick_finimizer(&sfs);
                //eprintln!("{}, {} {}", _f_len, f_colex, _f_pos);
                let our_class = &finimizer_to_kmers[&f_colex];
                //eprintln!("{}, {:?} {:?}", String::from_utf8_lossy(&kmer), colex, our_class);
                assert!(our_class.contains(&colex));
                n_kmers_checked += 1;
            }
        }
        bar.finish();

        log::info!("Checking that classes are disjoint and have total size equal to the number of k-mers in the sbwt");
        let mut seen_colex_ranks = bitvec::bitvec![0; sbwt.n_sets()];
        let mut total_class_size = 0;
        for (_, class) in finimizer_to_kmers.iter() {
            for r in class.iter() {
                assert!(!seen_colex_ranks[r]);
                seen_colex_ranks.set(r, true);
            }
            total_class_size += class.len();
        }
        assert_eq!(n_kmers_checked, sbwt.n_kmers());
        assert_eq!(total_class_size, sbwt.n_kmers());
    }
    */
}