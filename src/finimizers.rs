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

fn create_finimizer_function<'b>(index: StreamingIndex<'b, SbwtIndex<SubsetMatrix>, LcsArray>) -> impl (for<'a> Fn(&'a [u8]) -> (usize, usize)) + 'b {
    move |kmer: &[u8]| {
        assert!(kmer.len() == index.k());
        let sfs = index.shortest_freq_bound_suffixes(kmer, 1);
        let (len, _colex, end) = pick_finimizer(&sfs);
        (end+1-len, len)
    }
}

// The finimizer function should return a pair (start, len)
#[allow(clippy::collapsible_else_if)]
fn find_kmer_class_of_minimizer<'b>(sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, minimizer: &[u8], minimizer_fn: impl for<'a> Fn(&'a [u8]) -> (usize, usize)) -> Vec<usize> {
    let k = sbwt.k();

    let mut kmer_colex_with_same_minimizer = Vec::<usize>::new();
    
    let minimizer_colex_range = sbwt.search(minimizer).unwrap();
    for minimizer_colex in minimizer_colex_range {
        let mut initial_suffix_match = sbwt.access_kmer(minimizer_colex).to_vec();
        while *initial_suffix_match.first().unwrap() == b'$' {
            initial_suffix_match.remove(0);
        }

        let mut dfs_stack = Vec::<(usize, Vec<u8>, usize, bool)>::new(); // Depth, k-mer, colex, selected
        dfs_stack.push((0, initial_suffix_match, minimizer_colex, false));


        while let Some((depth, suffix_match, kmer_colex, selected_before)) = dfs_stack.pop() {
            if depth == k - minimizer.len() + 1 { continue } // Finimizer has fallen out of the k-mer
            //eprintln!("suffix match {}", String::from_utf8_lossy(&suffix_match));

            let mut selected_here = false;
            if suffix_match.len() == k {
                let (f_start, f_len) = minimizer_fn(&suffix_match);
                let new_minimizer  = &suffix_match[f_start..f_start+f_len];
                if new_minimizer == minimizer {
                    kmer_colex_with_same_minimizer.push(kmer_colex);
                    selected_here = true;
                } else if selected_before {
                    // The minimizer with was selected in a previous
                    // k-mer, but is not selected anymore. This means that there is now
                    // a smaller ,inimizer in the same window, so we must wait for that
                    // to fall our of the k-mer window before we can select the old minimizer again.
                    // But if this happens, then the minimizer is a suffix of the current k-mer,
                    // so it will processed on some other round of the out for-loop.
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
    }

    // Duplicate k-mers can happen e.g. if the DBG loops back to itself before the dfs depth limit
    kmer_colex_with_same_minimizer.sort();
    kmer_colex_with_same_minimizer.dedup();
    kmer_colex_with_same_minimizer
}


// Returns (n_correct, n_wrong, mean jaccard index)
fn evaluate_equivalence_class<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, class: &[usize]) -> (usize, usize, f64) {
    let storage = index.get_set_storage();
    let mut finimizer_colors = storage.get_empty_set();
    for &colex in class {
        storage.union(&mut finimizer_colors, &index.colex_to_set(colex)); 
    }
    
    let mut n_correct = 0_usize;
    let mut n_wrong = 0_usize;
    let mut sum_jaccard = 0_f64;
    for &colex in class {
        let kmer_colors_view = index.colex_to_set(colex);

        if kmer_colors_view.iter().eq(finimizer_colors.iter()) {
            n_correct += 1; // This k-mer has the same color set as the finimizer
        } else {
            n_wrong += 1;
        }

        let mut intersection = finimizer_colors.clone();
        storage.intersect(&mut intersection, &kmer_colors_view);

        let mut union = finimizer_colors.clone();
        storage.union(&mut union, &kmer_colors_view);

        sum_jaccard += intersection.len() as f64 / union.len() as f64;

    }
    (n_correct, n_wrong, sum_jaccard / class.len() as f64)
}

// Requires select support
pub fn finimizer_stats<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, n_threads: usize, verify: bool) {
    let sbwt = index.sbwt();
    let lcs = index.lcs();
    let si = StreamingIndex::new(sbwt, lcs);

    // Data for the critical section
    struct CriticalSectionData {
        visited_marks: bitvec::vec::BitVec,
        n_correct_by_finimizer_len: Vec<usize>,
        n_wrong_by_finimizer_len: Vec<usize>,
        class_size_by_finimizer_len: Vec<usize>,
        sum_mean_jaccard_by_finimizer_len: Vec<f64>,
        n_finimizers_by_len: Vec<usize>,
    }
    
    let crit = CriticalSectionData {
        visited_marks: bitvec::bitvec![0; sbwt.n_sets()],
        n_correct_by_finimizer_len: vec![0; sbwt.k()+1],
        n_wrong_by_finimizer_len: vec![0; sbwt.k()+1],
        class_size_by_finimizer_len: vec![0; sbwt.k()+1],
        sum_mean_jaccard_by_finimizer_len: vec![0.0; sbwt.k()+1],
        n_finimizers_by_len: vec![0; sbwt.k()+1],
    };
    let crit = Arc::new(Mutex::new(crit));

    let bar = indicatif::ProgressBar::new(sbwt.n_sets() as u64);

    /*
    let mut finimizer_to_kmers: Option<HashMap<usize, Vec<usize>>> = if verify { 
        Some(HashMap::new()) 
    } else { None };
     */

    //eprintln!("{}", String::from_utf8_lossy(&sbwt.access_kmer(364223)));

    let si = StreamingIndex::new(sbwt, lcs);
    let finimizer_fn = create_finimizer_function(si);

    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    pool.install(|| {
        (0..sbwt.n_sets()).into_par_iter().for_each(|colex| {
            if colex % 10000 == 0 {
                bar.inc(10000);
            }
            if crit.lock().unwrap().visited_marks[colex] { return } // Already visited
            let kmer = sbwt.access_kmer(colex);
            if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
                let (f_start, f_len) = finimizer_fn(&kmer);
                let finimizer = &kmer[f_start..f_start+f_len];
                let kmer_equivalence_class = find_kmer_class_of_minimizer(sbwt, lcs, finimizer, &finimizer_fn); 

                assert!(kmer_equivalence_class.len() > 0); // At least the k-mer itself should be here
                let (n_correct, n_wrong, mean_jaccard) = evaluate_equivalence_class(index, &kmer_equivalence_class);
                assert_eq!(n_correct + n_wrong, kmer_equivalence_class.len());

                // Critical section: we must make sure that each class is counted only once
                let cr = &mut *crit.lock().unwrap();
                // We need to check the visited bit again here because some other thread could have visited
                // this k-mer since we last checked the visited bit. The earlier check is redundant in the
                // sense that is does not change the result of the computation, but it does save unnecessary
                // computation.
                if !cr.visited_marks[colex] {
                    // Count this class
                    cr.n_correct_by_finimizer_len[f_len] += n_correct;
                    cr.n_wrong_by_finimizer_len[f_len] += n_wrong;
                    cr.class_size_by_finimizer_len[f_len] += kmer_equivalence_class.len();
                    cr.sum_mean_jaccard_by_finimizer_len[f_len] += mean_jaccard;
                    cr.n_finimizers_by_len[f_len] += 1;
                    for p in kmer_equivalence_class.iter() {
                        cr.visited_marks.set(p, true);
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

    let cr = &*crit.lock().unwrap();
    let n_correct_total: usize = cr.n_correct_by_finimizer_len.iter().sum();
    let n_wrong_total: usize = cr.n_wrong_by_finimizer_len.iter().sum();
    assert_eq!(n_correct_total + n_wrong_total, sbwt.n_kmers()); // No double counting
    eprintln!("Fraction correct: {:.2}%", n_correct_total as f64 / (n_correct_total + n_wrong_total) as f64 * 100.0);
    for f_len in 0..=sbwt.k() {
        let n_correct: usize = cr.n_correct_by_finimizer_len[f_len];
        let n_wrong: usize = cr.n_wrong_by_finimizer_len[f_len];
        let mean_class_size = cr.class_size_by_finimizer_len[f_len] as f64 / cr.n_finimizers_by_len[f_len] as f64;
        let mean_jaccard = cr.sum_mean_jaccard_by_finimizer_len[f_len] / cr.n_finimizers_by_len[f_len] as f64;
        println!("{}\t{:.5}\t{}\t{}\t{:.5}\t{:.5}\t{}", 
            f_len, 
            n_correct as f64 / (n_correct + n_wrong) as f64, 
            n_correct, 
            n_correct + n_wrong,
            mean_class_size,
            mean_jaccard,
            cr.n_finimizers_by_len[f_len]
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