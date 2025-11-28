use std::{cmp::min, collections::HashSet, marker::PhantomData, sync::{Arc, Mutex, atomic::{AtomicU8, AtomicU64, AtomicUsize, Ordering::{Acquire, Relaxed, Release, SeqCst}}}};

use rayon::{iter::{IntoParallelIterator, ParallelIterator}, slice::ParallelSliceMut};
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{ops::{BitVec, Rank}, raw_vector::AccessRaw};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView}, int_vec::AtomicCompactIntVec, set_of_sets_construction::SetElement, sparse_dense_storage::SparseDenseStorage};


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

fn create_minimizer_function<'b>(m: usize) -> impl (for<'a> Fn(&'a [u8]) -> (usize, usize)) + 'b {
    move |kmer: &[u8]| {
        let mut minimizer = &kmer[0..m];
        let mut min_pos = 0;
        for j in 1 .. (kmer.len() as i64) - (m as i64) + 1 {
            let j = j as usize;
            if kmer[j..j+m] < *minimizer {
                minimizer = &kmer[j..j+m];
                min_pos = j;
            }
        }
        (min_pos, m)
    }
}

// The finimizer function should return a pair (start, len)
#[allow(clippy::collapsible_else_if)]
fn find_kmer_class_of_minimizer<'b>(sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, minimizer: &[u8], minimizer_fn: impl for<'a> Fn(&'a [u8]) -> (usize, usize)) -> Vec<usize> {
    let k = sbwt.k();

    // Duplicate k-mers can happen e.g. if the DBG loops back to itself before the dfs depth limit
    let mut kmer_colex_with_same_minimizer = HashSet::<usize>::new();
    
    let minimizer_colex_range = sbwt.search(minimizer).unwrap();
    for minimizer_colex in minimizer_colex_range {
        let mut initial_suffix_match = sbwt.access_kmer(minimizer_colex).to_vec();
        while *initial_suffix_match.first().unwrap() == b'$' {
            initial_suffix_match.remove(0);
        }

        let mut dfs_stack = Vec::<(usize, Vec<u8>, usize, bool)>::new(); // Depth, k-mer, colex, selected
        dfs_stack.push((0, initial_suffix_match, minimizer_colex, false));

        let mut visited = HashSet::<usize>::new();
        while let Some((depth, suffix_match, kmer_colex, selected_before)) = dfs_stack.pop() {
            if depth == k - minimizer.len() + 1 { continue } // Finimizer has fallen out of the k-mer
            if visited.contains(&kmer_colex) { continue } // Already been here
            visited.insert(kmer_colex);

            let mut selected_here = false;
            if suffix_match.len() == k {
                let (f_start, f_len) = minimizer_fn(&suffix_match);
                let new_minimizer  = &suffix_match[f_start..f_start+f_len];
                if new_minimizer == minimizer {
                    kmer_colex_with_same_minimizer.insert(kmer_colex);
                    selected_here = true;
                    if kmer_colex_with_same_minimizer.len() % 1000000 == 0 {
                        log::warn!("Minimizer class of size at least {}", kmer_colex_with_same_minimizer.len());
                    }
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

            // Find the suffix group start of the k-mer. The outgoing DBG-edge labels are listed there.
            let mut suffix_group_start = kmer_colex;
            while suffix_group_start > 0 && lcs.access(suffix_group_start) == k - 1 {
                suffix_group_start -= 1;
            }
            for (c_idx, &c) in [b'A', b'C', b'G', b'T'].iter().enumerate() {
                if sbwt.sbwt().set_contains(suffix_group_start, c_idx as u8) {
                    // c is an outgoing edge from here
                    let mut new_suffix_match = suffix_match.clone();
                    new_suffix_match.push(c);
                    if new_suffix_match.len() > k {
                        assert!(new_suffix_match.len() == k+1);
                        new_suffix_match.remove(0); // Pop front to get back to length k
                    }

                    let new_colex = sbwt.lf_step(suffix_group_start, c_idx);
                    dfs_stack.push((depth+1, new_suffix_match, new_colex, selected_here));
                }
            }
        }
    }

    let mut ret: Vec<usize> = kmer_colex_with_same_minimizer.into_iter().collect();
    ret.sort();
    ret 
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

pub enum MinimizerType {
    Finimizer,
    Minimizer(usize), // The usize is the minimizer length
}

struct ElementGenerator<'a, 'c, F: for<'b> Fn(&'b [u8]) -> (usize, usize) + Sync + Send> {
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    minimizer_marks: &'c simple_sds_sbwt::bit_vector::BitVector,
    minimizer_fn: F,
} 

impl<'a,'c, F: for<'b> Fn(&'b [u8]) -> (usize, usize) + Sync + Send> crate::set_of_sets_construction::ParallelElementGenerator for ElementGenerator<'a, 'c, F> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        let bar = indicatif::ProgressBar::new(self.sbwt.n_kmers() as u64);
        let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        pool.install(|| {
            (0..self.sbwt.n_sets()).into_par_iter().for_each(|kmer_colex| {
                if kmer_colex % 1000000 == 0 {
                    bar.inc(1000000);
                }
                let kmer = self.sbwt.access_kmer(kmer_colex);
                if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
                    let (f_start, f_len) = (self.minimizer_fn)(kmer.as_slice());
                    let finimizer = &kmer[f_start..f_start+f_len];
                    let minimizer_colex = self.sbwt.search(finimizer).unwrap().start;
                    let minimizer_rank = self.minimizer_marks.rank(minimizer_colex);
                    callback(SetElement { set_id: minimizer_rank, color: kmer_colex });
                }
            });
        });
        bar.finish();
    }

    fn set_filter(&mut self, filter: simple_sds_sbwt::bit_vector::BitVector) {
        unimplemented!()
    }
}

pub fn generic_minimizer_stats_rewrite<CSS: ColorSetStorage + Sync, F: for<'a> Fn(&'a [u8]) -> (usize, usize) + Sync + Send> (index: &CompactColexKmers<CSS>, n_threads: usize, minimizer_fn: F) {
    let sbwt = index.sbwt();

    log::info!("Marking *inimizers");
    let minimizer_marks = AtomicBitmap::new(sbwt.n_sets());
    let bar = indicatif::ProgressBar::new(sbwt.n_kmers() as u64);
    (0..sbwt.n_sets()).into_par_iter().for_each(|kmer_colex| {
        if kmer_colex % 1000000 == 0 {
            bar.inc(1000000);
        }
        let kmer = sbwt.access_kmer(kmer_colex);
        if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
            let (f_start, f_len) = minimizer_fn(&kmer);
            let finimizer = &kmer[f_start..f_start+f_len];
            let minimizer_colex = sbwt.search(finimizer).unwrap().start;
            minimizer_marks.set(minimizer_colex, true);
        }
    });
    bar.finish();

    log::info!("Building rank support for *inimizer marks");
    let minimizer_marks = minimizer_marks.into_bitvec();
    let mut rv = simple_sds_sbwt::raw_vector::RawVector::with_len(minimizer_marks.len(), false);
    for b in minimizer_marks.iter_ones() {
        rv.set_bit(b, true);
    }
    let mut minimizer_marks = simple_sds_sbwt::bit_vector::BitVector::from(rv);
    minimizer_marks.enable_rank();
    let n_minimizers = minimizer_marks.rank(minimizer_marks.len());

    log::info!("Computing class sizes and *inimizer lengths");
    let class_sizes = (0..n_minimizers).map(|_| AtomicUsize::new(0)).collect::<Vec<_>>();
    let minimizer_lengths = (0..n_minimizers).map(|_| AtomicU8::new(0)).collect::<Vec<_>>();
    let bar = indicatif::ProgressBar::new(sbwt.n_kmers() as u64);
    (0..sbwt.n_sets()).into_par_iter().for_each(|kmer_colex| {
        if kmer_colex % 1000000 == 0 {
            bar.inc(1000000);
        }
        let kmer = sbwt.access_kmer(kmer_colex);
        if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
            let (f_start, f_len) = minimizer_fn(&kmer);
            let minimizer = &kmer[f_start..f_start+f_len];
            let minimizer_colex = sbwt.search(minimizer).unwrap().start;
            let minimizer_rank = minimizer_marks.rank(minimizer_colex);
            class_sizes[minimizer_rank].fetch_add(1, Release);
            minimizer_lengths[minimizer_rank].store(minimizer.len() as u8, Relaxed);
        }
    });
    bar.finish();
    let class_sizes = class_sizes.into_iter().map(|x| x.load(Relaxed)).collect::<Vec::<usize>>();

    let element_gen = ElementGenerator{
        sbwt,
        minimizer_marks: &minimizer_marks,
        minimizer_fn: Box::new(minimizer_fn),
    };
    
    log::info!("Storing minimizer kmer classes");
    let kmer_class_storage = SparseDenseStorage::new_parallel(element_gen, sbwt.n_kmers(), &class_sizes, n_threads);

    log::info!("Computing stats");
    struct Stats{
        n_correct_by_finimizer_len: Vec<usize>,
        n_wrong_by_finimizer_len: Vec<usize>,
        class_size_by_finimizer_len: Vec<usize>,
        sum_mean_jaccard_by_finimizer_len: Vec<f64>,
        n_finimizers_by_len: Vec<usize>,
    }
    
    let mut stats = Stats {
        n_correct_by_finimizer_len: vec![0; sbwt.k()+1],
        n_wrong_by_finimizer_len: vec![0; sbwt.k()+1],
        class_size_by_finimizer_len: vec![0; sbwt.k()+1],
        sum_mean_jaccard_by_finimizer_len: vec![0.0; sbwt.k()+1],
        n_finimizers_by_len: vec![0; sbwt.k()+1],
    };

    let bar = indicatif::ProgressBar::new(sbwt.n_kmers() as u64);
    for minimizer_id in 0..n_minimizers {
        let class: Vec<usize> = kmer_class_storage.get_set_view(minimizer_id).iter().collect();
        let (n_correct, n_wrong, mean_jaccard) = evaluate_equivalence_class(index, &class);
        assert_eq!(n_correct + n_wrong, class.len());

        let f_len = minimizer_lengths[minimizer_id].load(Relaxed) as usize;
        stats.n_correct_by_finimizer_len[f_len] += n_correct;
        stats.n_wrong_by_finimizer_len[f_len] += n_wrong;
        stats.class_size_by_finimizer_len[f_len] += class.len();
        stats.sum_mean_jaccard_by_finimizer_len[f_len] += mean_jaccard;
        stats.n_finimizers_by_len[f_len] += 1;
        bar.inc(1);
    }
    bar.finish();

    log::info!("Printing");
    let n_correct_total: usize = stats.n_correct_by_finimizer_len.iter().sum();
    let n_wrong_total: usize = stats.n_wrong_by_finimizer_len.iter().sum();
    assert_eq!(n_correct_total + n_wrong_total, sbwt.n_kmers()); // No double counting
    log::info!("Fraction correct: {:.2}%", n_correct_total as f64 / (n_correct_total + n_wrong_total) as f64 * 100.0);
    for f_len in 0..=sbwt.k() {
        let n_correct: usize = stats.n_correct_by_finimizer_len[f_len];
        let n_wrong: usize = stats.n_wrong_by_finimizer_len[f_len];
        let mean_class_size = stats.class_size_by_finimizer_len[f_len] as f64 / stats.n_finimizers_by_len[f_len] as f64;
        let mean_jaccard = stats.sum_mean_jaccard_by_finimizer_len[f_len] / stats.n_finimizers_by_len[f_len] as f64;
        println!("{}\t{:.5}\t{}\t{}\t{:.5}\t{:.5}\t{}", 
            f_len, 
            n_correct as f64 / (n_correct + n_wrong) as f64, 
            n_correct, 
            n_correct + n_wrong,
            mean_class_size,
            mean_jaccard,
            stats.n_finimizers_by_len[f_len]
        );
    }
}

pub fn generic_minimizer_stats<CSS: ColorSetStorage + Sync, F: for<'a> Fn(&'a [u8]) -> (usize, usize) + Sync + Send> (index: &CompactColexKmers<CSS>, n_threads: usize, minimizer_fn: F) {
    let sbwt = index.sbwt();
    let lcs = index.lcs();

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

    // When some class is processed, all colex positions of k-mers in the class
    // are marked here. If a bit is set, we know we do not have to process
    // that k-mer anymore. If a bit is not set, it might mean that it's in the
    // class of some k-mer that is currently being processed. So we *might* still
    // need to process it. The visited marks in the critical section will then
    // tell us whether we should record the result, or move on because some other
    // thread already recorded this work.
    let atomic_filter_bitmap = AtomicBitmap::new(sbwt.n_sets());

    let bar = indicatif::ProgressBar::new(sbwt.n_kmers() as u64);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    pool.install(|| {
        (0..sbwt.n_sets()).into_par_iter().for_each(|colex| {
            if atomic_filter_bitmap.get(colex) { return }
            let kmer = sbwt.access_kmer(colex);
            if kmer.iter().all(|&c| c != b'$') { // Not a dummy k-mer
                let (f_start, f_len) = minimizer_fn(&kmer);
                let finimizer = &kmer[f_start..f_start+f_len];
                let kmer_equivalence_class = find_kmer_class_of_minimizer(sbwt, lcs, finimizer, &minimizer_fn); 

                assert!(kmer_equivalence_class.len() > 0); // At least the k-mer itself should be here
                let (n_correct, n_wrong, mean_jaccard) = evaluate_equivalence_class(index, &kmer_equivalence_class);
                assert_eq!(n_correct + n_wrong, kmer_equivalence_class.len());

                // Signal to other threads that these k-mers do not need to be processed anymore
                for p in kmer_equivalence_class.iter() {
                    atomic_filter_bitmap.set(p, true);
                }

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
                        //assert!(!cr.visited_marks[p]);
                        cr.visited_marks.set(p, true);
                    }
                    bar.inc(kmer_equivalence_class.len() as u64);
                }
                // End of critical section
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
}

// Requires select support
pub fn minimizer_stats<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, n_threads: usize, minimizer_type: MinimizerType) {
    let sbwt = index.sbwt();
    let lcs = index.lcs();
    let si = StreamingIndex::new(sbwt, lcs);
    match minimizer_type {
        MinimizerType::Finimizer => {
            let finimizer_fn = create_finimizer_function(si);
            generic_minimizer_stats_rewrite(index, n_threads, finimizer_fn);
        },
        MinimizerType::Minimizer(m) => {
            let minimizer_fn = create_minimizer_function(m);
            generic_minimizer_stats_rewrite(index, n_threads, minimizer_fn);
        },
    }
}
