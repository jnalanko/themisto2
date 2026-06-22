use bitvec::order::Lsb0;
use crossbeam::channel::bounded;
use jseqio::reverse_complement;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use sbwt::dbg::Dbg;
use sbwt::reverse_complement_in_place;
use sbwt::{dbg::Node, SbwtIndex, SubsetMatrix};
use simple_sds_sbwt::ops::BitVec;
use std::io::Write;
use std::ops::Range;
use std::sync::Mutex;

use crate::colex_colored_kmers::ColexToColorSetMap;
use crate::coloring_interface::{ColorSetStorage, ColorSetView};


// Bit vectors accumulated during GFA segment export (phase 1) and consumed during link
// computation (phase 2).  Each vector is indexed by SBWT colex rank.
//   is_fw_first: set for the first forward k-mer of every output segment
//   is_fw_last:  set for the last forward k-mer of every output segment (= segment name)
//   is_rc_last:  set for rc_colex[last] = RC(fw_last) of every segment, i.e. the entry
//                point of the segment's reverse-complement traversal
struct GfaBitVecs {
    is_fw_first: bitvec::vec::BitVec<usize, Lsb0>,
    is_fw_last: bitvec::vec::BitVec<usize, Lsb0>,
    is_rc_first: bitvec::vec::BitVec<usize, Lsb0>,
    is_rc_last: bitvec::vec::BitVec<usize, Lsb0>,
}


/// Same format as [crate::index_import].
/// Select support must be built before calling this!
/// Returns the number of exported unitigs
pub fn export_colored_unitigs(sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<SubsetMatrix>, map: &ColexToColorSetMap, unitigs_out: impl Write + Sync + Send, n_threads: usize) -> usize {
    let k = sbwt.k();
    let sbwt_len = sbwt.n_sets();
    log::info!("Exporting unitigs");
    let n_unitigs = export_canonical_unitigs_with_shared_color_set(dbg, map, sbwt, k, sbwt_len, unitigs_out, n_threads);

    #[allow(clippy::let_and_return)] // Is clearer
    n_unitigs
}

// Like export_canonical_unitigs_with_shared_color_set but writes GFA S-lines and
// accumulates the three bit vectors needed for link computation.
fn export_canonical_unitigs_for_gfa(mut gfa_out: impl Write + Sync + Send, dbg_ref: &Dbg<SubsetMatrix>, map: &ColexToColorSetMap, sbwt: &SbwtIndex<SubsetMatrix>, k: usize, sbwt_len: usize, n_threads: usize) -> (usize, GfaBitVecs) {

    log::info!("Computing unitigs");
    let n_unitig_searches = std::sync::atomic::AtomicUsize::new(0);
    let n_unitig_searches_ref = &n_unitig_searches;

    let bar = indicatif::ProgressBar::new(sbwt_len as u64);
    let (n_segments, bvs) = std::thread::scope(|scope| {

        let (worker_out, collector_in) = bounded::<(Vec<usize>, Vec<usize>, Vec<u8>, Vec<Range<usize>>, Vec<usize>)>(n_threads);

        let mut worker_handles = Vec::<_>::new();
        let bar_ref = &bar;
        for thread_id in 0..n_threads {
            let worker_out_clone = worker_out.clone();
            let handle = scope.spawn(move || {
                let mut colex = thread_id;
                while colex < sbwt_len {
                    let v = Node { id: colex };
                    if !dbg_ref.is_dummy_colex_position(colex) && dbg_ref.is_first_kmer_of_unitig(v) {
                        n_unitig_searches_ref.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        worker_out_clone.send(search_unitig_from(v, dbg_ref, k, map, sbwt)).unwrap();
                    }
                    colex += n_threads;
                    if ((colex - thread_id)/n_threads) % 10000 == 0 {
                        bar_ref.inc(10000);
                    }
                }
                log::info!("Thread {} finished", thread_id);
            });
            worker_handles.push(handle);
        }

        let collector_handle = scope.spawn(move || {
            let mut unitig_id = 0_usize;
            let mut visited = bitvec::bitvec![usize, Lsb0; 0; sbwt_len];
            let mut gfa_bvs = Some(GfaBitVecs {
                is_fw_first: bitvec::bitvec![usize, Lsb0; 0; sbwt_len],
                is_fw_last:  bitvec::bitvec![usize, Lsb0; 0; sbwt_len],
                is_rc_first:  bitvec::bitvec![usize, Lsb0; 0; sbwt_len],
                is_rc_last:  bitvec::bitvec![usize, Lsb0; 0; sbwt_len],
            });

            while let Ok((fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subunitig_color_set_ids)) = collector_in.recv() {
                visit_and_output_kmers(&unitig_string, &subunitig_kmer_ranges, &subunitig_color_set_ids, &fw_colex, &rc_colex, &mut visited, &mut gfa_out, &mut unitig_id, &mut gfa_bvs, k);
            }

            log::info!("Processing remaining cyclic unitigs");
            let n_acyclic = unitig_id;
            let mut colex = 0_usize;
            while colex < visited.len() {
                colex = match visited[colex..].first_zero() {
                    Some(i) => colex + i,
                    None => break,
                };
                if !dbg_ref.is_dummy_colex_position(colex) {
                    let (fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subunitig_color_set_ids)
                        = search_unitig_from(Node { id: colex }, dbg_ref, k, map, sbwt);
                    assert!(unitig_string[..k-1] == unitig_string[unitig_string.len()-(k-1)..]);
                    visit_and_output_kmers(&unitig_string, &subunitig_kmer_ranges, &subunitig_color_set_ids, &fw_colex, &rc_colex, &mut visited, &mut gfa_out, &mut unitig_id, &mut gfa_bvs, k);
                }
                colex += 1;
            }
            gfa_out.flush().unwrap();
            log::info!("Found {} cyclic unitigs", unitig_id - n_acyclic);
            (unitig_id, gfa_bvs.unwrap())
        });

        for h in worker_handles {
            h.join().unwrap();
        }
        drop(worker_out);
        collector_handle.join().unwrap()
    });
    bar.finish();

    log::info!("Wrote {} segments", n_segments);
    (n_segments, bvs)
}

/// Export subunitigs as GFA 1.0.  Segment names are "u{colex}" where colex is the
/// colex rank of the last forward k-mer of the segment.  Each S-line carries a
/// "cs:i:{id}" tag with the color-set id.  L-lines encode all inter-segment
/// overlaps of type ++, +- of length k-1.  Also writes the metadata and color sets files in the
/// same format as [Self::export_colored_unitigs].
/// Select support must be built before calling this.
/// Returns the number of exported unitigs.
pub fn export_gfa(sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<SubsetMatrix>, map: &ColexToColorSetMap, mut gfa_out: impl Write + Sync + Send, n_threads: usize) -> usize {
    let k = sbwt.k();

    writeln!(gfa_out, "H\tVN:Z:1.0").unwrap();
    
    log::info!("Exporting unitigs (GFA format)");
    log::info!("Computing GFA segments");
    let (n_segments, bvs) = export_canonical_unitigs_for_gfa(&mut gfa_out, &dbg, map, sbwt, k, sbwt.n_sets(), n_threads);

    log::info!("Computing GFA links");

    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    pool.install(|| {

        // Put the output writer into a mutex in order to be
        // able to share it between threads.
        let gfa_out_mutex = Mutex::new(gfa_out);
        let gfa_out_mutex_ref = &gfa_out_mutex;

        (0..sbwt.n_sets()).into_par_iter().for_each(|colex| {
            if dbg.is_dummy_colex_position(colex) || (!bvs.is_fw_last[colex] && !bvs.is_rc_last[colex]) {
                return;
            }

            // Todo: make these into thread-local state so they are not
            // reallocated every time.
            let mut out_nbuf  = Vec::<(Node, u8)>::new();
            let mut walk_nbuf = Vec::<(Node, u8)>::new();
            let mut kmer_buf  = Vec::<u8>::new();

            let left_orientation = if bvs.is_fw_last[colex] { '+' } else { '-' };
            let left_name = if left_orientation == '+' { 
                colex // Colex rank of the last k-mer
            } else {
                // Now the name is the colex of the last k-mer of the other direction
                kmer_buf.clear();
                sbwt.push_kmer_to_vec(colex, &mut kmer_buf);
                reverse_complement_in_place(&mut kmer_buf);
                let other = sbwt.search(&kmer_buf)
                    .expect("RC of an indexed k-mer must itself be indexed")
                    .start;

                // Walk to the opposite end
                let mut v = other;
                while !bvs.is_fw_last[v] {
                    walk_nbuf.clear();
                    dbg.push_out_neighbors(Node { id: v }, &mut walk_nbuf);
                    v = walk_nbuf[0].0.id;
                }
                v
            };

            dbg.push_out_neighbors(Node { id: colex }, &mut out_nbuf);
            for i in 0..out_nbuf.len() {
                let nc = out_nbuf[i].0.id;

                if bvs.is_fw_first[nc] {
                    // ++ link: walk forward from nc within its segment to find fw_last_B
                    let mut v = nc;
                    while !bvs.is_fw_last[v] {
                        walk_nbuf.clear();
                        dbg.push_out_neighbors(Node { id: v }, &mut walk_nbuf);
                        v = walk_nbuf[0].0.id;
                    }
                    let out = &mut (*gfa_out_mutex_ref.lock().unwrap());

                    writeln!(out, "L\tu{}\t{}\tu{}\t+\t{}M", left_name, left_orientation, v, k - 1).unwrap();
                }

                if bvs.is_rc_first[nc] {
                    kmer_buf.clear();
                    sbwt.push_kmer_to_vec(nc, &mut kmer_buf);
                    reverse_complement_in_place(&mut kmer_buf);
                    let fw_last = sbwt.search(&kmer_buf) // Name of the target unitig
                        .expect("RC of an indexed k-mer must itself be indexed")
                        .start;

                    let out = &mut (*gfa_out_mutex_ref.lock().unwrap());
                    writeln!(out, "L\tu{}\t{}\tu{}\t-\t{}M", left_name, left_orientation, fw_last, k - 1).unwrap();
                }
            }
        });
        gfa_out_mutex_ref.lock().unwrap().flush().unwrap();
    });

    n_segments
}

pub fn break_to_colored_subunitigs(unitig_colex_ranks: &[usize], _unitig_string: &[u8], map: &ColexToColorSetMap, sbwt: &SbwtIndex<SubsetMatrix>) -> (Vec<usize>, Vec<Range<usize>>){
    if unitig_colex_ranks.len() == 0 {
        // Make this a special case to ensure that there is always at least
        // one run to avoid a special case at the end.
        return (vec![], vec![]);
    }
    let mut subunitig_color_set_ids: Vec<usize> = vec![];
    let mut subunitigs: Vec<Range<usize>> = vec![]; // Ranges of k-mers (= starts of k-mers)
    let mut current_run_set_id = usize::MAX; // Will be set at the start of the first iteration
    let mut current_run_end = unitig_colex_ranks.len(); 

    // Iterate from end to start, updating the color set when the current
    // node is marked.
    for (pos, &colex) in unitig_colex_ranks.iter().enumerate().rev() {
        if pos == unitig_colex_ranks.len()-1 {
            assert!(map.sampling.get(colex)); // Last position of a unitig should always be marked
        }

        if map.sampling.get(colex) {
            // Update the set id
            let new_set_id = map.colex_to_color_set_id(colex, sbwt);

            if new_set_id != current_run_set_id {
                // Close the active run (if exists)
                let start = pos + 1;
                if current_run_end > start { // Active run exists
                    subunitigs.push(start..current_run_end);
                    subunitig_color_set_ids.push(current_run_set_id);
                    current_run_end = pos + 1;
                }
            }
            current_run_set_id = new_set_id;
        }
    }

    // Close the active run (exists because of the assert at the start)
    assert!(current_run_set_id != usize::MAX);
    subunitigs.push(0..current_run_end);
    subunitig_color_set_ids.push(current_run_set_id);

    subunitigs.reverse();
    subunitig_color_set_ids.reverse();

    (subunitig_color_set_ids, subunitigs)
}

#[allow(clippy::type_complexity)] // Yeah yeah I know
// Todo: add a function to dbg to get the k from there
fn search_unitig_from(v: Node, dbg: &Dbg<'_, SubsetMatrix>, k: usize, map: &ColexToColorSetMap, sbwt: &SbwtIndex<SubsetMatrix>) -> (Vec<usize>, Vec<usize>, Vec<u8>, Vec<Range<usize>>, Vec<usize>) {
    // Walk the unitig in forward orientation, and then backwards
    let mut workspace = Vec::<u8>::new();
    let nodes = dbg.walk_unitig_from(v, &mut workspace);
    workspace.clear();
    let mut unitig_string = Vec::<u8>::new();
    dbg.push_unitig_string(&nodes, &mut unitig_string);

    let string_len = unitig_string.len();
    assert!(string_len >= k);
    let last_kmer = &unitig_string[string_len-k..];
    let last_kmer_rc = reverse_complement(last_kmer);
    let last_kmer_rc_colex = dbg.get_node(&last_kmer_rc).unwrap_or_else(|| panic!(
        "Reverse complement of k-mer {} not found in index", 
        String::from_utf8_lossy(last_kmer))
    ).id;
    let rc_nodes = dbg.walk_unitig_from(sbwt::dbg::Node{id: last_kmer_rc_colex}, &mut workspace);

    let fw_colex: Vec<usize> = nodes.into_iter().map(|v| v.id).collect();
    let rc_colex: Vec<usize> = rc_nodes.into_iter().rev().map(|v| v.id).collect();
    assert_eq!(fw_colex.len(), rc_colex.len());

    // Figure out color set id runs in the forward strand 
    let (subuniting_color_set_ids, subunitig_kmer_ranges) = break_to_colored_subunitigs(&fw_colex, &unitig_string, map, sbwt);

    (fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subuniting_color_set_ids)
}

#[allow(clippy::too_many_arguments)] // Yeah yeah I know
fn visit_and_output_kmers(unitig_string: &[u8], subunitig_kmer_ranges: &[Range<usize>], subunitig_color_set_ids: &[usize], fw_colex: &[usize], rc_colex: &[usize], visited: &mut bitvec::vec::BitVec, unitigs_out: &mut impl Write, unitig_id: &mut usize, gfa_bvs: &mut Option<GfaBitVecs>, k: usize) {

    for (subunitig_idx, r) in subunitig_kmer_ranges.iter().enumerate() {
        // All k-mers in this subunitig have the same color set id.
        // It would be nice if we could just figure out the unvisited
        // runs of k-mers and visit and output those, but there is a subtle problem:
        // A subunitig may loop back to itself in reverse complement orientation.
        // Printing the subunitig would print the same k-mer in both orientations.
        // So, we need to keep track of the visited bit vector also while processing
        // a subunitig, and end the subunitig when we encounter a visited k-mer.
        let subunitig = &unitig_string[r.start..r.end+k-1];
        let color_set_id = subunitig_color_set_ids[subunitig_idx];
        let fw_colex_slice = &fw_colex[r.start..r.end];
        let rc_colex_slice = &rc_colex[r.start..r.end];

        let mut subsubunitig_start: Option<usize> = None;
        for kmer_idx in 0..fw_colex_slice.len() {
            if !visited[fw_colex_slice[kmer_idx]] {
                // Extend the current subunitig and visit this k-mer
                if subsubunitig_start.is_none() {
                    subsubunitig_start = Some(kmer_idx);
                }
                visited.set(fw_colex_slice[kmer_idx], true);
                visited.set(rc_colex_slice[kmer_idx], true);
            } else {
                // Already visited! Output the current subunitig
                if let Some(s) = subsubunitig_start {
                    let e = kmer_idx + k - 1;
                    if let Some(bvs) = gfa_bvs.as_mut() {
                        let fw_last = fw_colex_slice[kmer_idx - 1];
                        write!(unitigs_out, "S\tu{}\t", fw_last).unwrap();
                        unitigs_out.write_all(&subunitig[s..e]).unwrap();
                        writeln!(unitigs_out, "\tcs:i:{}", color_set_id).unwrap();
                        bvs.is_fw_first.set(fw_colex_slice[s], true);
                        bvs.is_fw_last.set(fw_colex_slice[kmer_idx - 1], true);
                        bvs.is_rc_last.set(rc_colex_slice[s], true);
                        bvs.is_rc_first.set(rc_colex_slice[kmer_idx - 1], true);
                    } else {
                        writeln!(unitigs_out, "> unitig_id={} color_set_id={}", unitig_id, color_set_id).unwrap();
                        unitigs_out.write_all(&subunitig[s..e]).unwrap();
                        unitigs_out.write_all(b"\n").unwrap();
                        *unitig_id += 1;
                    }
                }
                subsubunitig_start = None;
            }
        }

        // Write the last subunitig if it's still open
        if let Some(s) = subsubunitig_start {
            let e = fw_colex_slice.len() + k - 1;
            let last = fw_colex_slice.len() - 1;
            if let Some(bvs) = gfa_bvs.as_mut() {
                let fw_last = fw_colex_slice[last];
                write!(unitigs_out, "S\tu{}\t", fw_last).unwrap();
                unitigs_out.write_all(&subunitig[s..e]).unwrap();
                writeln!(unitigs_out, "\tcs:i:{}", color_set_id).unwrap();
                bvs.is_fw_first.set(fw_colex_slice[s], true);
                bvs.is_fw_last.set(fw_colex_slice[last], true);
                bvs.is_rc_last.set(rc_colex_slice[s], true);
                bvs.is_rc_first.set(rc_colex_slice[last], true);
            } else {
                writeln!(unitigs_out, "> unitig_id={} color_set_id={}", unitig_id, color_set_id).unwrap();
                unitigs_out.write_all(&subunitig[s..e]).unwrap();
                unitigs_out.write_all(b"\n").unwrap();
                *unitig_id += 1;
            }
        }
    }
}

/// Canonical here means whichever strand is visited first.
/// This assumes that the color set of a forward k-mer and a reverse k-mer is the same.
/// Returns the number of unitigs written
fn export_canonical_unitigs_with_shared_color_set(dbg_ref: &Dbg<SubsetMatrix>, map: &ColexToColorSetMap, sbwt: &SbwtIndex<SubsetMatrix>, k: usize, sbwt_len: usize, mut unitigs_out: impl Write + Sync + Send, n_threads: usize) -> usize {

    log::info!("Computing unitigs");
    let n_unitig_searches = std::sync::atomic::AtomicUsize::new(0);
    let n_unitig_searches_ref = &n_unitig_searches;

    let bar = indicatif::ProgressBar::new(sbwt_len as u64);
    let n_unitigs = std::thread::scope(|scope| {

        // Channels of tuples of with these fields: 
        //   * forward colex ranks 
        //   * reverse complement colex ranks
        //   * unitig string
        //   * colored subunitig k-mer ranges
        //   * color set ids of the colored subunitig ranges
        // TODO: less heap allocation
        let (worker_out, collector_in) = bounded::<(Vec<usize>, Vec<usize>, Vec<u8>, Vec<Range<usize>>, Vec<usize>)>(n_threads);

        // Create unitig search threads 
        let mut worker_handles = Vec::<_>::new();
        let bar_ref = &bar;
        for thread_id in 0..n_threads { 
            let worker_out_clone = worker_out.clone();
            let handle = scope.spawn(move || {
                // Iterating all colex positions that have remainder thread_id modulo number of threads
                let mut colex = thread_id;
                while colex < sbwt_len {
                    let v = Node { id: colex };
                    if !dbg_ref.is_dummy_colex_position(colex) && dbg_ref.is_first_kmer_of_unitig(v) {
                        n_unitig_searches_ref.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        worker_out_clone.send(search_unitig_from(v, dbg_ref, k, map, sbwt)).unwrap();
                    }
                    colex += n_threads;
                    if ((colex - thread_id)/n_threads) % 10000 == 0 {
                        bar_ref.inc(10000);
                        //eprintln!("number of unitig searches: {}", n_unitig_searches_ref.load(std::sync::atomic::Ordering::Relaxed));
                    }
                }
                log::info!("Thread {} finished", thread_id);
            });
            worker_handles.push(handle);
        }

        let collector_handle = scope.spawn(move || {
            // We maintain the visited bit vector so that when we mark a k-mer, we also mark its
            // reverse complement.
            let mut unitig_id = 0_usize;

            // Bitvector marking visited colex ranks 
            let mut visited = bitvec::bitvec![usize, Lsb0; 0; sbwt_len];

            // Process all non-cyclic unitigs
            while let Ok((fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subunitig_color_set_ids)) = collector_in.recv() {
                visit_and_output_kmers(&unitig_string, &subunitig_kmer_ranges, &subunitig_color_set_ids, &fw_colex, &rc_colex, &mut visited, &mut unitigs_out, &mut unitig_id, &mut None, k);
            }

            // Process remaining cyclic unitigs
            log::info!("Processing remaining cyclic unitigs");
            let n_acyclic = unitig_id; // This many unitigs have been written so far
            let mut colex = 0_usize;
            while colex < visited.len() {
                colex = match visited[colex..].first_zero() {
                    Some(i) => colex + i,
                    None => break,
                };
                if !dbg_ref.is_dummy_colex_position(colex) {
                    let (fw_colex, rc_colex, unitig_string, subunitig_kmer_ranges, subunitig_color_set_ids)
                    = search_unitig_from(Node { id: colex }, dbg_ref, k, map, sbwt);

                    // Make sure it's really cyclic
                    assert!(unitig_string[..k-1] == unitig_string[unitig_string.len()-(k-1)..]);

                    visit_and_output_kmers(&unitig_string, &subunitig_kmer_ranges, &subunitig_color_set_ids, &fw_colex, &rc_colex, &mut visited, &mut unitigs_out, &mut unitig_id, &mut None, k);
                }
                colex += 1;
            }
            unitigs_out.flush().unwrap();
            log::info!("Found {} cyclic unitigs", unitig_id - n_acyclic);
            unitig_id
        });

        for h in worker_handles { // Wait for the workers to finish
            h.join().unwrap();
        }

        drop(worker_out);

        // Wait for the collector to finish
        let n_unitigs = collector_handle.join().unwrap();

        #[allow(clippy::let_and_return)] // It's renaming of the variable. Clearer this way.
        n_unitigs
    });
    bar.finish();

    log::info!("Wrote {} unitigs", n_unitigs);
    n_unitigs
}


pub fn write_color_sets<CSS: ColorSetStorage + Sync + Send>(mut colors_out: impl Write, sets: &CSS) {
    // Each line should look something like this:
    // color_set_id=9 size=7 3 4 9 12 14 15 16
    for set_id in 0..sets.n_sets() {
        let set_view = sets.get_set_view(set_id);
        write!(colors_out, "color_set_id={} size={}", set_id, set_view.len()).unwrap();
        for color in set_view.iter() {
            write!(colors_out, " {}", color).unwrap(); // TODO: faster IO
        }
        writeln!(colors_out).unwrap();
    }

}

pub fn write_metadata<CSS: ColorSetStorage + Sync + Send>(mut metadata_out: impl Write, n_unitigs: Option<usize>, sets: &CSS, k: usize) {

    // Should look something like this
    // num_colors=3682
    // num_unitigs=9314735
    // num_color_sets=5591009
    // k=31
    metadata_out.write_all(format!("num_colors={}\n", sets.n_colors()).as_bytes()).unwrap();
    if let Some(n_unitigs) = n_unitigs {
        metadata_out.write_all(format!("num_unitigs={}\n", n_unitigs).as_bytes()).unwrap();
    } else {
        log::info!("Number of unitigs missing -> not written to metadata");
    }
    metadata_out.write_all(format!("num_color_sets={}\n", sets.n_sets()).as_bytes()).unwrap();
    metadata_out.write_all(format!("k={}\n", k).as_bytes()).unwrap();
}