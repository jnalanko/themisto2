// This is dead code for now, might re-integrate later.

//! Direct construction of a Themisto 2 index from a ggcat colored graph,
//! reading through the ggcat API instead of via an intermediate text dump.
//!
//! Architecture: one `dump_unitigs` pass on a driver thread, with two
//! consumers:
//!   * `GgcatColorSetGenerator` feeds new color subsets into
//!     `ColorSetStorage::new` on a CSS thread.
//!   * A worker pool consumes batched `(unitig, color_set_id)` pairs and
//!     accumulates `(colex, color_set_id)` pairs against the SBWT, mirroring
//!     `CompactColexKmers::new_from_colored_unitig_dump`.

use std::path::Path;

use crossbeam::channel::{bounded, Receiver};
use ggcat_api::{GGCATConfig, GGCATInstance};
use parking_lot::Mutex;
use rayon::slice::ParallelSliceMut;
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};
use simple_sds_sbwt::raw_vector::AccessRaw;

use crate::colex_colored_kmers::{
    ColexToColorSetMap, CompactColexKmers, UnitigImportSeqBatch,
};
use crate::coloring_interface::ColorSetStorage;
use crate::int_vec::CompactIntVec;
use crate::iterators::{USizeIterator, USizeIteratorGenerator};

const UNITIG_BATCH_BYTES: usize = 1 << 20; // 1 MiB; matches the on-disk-import producer

/// Borrowed iterator over a `Vec<usize>` exposing the `USizeIterator` trait.
pub struct VecIterRef<'a> {
    slice: &'a [usize],
    pos: usize,
}

impl<'a> VecIterRef<'a> {
    fn new(slice: &'a [usize]) -> Self {
        Self { slice, pos: 0 }
    }
}

impl<'a> USizeIterator<'a> for VecIterRef<'a> {
    fn next(&mut self) -> Option<usize> {
        if self.pos == self.slice.len() {
            None
        } else {
            let x = self.slice[self.pos];
            self.pos += 1;
            Some(x)
        }
    }
}

/// `USizeIteratorGenerator` backed by a channel of color sets. Each `recv`
/// yields the next distinct color subset (in ggcat output order).
pub struct GgcatColorSetGenerator {
    rx: Receiver<Vec<usize>>,
    buf: Vec<usize>,
}

impl GgcatColorSetGenerator {
    pub fn new(rx: Receiver<Vec<usize>>) -> Self {
        Self { rx, buf: vec![] }
    }
}

impl USizeIteratorGenerator for GgcatColorSetGenerator {
    type Iter<'a> = VecIterRef<'a>;

    fn next<'b>(&'b mut self) -> Option<Self::Iter<'b>> {
        match self.rx.recv() {
            Ok(v) => {
                self.buf = v;
                Some(VecIterRef::new(&self.buf))
            }
            Err(_) => None, // Sender dropped — stream complete
        }
    }
}

struct DriverState {
    color_set_counter: usize,
    last_color_set_id: usize,
    cur_concat: Vec<u8>,
    cur_starts: Vec<usize>,
    cur_color_set_ids: Vec<usize>,
}

impl DriverState {
    fn new() -> Self {
        Self {
            color_set_counter: 0,
            last_color_set_id: 0,
            cur_concat: Vec::with_capacity(UNITIG_BATCH_BYTES),
            cur_starts: Vec::new(),
            cur_color_set_ids: Vec::new(),
        }
    }

    fn take_batch(&mut self) -> UnitigImportSeqBatch {
        let mut concat = std::mem::replace(
            &mut self.cur_concat,
            Vec::with_capacity(UNITIG_BATCH_BYTES),
        );
        let mut starts = std::mem::take(&mut self.cur_starts);
        let color_set_ids = std::mem::take(&mut self.cur_color_set_ids);
        starts.push(concat.len()); // End sentinel, as required by UnitigImportSeqBatch
        concat.shrink_to_fit();
        starts.shrink_to_fit();
        UnitigImportSeqBatch { concat, starts, color_set_ids }
    }
}

/// Build a Themisto 2 index from a ggcat colored graph using a single
/// `dump_unitigs` pass.
///
/// `sbwt`/`lcs` must already cover the same unitigs (build them directly from
/// the unitig FASTA via the normal `get_sbwt_and_lcs` path).
pub fn build_index_from_ggcat<CSS: ColorSetStorage + Send + Sync + 'static>(
    sbwt: SbwtIndex<SubsetMatrix>,
    lcs: LcsArray,
    ggcat_unitigs: &Path,
    ggcat_colors: &Path,
    temp_dir: &Path,
    k: usize,
    sample_distance: usize,
    n_threads: usize,
) -> CompactColexKmers<CSS> {
    assert!(sample_distance > 0);

    log::info!("Initializing GGCAT instance");
    // ggcat runs as the single producer feeding our colex workers — analogous
    // to `unitig_import_parser_thread` in `new_from_colored_unitig_dump`. Cap
    // its thread pool at 1 so total concurrent worker threads stays at
    // `n_threads` (the colex pool) instead of `2 * n_threads`.
    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(temp_dir.to_path_buf()),
        memory: 2.0,
        prefer_memory: true,
        total_threads_count: 1,
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .expect("GGCATInstance::create failed");

    log::info!("Reading color names from {}", ggcat_colors.display());
    let color_names: Vec<String> = GGCATInstance::dump_colors(ggcat_colors)
        .expect("Failed to read ggcat colormap")
        .collect();
    let n_colors = color_names.len();
    log::info!("ggcat colormap has {} colors", n_colors);

    let (css_tx, css_rx) = bounded::<Vec<usize>>(64);
    let (batch_tx, batch_rx) = bounded::<UnitigImportSeqBatch>(n_threads * 2);

    let sbwt_ref = &sbwt;
    let lcs_ref = &lcs;

    let driver_state = Mutex::new(DriverState::new());

    let (css_boxed, colex_pairs) = std::thread::scope(|s| {
        // Move owning clones into the driver closure so the channels close
        // automatically when the driver thread exits.
        let driver_css_tx = css_tx.clone();
        let driver_batch_tx = batch_tx.clone();
        let driver_state_ref = &driver_state;
        let driver = s.spawn(move || {
            log::info!("Streaming unitigs from {}", ggcat_unitigs.display());
            instance
                .dump_unitigs(
                    ggcat_unitigs,
                    k,
                    None,
                    true, // colors required
                    1,    // ggcat is the producer; parallelism lives in our colex worker pool
                    true, // single_thread_output_function: callback is serialized for us
                    |seq, colors, same_colors| {
                        // `single_thread_output_function=true` serializes our
                        // callback; the mutex is only here for interior
                        // mutability of the captured driver state. Holding it
                        // across the channel sends just provides natural
                        // backpressure when the consumers fall behind.
                        let mut st = driver_state_ref.lock();
                        let id = if same_colors {
                            st.last_color_set_id
                        } else {
                            let new_id = st.color_set_counter;
                            st.color_set_counter += 1;
                            let v: Vec<usize> =
                                colors.iter().map(|c| *c as usize).collect();
                            driver_css_tx.send(v).expect("CSS channel closed");
                            new_id
                        };
                        st.last_color_set_id = id;

                        let cur_len = st.cur_concat.len();
                        st.cur_starts.push(cur_len);
                        st.cur_concat.extend_from_slice(seq);
                        st.cur_color_set_ids.push(id);

                        if st.cur_concat.len() >= UNITIG_BATCH_BYTES {
                            let batch = st.take_batch();
                            driver_batch_tx
                                .send(batch)
                                .expect("Unitig batch channel closed");
                        }
                    },
                )
                .expect("dump_unitigs failed");

            // Flush trailing partial batch, if any.
            let mut st = driver_state_ref.lock();
            if !st.cur_concat.is_empty() {
                let batch = st.take_batch();
                drop(st);
                driver_batch_tx
                    .send(batch)
                    .expect("Unitig batch channel closed");
            }
        });

        // Drop the local Sender originals so consumer threads see disconnect
        // as soon as the driver (the only remaining owner) finishes.
        drop(css_tx);
        drop(batch_tx);

        // CSS consumer thread.
        let css_handle = s.spawn(move || {
            let color_set_gen = GgcatColorSetGenerator::new(css_rx);
            CSS::new(color_set_gen, n_colors)
        });

        // Colex worker pool.
        let mut workers = Vec::with_capacity(n_threads);
        for _ in 0..n_threads {
            let rx = batch_rx.clone();
            workers.push(s.spawn(move || {
                let index = StreamingIndex::new(sbwt_ref, lcs_ref);
                let mut local: Vec<(usize, usize)> = Vec::new();
                while let Ok(batch) = rx.recv() {
                    batch.process(&mut local, &index, sample_distance);
                }
                local
            }));
        }
        // Drop our local Receiver clone so the workers see disconnect once the
        // driver closes its end.
        drop(batch_rx);

        driver.join().expect("driver thread panicked");

        let css_boxed = css_handle.join().expect("CSS thread panicked");

        let mut colex_pairs: Vec<(usize, usize)> = Vec::new();
        for w in workers {
            colex_pairs.extend(w.join().expect("colex worker panicked"));
        }
        (css_boxed, colex_pairs)
    });

    let css = *css_boxed;
    let n_color_sets = css.n_sets();
    log::info!(
        "Collected {} (colex, color_set_id) pairs across {} color sets",
        colex_pairs.len(),
        n_color_sets
    );

    log::info!("Sorting (colex, color set id) pairs");
    let mut colex_pairs = colex_pairs;
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap().install(|| {
        colex_pairs.par_sort_unstable();
    });

    let bit_width = n_color_sets.next_power_of_two().trailing_zeros() as usize;
    log::info!("Building compressed representation for color set ids");
    let mut stored_color_set_ids = CompactIntVec::new(colex_pairs.len(), bit_width);
    let mut sample_marks =
        simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
    for (rank, (colex, id)) in colex_pairs.into_iter().enumerate() {
        stored_color_set_ids.set(rank, id);
        sample_marks.set_bit(colex, true);
    }
    let mut sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(sample_marks);
    sample_marks.enable_rank();
    log::info!(
        "Marked {:.2} % of all k-mers",
        sample_marks.count_ones() as f64 / sbwt.n_kmers() as f64 * 100.0
    );

    let colex_map = ColexToColorSetMap {
        sampling: sample_marks,
        color_set_ids: stored_color_set_ids,
    };

    CompactColexKmers::<CSS>::new(sbwt, lcs, colex_map, css, Some(&color_names))
}
