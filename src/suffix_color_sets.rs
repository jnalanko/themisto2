// Counts the number of distinct color sets of k'-mers for k' = k_min..=k, where
// the color set of a k'-mer is the union of the color sets of the k-mers that
// have that k'-mer as a suffix. For each k', we also compute the size of a
// sparse-dense color set storage holding the distinct sets.
//
// The k-mers sharing a k'-suffix form a contiguous colex interval: a maximal run
// where lcs[i] >= k'. These intervals nest like LCP intervals of a suffix array,
// so we compute the unions for all k' with a single bottom-up stack traversal
// of the LCS array. An interval at depth d whose parent interval is at depth p
// is the interval of a k'-mer for every k' in (p, d], with the same union.
//
// k'-mers whose union is empty (only dummy k-mers have it as a suffix, or the
// k'-mer contains $) are not counted.
//
// The computation runs in three phases:
//
// A. Set ids. We walk the unitigs of the de Bruijn graph and carry the color
//    set id backward from the sampled k-mers, which gives the id of every k-mer
//    with about one graph step per k-mer. The (colex, id) pairs are written to
//    one bucket file per colex range. After this, the SBWT and the colex-to-set
//    map are dropped.
//
// B. Traversal. Each colex range is loaded from its bucket file and traversed
//    with the LCS stack. Each interval updates a single map entry
//    fingerprint -> (bitmask of k' values, set size). The map is sharded by the
//    fingerprint. If it would exceed the memory budget, shards are spilled to
//    disk: their entries and all later updates are appended to a file instead.
//
// C. Aggregation. Each shard (from memory, or rebuilt from its spill file) adds
//    each of its sets to the counters of the k' values in its bitmask.

use std::collections::hash_map::RandomState;
use std::hash::BuildHasher;
use std::io::Write;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use rustc_hash::FxHashMap;
use sbwt::LcsArray;
use simple_sds_sbwt::ops::BitVec;

use crate::colex_colored_kmers::CompactColexKmers;
use crate::coloring_interface::{ColorSetStorage, ColorSetView};
use crate::sparse_dense_storage::{is_dense_formula, SparseDenseStorage, StorageSizeEstimate};

/// Maximum number of k' values in one run, because the k' values of a set are stored in a u32 bitmask.
pub const MAX_N_K_VALUES: usize = 32;

const NO_SET: u32 = u32::MAX; // Set id of positions that have no entry in the bucket files (dummy k-mers)
const ID_BUF_LEN: usize = 2048; // Per thread and bucket, in phase A
const UPDATE_BUF_LEN: usize = 1024; // Per thread and shard, in phase B
const SPILL_BUF_BYTES: usize = 1 << 18; // Per spilled shard
const DEDUP_CACHE_LEN: usize = 1 << 14; // Per thread, in phase B
const PROGRESS_STEP: usize = 1 << 16;

pub struct Config {
    pub k_min: usize,
    pub n_threads: usize,
    pub temp_dir: PathBuf, // Bucket and spill files go into a fresh subdirectory of this
    pub mem_budget_bytes: usize, // Target peak memory of phases B and C. usize::MAX for no limit.
    pub n_shards: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SuffixColorSetStats {
    pub k_prime: usize,
    pub n_distinct_color_sets: usize,
    pub n_sparse: usize,
    pub n_dense: usize,
    pub total_sparse_elements: usize,
    pub storage_size: StorageSizeEstimate, // Of a sparse-dense storage of the distinct sets
}

impl SuffixColorSetStats {
    // Aggregates the stats from the sizes of the distinct sets
    #[cfg(test)]
    fn from_set_sizes(k_prime: usize, n_colors: usize, sizes: impl Iterator<Item = usize>) -> Self {
        let mut acc = KAccumulator::default();
        for size in sizes {
            acc.add(size, n_colors);
        }
        acc.into_stats(k_prime, n_colors)
    }
}

/// Information about the run, for logging and tests.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)] // Not all fields are read outside tests
pub struct RunInfo {
    pub n_ranges: usize,
    pub n_updates: usize, // Map updates after per-thread deduplication
    pub n_distinct_sets_all_k: usize, // Distinct sets over all k' values together
    pub n_spilled_shards: usize,
    pub spilled_bytes: usize,
    pub peak_map_bytes: usize,
}

// Fingerprint of a color set: a 64-bit SipHash with a random key.
struct Fingerprinter {
    h: RandomState,
}

impl Fingerprinter {
    fn fingerprint(&self, set: &[u32]) -> u64 {
        self.h.hash_one(set)
    }
}

// a := a ∪ b. Both must be sorted. Uses buf as scratch space.
fn union_into(a: &mut Vec<u32>, b: &[u32], buf: &mut Vec<u32>) {
    if b.is_empty() || a.as_slice() == b {
        return;
    }
    if a.is_empty() {
        a.extend_from_slice(b);
        return;
    }
    buf.clear();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => { buf.push(a[i]); i += 1; },
            std::cmp::Ordering::Greater => { buf.push(b[j]); j += 1; },
            std::cmp::Ordering::Equal => { buf.push(a[i]); i += 1; j += 1; },
        }
    }
    buf.extend_from_slice(&a[i..]);
    buf.extend_from_slice(&b[j..]);
    std::mem::swap(a, buf);
}

// Splits 0..n into roughly n_pieces ranges such that every split point i has lcs[i] < k_min.
fn split_range(lcs: &LcsArray, k_min: usize, n_pieces: usize) -> Vec<Range<usize>> {
    let n = lcs.len();
    let piece_len = n.div_ceil(n_pieces.max(1)).max(1);
    let mut ranges = vec![];
    let mut start = 0;
    while start < n {
        let mut end = std::cmp::min(start + piece_len, n);
        while end < n && lcs.access(end) >= k_min {
            end += 1;
        }
        ranges.push(start..end);
        start = end;
    }
    ranges
}

fn append_to_file(path: &Path, bytes: &[u8]) {
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)
        .unwrap_or_else(|e| panic!("Could not open {}: {}", path.display(), e));
    f.write_all(bytes).unwrap_or_else(|e| panic!("Could not write to {}: {}", path.display(), e));
}

// Index of the calling thread in the current rayon pool, with a spare slot for other threads
fn thread_slot(n_threads: usize) -> usize {
    rayon::current_thread_index().map_or(n_threads, |i| std::cmp::min(i, n_threads))
}

fn new_progress_bar(len: usize) -> indicatif::ProgressBar {
    let bar = indicatif::ProgressBar::new(len as u64);
    bar.enable_steady_tick(std::time::Duration::from_secs(1));
    bar
}

/*
 * Phase A: set ids into bucket files
 */

fn bucket_path(dir: &Path, bucket: usize) -> PathBuf {
    dir.join(format!("ids-{}.bin", bucket))
}

// Writes for each non-dummy k-mer the pair (offset in its range, set id) into the bucket file of its range.
fn write_set_id_buckets<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, ranges: &[Range<usize>], dir: &Path, n_threads: usize) {
    let starts: Vec<usize> = ranges.iter().map(|r| r.start).collect();
    let bucket_locks: Vec<Mutex<()>> = (0..ranges.len()).map(|_| Mutex::new(())).collect();

    // Buffers per thread and bucket. Entries are (offset << 32) | id.
    let thread_bufs: Vec<Mutex<Vec<Vec<u64>>>> = (0..=n_threads).map(|_| Mutex::new(vec![vec![]; ranges.len()])).collect();

    let flush = |bucket: usize, buf: &mut Vec<u64>| {
        let bytes: &[u8] = bytemuck::cast_slice(buf.as_slice());
        let _guard = bucket_locks[bucket].lock().unwrap();
        append_to_file(&bucket_path(dir, bucket), bytes);
        buf.clear();
    };

    log::info!("Initializing the de Bruijn graph");
    let dbg = sbwt::dbg::Dbg::new(index.sbwt(), Some(index.lcs()), n_threads);

    log::info!("Computing set ids along unitigs");
    let bar = new_progress_bar(index.sbwt().n_kmers());
    dbg.iter_unitigs_with_callback(|nodes| {
        let mut bufs = thread_bufs[thread_slot(n_threads)].lock().unwrap();
        let mut cur_id: Option<usize> = None;
        // The last k-mer of a unitig and every k-mer where the set changes are sampled,
        // so the id of a non-sampled k-mer is the id of the next sampled k-mer in the unitig.
        for v in nodes.iter().rev() {
            if cur_id.is_none() || index.get_map().sampling.get(v.id) {
                cur_id = Some(index.colex_to_set_id(v.id));
            }
            let id = cur_id.unwrap();
            assert!(id < NO_SET as usize, "Too many distinct color sets");
            let bucket = starts.partition_point(|&s| s <= v.id) - 1;
            let offset = v.id - starts[bucket];
            let buf = &mut bufs[bucket];
            buf.push(((offset as u64) << 32) | id as u64);
            if buf.len() >= ID_BUF_LEN {
                flush(bucket, buf);
            }
        }
        bar.inc(nodes.len() as u64);
    }, n_threads);
    bar.finish();

    for bufs in thread_bufs.iter() {
        let mut bufs = bufs.lock().unwrap();
        for (bucket, buf) in bufs.iter_mut().enumerate() {
            if !buf.is_empty() {
                flush(bucket, buf);
            }
        }
    }
}

// Reads the set ids of the range from its bucket file and deletes the file.
fn read_set_id_bucket(dir: &Path, bucket: usize, range_len: usize, ids: &mut Vec<u32>) {
    ids.clear();
    ids.resize(range_len, NO_SET);
    let path = bucket_path(dir, bucket);
    if !path.exists() {
        return; // Only dummies in this range
    }
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("Could not read {}: {}", path.display(), e));
    for entry in bytes.chunks_exact(8) {
        let x = u64::from_le_bytes(entry.try_into().unwrap());
        ids[(x >> 32) as usize] = x as u32;
    }
    std::fs::remove_file(&path).unwrap();
}

/*
 * Phase B: sharded map with spilling
 */

#[derive(Clone, Copy)]
struct Update {
    fp: u64,
    mask: u32, // Bit i means k' = k_min + i
    size: u32, // Set size, clamped to u32::MAX
}

impl Update {
    fn to_bytes(self) -> [u8; 16] {
        let mut b = [0u8; 16];
        b[0..8].copy_from_slice(&self.fp.to_le_bytes());
        b[8..12].copy_from_slice(&self.mask.to_le_bytes());
        b[12..16].copy_from_slice(&self.size.to_le_bytes());
        b
    }

    fn from_bytes(b: &[u8]) -> Self {
        Update {
            fp: u64::from_le_bytes(b[0..8].try_into().unwrap()),
            mask: u32::from_le_bytes(b[8..12].try_into().unwrap()),
            size: u32::from_le_bytes(b[12..16].try_into().unwrap()),
        }
    }
}

type ShardMap = FxHashMap<u64, (u32, u32)>; // fingerprint -> (mask, size)
type IndexedShards = Vec<(usize, Shard)>;

enum Shard {
    InMemory(ShardMap),
    Spilled(Vec<u8>), // Pending bytes not yet appended to the spill file
}

// Approximate heap size of a hash table with the given capacity: 16-byte slots plus one control byte each.
fn table_bytes(capacity: usize) -> usize {
    if capacity == 0 {
        return 0;
    }
    let buckets = if capacity < 8 { (capacity + 1).next_power_of_two() } else { (capacity * 8 / 7).next_power_of_two() };
    buckets * 17 + 16
}

struct ShardedMap {
    shards: Vec<Mutex<Shard>>,
    used_bytes: AtomicUsize,
    peak_bytes: AtomicUsize,
    budget_bytes: usize,
    n_spilled: AtomicUsize,
    spilled_bytes: AtomicUsize,
    n_updates: AtomicUsize,
    dir: PathBuf,
}

impl ShardedMap {
    fn new(n_shards: usize, budget_bytes: usize, dir: PathBuf) -> Self {
        Self {
            shards: (0..n_shards).map(|_| Mutex::new(Shard::InMemory(ShardMap::default()))).collect(),
            used_bytes: AtomicUsize::new(0),
            peak_bytes: AtomicUsize::new(0),
            budget_bytes,
            n_spilled: AtomicUsize::new(0),
            spilled_bytes: AtomicUsize::new(0),
            n_updates: AtomicUsize::new(0),
            dir,
        }
    }

    fn shard_of(&self, fp: u64) -> usize {
        // The fingerprint is uniformly random, so any bits will do. The hash table hashes the
        // whole fingerprint again, so using these bits here does not hurt the table.
        ((fp >> 32) as usize) % self.shards.len()
    }

    fn spill_path(&self, shard: usize) -> PathBuf {
        self.dir.join(format!("spill-{}.bin", shard))
    }

    fn flush_spill_buf(&self, shard: usize, buf: &mut Vec<u8>) {
        if !buf.is_empty() {
            append_to_file(&self.spill_path(shard), buf);
            self.spilled_bytes.fetch_add(buf.len(), Ordering::Relaxed);
            buf.clear();
        }
    }

    // Applies the updates, which must all belong to the given shard.
    fn apply(&self, shard_idx: usize, updates: &[Update]) {
        self.n_updates.fetch_add(updates.len(), Ordering::Relaxed);
        let mut shard = self.shards[shard_idx].lock().unwrap();

        if let Shard::InMemory(map) = &mut *shard {
            if map.len() + updates.len() > map.capacity() {
                // The table may grow. Reserve the memory of the new table up front, since the old
                // and the new table coexist while rehashing.
                let old_bytes = table_bytes(map.capacity());
                let new_bytes = table_bytes(map.len() + updates.len());
                let prev = self.used_bytes.fetch_add(new_bytes, Ordering::SeqCst);
                if prev.saturating_add(new_bytes) > self.budget_bytes {
                    self.used_bytes.fetch_sub(new_bytes, Ordering::SeqCst);
                    // Spill: move the entries of this shard to disk and free the table
                    let mut buf = Vec::<u8>::with_capacity(map.len() * 16);
                    for (&fp, &(mask, size)) in map.iter() {
                        buf.extend_from_slice(&Update { fp, mask, size }.to_bytes());
                    }
                    self.flush_spill_buf(shard_idx, &mut buf);
                    self.used_bytes.fetch_sub(old_bytes, Ordering::SeqCst);
                    *shard = Shard::Spilled(Vec::new());
                    self.n_spilled.fetch_add(1, Ordering::Relaxed);
                } else {
                    self.peak_bytes.fetch_max(prev + new_bytes, Ordering::Relaxed); // Old and new table both alive
                    map.reserve(updates.len());
                    // Replace the reservation and the old table by the actual new table
                    self.used_bytes.fetch_add(table_bytes(map.capacity()), Ordering::SeqCst);
                    self.used_bytes.fetch_sub(new_bytes + old_bytes, Ordering::SeqCst);
                }
            }
        }

        match &mut *shard {
            Shard::InMemory(map) => {
                for u in updates {
                    map.entry(u.fp).or_insert((0, u.size)).0 |= u.mask;
                }
            },
            Shard::Spilled(buf) => {
                for u in updates {
                    buf.extend_from_slice(&u.to_bytes());
                }
                if buf.len() >= SPILL_BUF_BYTES {
                    self.flush_spill_buf(shard_idx, buf);
                }
            },
        }
    }
}

// Per-thread state of phase B
struct Worker {
    update_bufs: Vec<Vec<Update>>, // Per shard
    dedup_cache: Vec<(u64, u32)>, // Direct-mapped (fingerprint, mask) cache of recent updates
    ids: Vec<u32>, // Set ids of the current range
}

impl Worker {
    fn new(n_shards: usize) -> Self {
        Self { update_bufs: vec![vec![]; n_shards], dedup_cache: vec![(0, 0); DEDUP_CACHE_LEN], ids: vec![] }
    }

    fn push(&mut self, map: &ShardedMap, u: Update) {
        // Skip the update if a recent update of the same set already covered these k' values.
        // Updates are idempotent and commutative (OR of masks), so skipping is safe.
        let slot = &mut self.dedup_cache[(u.fp as usize) & (DEDUP_CACHE_LEN - 1)];
        if slot.0 == u.fp {
            if slot.1 | u.mask == slot.1 {
                return;
            }
            slot.1 |= u.mask;
        } else {
            *slot = (u.fp, u.mask);
        }

        let shard = map.shard_of(u.fp);
        let buf = &mut self.update_bufs[shard];
        buf.push(u);
        if buf.len() >= UPDATE_BUF_LEN {
            map.apply(shard, buf);
            buf.clear();
        }
    }

    fn flush(&mut self, map: &ShardedMap) {
        for (shard, buf) in self.update_bufs.iter_mut().enumerate() {
            if !buf.is_empty() {
                map.apply(shard, buf);
                buf.clear();
            }
        }
    }
}

struct Traversal<'a, CSS: ColorSetStorage> {
    lcs: &'a LcsArray,
    sets: &'a CSS,
    k: usize,
    k_min: usize,
    fingerprinter: Fingerprinter,
    map: &'a ShardedMap,
}

impl<CSS: ColorSetStorage> Traversal<'_, CSS> {

    // Records the union of an interval at depth d with parent depth p, i.e.
    // the union of the k'-mers for k' in (p, d].
    fn report(&self, worker: &mut Worker, union: &[u32], p: usize, d: usize) {
        let first = std::cmp::max(p + 1, self.k_min);
        if union.is_empty() || first > d {
            return;
        }
        // Bits first - k_min ..= d - k_min. At most 32 bits, so u64 arithmetic cannot overflow.
        let mask = (((1_u64 << (d - first + 1)) - 1) << (first - self.k_min)) as u32;
        let size = std::cmp::min(union.len(), u32::MAX as usize) as u32; // Clamps only the full set of 2^32 colors
        worker.push(self.map, Update { fp: self.fingerprinter.fingerprint(union), mask, size });
    }

    // Processes the colex range using the set ids in worker.ids. The caller must guarantee that
    // no interval of a k'-mer with k' >= k_min crosses the boundaries of the range, i.e.
    // lcs[range.start] < k_min and lcs[range.end] < k_min (where defined).
    fn process(&self, range: Range<usize>, worker: &mut Worker, bar: &indicatif::ProgressBar) {
        let ids = std::mem::take(&mut worker.ids);

        // Stack of (depth, union). The bottom is the root at depth 0, which is never reported.
        let mut stack: Vec<(usize, Vec<u32>)> = vec![(0, vec![])];
        let mut buf = Vec::<u32>::new();
        let mut prev_leaf: Option<(u32, Vec<u32>)> = None; // (set id, set) of the previous non-dummy k-mer

        for i in range.start..=range.end {
            // Treat the boundaries of the range as LCS 0
            let l = if i == range.start || i == range.end { 0 } else { self.lcs.access(i) };

            // Close the intervals deeper than l
            while stack.last().unwrap().0 > l {
                let (d, u) = stack.pop().unwrap();
                let top = stack.last_mut().unwrap();
                let p = std::cmp::max(l, top.0);
                self.report(worker, &u, p, d);
                if top.0 >= l {
                    if top.1.is_empty() {
                        top.1 = u;
                    } else {
                        union_into(&mut top.1, &u, &mut buf);
                    }
                } else {
                    // New internal node at depth l that starts with the popped child
                    stack.push((l, u));
                }
            }

            if i == range.end {
                break;
            }

            // Push the k-mer at colex i as a leaf at depth k
            let id = ids[i - range.start];
            let leaf = if id == NO_SET {
                vec![] // Dummy k-mer
            } else {
                match &prev_leaf {
                    Some((prev_id, prev_set)) if *prev_id == id => prev_set.clone(),
                    _ => {
                        let set: Vec<u32> = self.sets.get_set_view(id as usize).iter().map(|c| c as u32).collect();
                        prev_leaf = Some((id, set.clone()));
                        set
                    }
                }
            };
            stack.push((self.k, leaf));

            if (i - range.start + 1) % PROGRESS_STEP == 0 {
                bar.inc(PROGRESS_STEP as u64);
            }
        }
        bar.inc((range.len() % PROGRESS_STEP) as u64);

        worker.ids = ids;
    }
}

/*
 * Phase C: aggregation
 */

#[derive(Debug, Clone, Default)]
struct KAccumulator {
    n_sparse: usize,
    n_dense: usize,
    total_sparse_elements: usize,
}

impl KAccumulator {
    fn add(&mut self, size: usize, n_colors: usize) {
        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;
        if is_dense_formula(size, color_id_bit_width, n_colors) {
            self.n_dense += 1;
        } else {
            self.n_sparse += 1;
            self.total_sparse_elements += size;
        }
    }

    fn merge(&mut self, other: &Self) {
        self.n_sparse += other.n_sparse;
        self.n_dense += other.n_dense;
        self.total_sparse_elements += other.total_sparse_elements;
    }

    fn into_stats(self, k_prime: usize, n_colors: usize) -> SuffixColorSetStats {
        let storage_size = SparseDenseStorage::serialized_size_estimate(n_colors, self.n_sparse, self.n_dense, self.total_sparse_elements);
        SuffixColorSetStats {
            k_prime,
            n_distinct_color_sets: self.n_sparse + self.n_dense,
            n_sparse: self.n_sparse,
            n_dense: self.n_dense,
            total_sparse_elements: self.total_sparse_elements,
            storage_size,
        }
    }
}

// Adds the sets of the map to the accumulators. Returns the number of sets.
fn aggregate_map(map: &ShardMap, accs: &mut [KAccumulator], n_colors: usize) -> usize {
    for &(mask, size) in map.values() {
        let mut bits = mask;
        while bits != 0 {
            accs[bits.trailing_zeros() as usize].add(size as usize, n_colors);
            bits &= bits - 1;
        }
    }
    map.len()
}

fn read_spill_file(path: &Path) -> ShardMap {
    let mut map = ShardMap::default();
    if path.exists() {
        let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("Could not read {}: {}", path.display(), e));
        for entry in bytes.chunks_exact(16) {
            let u = Update::from_bytes(entry);
            map.entry(u.fp).or_insert((0, u.size)).0 |= u.mask;
        }
        std::fs::remove_file(path).unwrap();
    }
    map
}

/*
 * Driver
 */

struct CountingWriter(usize);

impl Write for CountingWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0 += buf.len();
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// Memory of the per-thread buffers of phase B and the pending bytes of spilled shards
fn buffer_bytes(n_threads: usize, n_shards: usize, max_range_len: usize) -> usize {
    (n_threads + 1) * (n_shards * UPDATE_BUF_LEN * 16 + DEDUP_CACHE_LEN * 16 + max_range_len * 4) + n_shards * SPILL_BUF_BYTES
}

type Accumulators = (Vec<KAccumulator>, usize); // Per k', and the number of distinct sets over all k'

fn combine(mut a: Accumulators, b: Accumulators) -> Accumulators {
    for (x, y) in a.0.iter_mut().zip(b.0.iter()) {
        x.merge(y);
    }
    (a.0, a.1 + b.1)
}

/// Returns the stats of the distinct non-empty color sets of k'-mers for k' = k_min..=k. Consumes
/// the index so that the SBWT and the colex-to-set map can be freed after phase A.
/// The SBWT must have select support (required by the de Bruijn graph).
pub fn count_distinct_suffix_color_sets<CSS: ColorSetStorage + Sync + Send>(index: CompactColexKmers<CSS>, config: &Config) -> (Vec<SuffixColorSetStats>, RunInfo) {
    let k = index.get_k();
    let k_min = config.k_min;
    let n_threads = config.n_threads;
    assert!(k_min >= 1 && k_min <= k, "k_min must be in the range [1, k]");
    let n_k_values = k - k_min + 1;
    assert!(n_k_values <= MAX_N_K_VALUES, "At most {} values of k' are supported per run", MAX_N_K_VALUES);
    let n_colors = index.get_set_storage().n_colors();

    let dir = config.temp_dir.join(format!("suffix-color-sets-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("Could not create {}: {}", dir.display(), e));

    let ranges = split_range(index.lcs(), k_min, n_threads * 64);
    let max_range_len = ranges.iter().map(|r| r.len()).max().unwrap_or(0);
    assert!(max_range_len < u32::MAX as usize, "A colex range is too long ({} positions). Try a larger k_min.", max_range_len);
    log::info!("Split the colex space into {} ranges (longest: {} positions)", ranges.len(), max_range_len);

    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    // Phase A
    let start = Instant::now();
    log::info!("Phase A: computing set ids into bucket files in {}", dir.display());
    write_set_id_buckets(&index, &ranges, &dir, n_threads);
    log::info!("Phase A took {:.1} s", start.elapsed().as_secs_f64());

    let (sbwt, lcs, colex_map, sets, _names) = index.into_parts();
    drop(sbwt);
    drop(colex_map);

    // Memory budget of the map: the total budget minus the LCS, the color sets, and buffers
    let lcs_bytes = lcs.size_in_bytes();
    let sets_bytes = { let mut cw = CountingWriter(0); sets.serialize(&mut cw); cw.0 };
    let buffer_bytes = buffer_bytes(n_threads, config.n_shards, max_range_len);
    let map_budget = config.mem_budget_bytes.saturating_sub(lcs_bytes + sets_bytes + buffer_bytes);
    if config.mem_budget_bytes != usize::MAX {
        log::info!("Memory: LCS {:.2} GB, color sets {:.2} GB, buffers up to {:.2} GB, leaving {:.2} GB for the map",
            lcs_bytes as f64 / 1e9, sets_bytes as f64 / 1e9, buffer_bytes as f64 / 1e9, map_budget as f64 / 1e9);
    }

    // Phase B
    let start = Instant::now();
    log::info!("Phase B: traversing {} colex ranges", ranges.len());
    let map = ShardedMap::new(config.n_shards, map_budget, dir.clone());
    let traversal = Traversal { lcs: &lcs, sets: &sets, k, k_min, fingerprinter: Fingerprinter { h: RandomState::new() }, map: &map };
    let workers: Vec<Mutex<Worker>> = (0..=n_threads).map(|_| Mutex::new(Worker::new(config.n_shards))).collect();
    let bar = new_progress_bar(lcs.len());
    pool.install(|| {
        ranges.par_iter().enumerate().for_each(|(bucket, range)| {
            let mut worker = workers[thread_slot(n_threads)].lock().unwrap();
            read_set_id_bucket(&dir, bucket, range.len(), &mut worker.ids);
            traversal.process(range.clone(), &mut worker, &bar);
        });
    });
    bar.finish();
    for worker in workers.iter() {
        worker.lock().unwrap().flush(&map);
    }
    drop(workers);
    for (shard_idx, shard) in map.shards.iter().enumerate() {
        if let Shard::Spilled(buf) = &mut *shard.lock().unwrap() {
            map.flush_spill_buf(shard_idx, buf);
        }
    }
    let info = RunInfo {
        n_ranges: ranges.len(),
        n_updates: map.n_updates.load(Ordering::Relaxed),
        n_distinct_sets_all_k: 0, // Filled in after phase C
        n_spilled_shards: map.n_spilled.load(Ordering::Relaxed),
        spilled_bytes: map.spilled_bytes.load(Ordering::Relaxed),
        peak_map_bytes: map.peak_bytes.load(Ordering::Relaxed),
    };
    log::info!("Phase B took {:.1} s. {} map updates. Map peak {:.2} GB. {} of {} shards spilled ({:.2} GB on disk)",
        start.elapsed().as_secs_f64(), info.n_updates, info.peak_map_bytes as f64 / 1e9,
        info.n_spilled_shards, config.n_shards, info.spilled_bytes as f64 / 1e9);
    drop(lcs);

    // Phase C. In-memory shards first, freeing each one, so that the spilled shards have room.
    let start = Instant::now();
    log::info!("Phase C: aggregating");
    let empty = || -> Accumulators { (vec![KAccumulator::default(); n_k_values], 0) };
    let shards: IndexedShards = map.shards.into_iter().map(|s| s.into_inner().unwrap()).enumerate().collect();
    let (in_memory, spilled): (IndexedShards, IndexedShards) = shards.into_iter().partition(|(_, s)| matches!(s, Shard::InMemory(_)));
    let (accs, n_distinct_all_k) = pool.install(|| {
        let from_memory = in_memory.into_par_iter().fold(empty, |mut acc, (_, shard)| {
            if let Shard::InMemory(m) = shard {
                acc.1 += aggregate_map(&m, &mut acc.0, n_colors);
            }
            acc
        }).reduce(empty, combine);
        let from_disk = spilled.into_par_iter().fold(empty, |mut acc, (shard_idx, _)| {
            let m = read_spill_file(&dir.join(format!("spill-{}.bin", shard_idx)));
            acc.1 += aggregate_map(&m, &mut acc.0, n_colors);
            acc
        }).reduce(empty, combine);
        combine(from_memory, from_disk)
    });
    log::info!("Phase C took {:.1} s. {} distinct sets over all k' values", start.elapsed().as_secs_f64(), n_distinct_all_k);

    std::fs::remove_dir_all(&dir).ok();

    let info = RunInfo { n_distinct_sets_all_k: n_distinct_all_k, ..info };
    let stats = accs.into_iter().enumerate().map(|(i, acc)| acc.into_stats(k_min + i, n_colors)).collect();
    (stats, info)
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeSet, HashMap, HashSet};

    use rand::{Rng, SeedableRng};

    use super::*;

    fn brute_force(index: &CompactColexKmers<SparseDenseStorage>, k_min: usize) -> Vec<SuffixColorSetStats> {
        let k = index.get_k();
        let mut kmers = vec![];
        for colex in 0..index.sbwt().n_sets() {
            let kmer = index.sbwt().access_kmer(colex);
            if kmer.contains(&b'$') {
                continue;
            }
            let set: BTreeSet<usize> = index.colex_to_set(colex).iter().collect();
            kmers.push((kmer, set));
        }

        (k_min..=k).map(|k_prime| {
            let mut unions = HashMap::<Vec<u8>, BTreeSet<usize>>::new();
            for (kmer, set) in kmers.iter() {
                unions.entry(kmer[k - k_prime..].to_vec()).or_default().extend(set.iter().copied());
            }
            let distinct: HashSet<&BTreeSet<usize>> = unions.values().collect();
            SuffixColorSetStats::from_set_sizes(k_prime, index.get_set_storage().n_colors(), distinct.iter().map(|s| s.len()))
        }).collect()
    }

    fn random_seq(rng: &mut impl Rng, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
    }

    // Mutates random positions to get related sequences with partially shared k-mers
    fn mutate(rng: &mut impl Rng, seq: &[u8], n_mutations: usize) -> Vec<u8> {
        let mut seq = seq.to_vec();
        for _ in 0..n_mutations {
            let i = rng.gen_range(0..seq.len());
            seq[i] = b"ACGT"[rng.gen_range(0..4)];
        }
        seq
    }

    fn test_dir() -> PathBuf {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        std::env::temp_dir().join(format!("themisto2-suffix-color-sets-test-{}-{}", std::process::id(), COUNTER.fetch_add(1, Ordering::Relaxed)))
    }

    fn run(index: CompactColexKmers<SparseDenseStorage>, k_min: usize, n_threads: usize, mem_budget_bytes: usize, n_shards: usize) -> (Vec<SuffixColorSetStats>, RunInfo) {
        let dir = test_dir();
        let config = Config { k_min, n_threads, temp_dir: dir.clone(), mem_budget_bytes, n_shards };
        let result = count_distinct_suffix_color_sets(index, &config);
        std::fs::remove_dir_all(&dir).ok();
        result
    }

    #[test]
    fn matches_brute_force() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1234);
        // With more than 64 colors, small sets are sparse and large sets are dense
        for &(k, n_colors, seq_len) in &[(5, 3, 50), (7, 6, 200), (11, 10, 300), (31, 8, 500), (9, 100, 150)] {
            let base = random_seq(&mut rng, seq_len);
            let seqs: Vec<Vec<u8>> = (0..n_colors).map(|_| mutate(&mut rng, &base, seq_len / 20)).collect();
            let colored_seqs: Vec<(&[u8], usize)> = seqs.iter().enumerate().map(|(c, s)| (s.as_slice(), c)).collect();
            let build = || CompactColexKmers::<SparseDenseStorage>::new_from_small_input(&colored_seqs, k, 3, 2);
            let reference = build();

            for k_min in [1, 2, k / 2, k] {
                if k - k_min + 1 > MAX_N_K_VALUES {
                    continue;
                }
                let expected = brute_force(&reference, k_min);
                // Unlimited memory, no memory (every shard spills), and a small budget (some shards spill)
                for (budget, n_threads) in [(usize::MAX, 1), (usize::MAX, 3), (0, 2), (2000, 3)] {
                    let (got, info) = run(build(), k_min, n_threads, budget, 8);
                    assert_eq!(got, expected, "k = {}, k_min = {}, n_threads = {}, budget = {}", k, k_min, n_threads, budget);
                    if budget == usize::MAX {
                        assert_eq!(info.n_spilled_shards, 0);
                    }
                    if budget == 0 && info.n_distinct_sets_all_k > 0 {
                        assert!(info.n_spilled_shards > 0);
                    }
                }
            }

            // At k' = k, the distinct sets are the non-empty sets of the storage of the index.
            // The storage built by new_from_small_input also contains the empty set of the dummy
            // k-mers, so we include it when checking the size estimate against the real storage.
            let storage = reference.get_set_storage();
            let sizes: Vec<usize> = (0..storage.n_sets()).map(|i| storage.get_set_view(i).len()).collect();
            let (got, _) = run(build(), k, 2, usize::MAX, 8);
            if n_colors > 64 {
                assert!(got[0].n_sparse > 0 && got[0].n_dense > 0, "Expected both sparse and dense sets: {:?}", got[0]);
            }
            assert_eq!(got[0], SuffixColorSetStats::from_set_sizes(k, storage.n_colors(), sizes.iter().copied().filter(|&s| s > 0)), "k = {}", k);

            let mut buf = Vec::<u8>::new();
            storage.serialize(&mut buf);
            let with_empty = SuffixColorSetStats::from_set_sizes(k, storage.n_colors(), sizes.iter().copied());
            assert_eq!(with_empty.storage_size.total(), buf.len(), "k = {}", k);
        }
    }

    #[test]
    fn partial_spill() {
        // A budget that fits some shards but not all, to exercise the mix of in-memory and spilled shards
        let mut rng = rand::rngs::StdRng::seed_from_u64(99);
        let base = random_seq(&mut rng, 3000);
        let seqs: Vec<Vec<u8>> = (0..20).map(|_| mutate(&mut rng, &base, 150)).collect();
        let colored_seqs: Vec<(&[u8], usize)> = seqs.iter().enumerate().map(|(c, s)| (s.as_slice(), c)).collect();
        let build = || CompactColexKmers::<SparseDenseStorage>::new_from_small_input(&colored_seqs, 15, 3, 2);

        let (expected, full) = run(build(), 6, 2, usize::MAX, 16);
        assert_eq!(full.n_spilled_shards, 0);
        assert_eq!(expected, brute_force(&build(), 6));

        // The full map needs about peak_map_bytes. Give the map roughly half of it, on top of
        // what count_distinct_suffix_color_sets reserves for the LCS, the sets and the buffers.
        let index = build();
        let mut sets_buf = Vec::<u8>::new();
        index.get_set_storage().serialize(&mut sets_buf);
        let max_range_len = split_range(index.lcs(), 6, 2 * 64).iter().map(|r| r.len()).max().unwrap();
        let reserved = index.lcs().size_in_bytes() + sets_buf.len() + buffer_bytes(2, 16, max_range_len);
        let (got, info) = run(index, 6, 2, reserved + full.peak_map_bytes / 2, 16);
        assert_eq!(got, expected);
        assert!(info.n_spilled_shards > 0 && info.n_spilled_shards < 16, "Expected some but not all shards to spill: {:?}", info);
    }
}
