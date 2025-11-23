use bitvec::order::Lsb0;
use bitvec::{field::BitField, slice::BitSlice};
use bitvec::{bitvec, vec};
use crossbeam::channel::{Sender, unbounded};
use std::io::BufRead;
use sbwt::{MergeInterleaving, reverse_complement_in_place};
use sbwt::LcsArray;
use sbwt::{dbg::Node, SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::serialize::Serialize;
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Push, Rank, Resize, Vector}, raw_vector::AccessRaw};
use rustc_hash::FxHasher;
use std::cmp::max;
use std::sync::Arc;
use std::{cmp::min, collections::HashMap, hash::BuildHasherDefault, sync::Mutex};
use std::hash::{Hash, Hasher};

use crate::int_vec::{AtomicCompactIntVec, CompactIntVec};
use crate::coloring_interface::{self, ColorSetOwned, ColorSetStorage, ColorSetView};
use crate::index_import::{self, parse_color_set_line};
use crate::iterators::VecVecUsizeIteratorGenerator;

/// This is the main data structure in this file: a set of compressed color sets, and a mapping
/// from SBWT colex ranks to color sets such that we can look up the color set of a k-mer by its
/// colex rank in the SBWT. 
pub struct CompactColexKmers<CSS: coloring_interface::ColorSetStorage> {
    // This is on the heap to allow map to refer to it (otherwise assuring lifetime 
    // guarantees becomes problematic). It's reference counted because this struct
    // will have two references to it, the one in self.sbwt, and one in self.map.sbwt.
    // Note that this means that if we replace sbwt here with a new Arc pointing to a new
    // sbwt, then, the map will continue to point to the old sbwt. So don't do that!
    // It's atomic (Arc) because we want to pass this struct to multiple threads.
    sbwt: Arc<SbwtIndex<SubsetMatrix>>, 

    lcs: LcsArray,
    sets: CSS, // Distinct color sets
    map: ColexToColorSetMap, // A mapping from the colex rank of a k-mer in the SBWT into a color set id in `sets`
}

/// A data structure that stores the color set ids for a subset of sampled k-mers in the SBWT such that
/// the color sets of the rest can be obtained by walking forward in the de Bruijn graph to the
/// closest sampled node.
pub struct ColexToColorSetMap {

    // See the comments inside CompactcolexColoring
    sbwt: Arc<SbwtIndex<SubsetMatrix>>,

    sampling: simple_sds_sbwt::bit_vector::BitVector, // Marks colex ranks that have a color set stored. Has rank support.
    color_set_ids: CompactIntVec, // One color set id for every 1-bit in the sampling
}

impl ColexToColorSetMap {

    // sets maps from color set to its index in the distinct color sets
    fn new(sbwt: Arc<SbwtIndex<SubsetMatrix>>, lcs: Option<&LcsArray>, sample_distance: usize, colex_to_color_set_id: Vec<usize>, n_distinct_color_sets: usize, n_threads: usize) -> Self {

        let get_colorset_fn = |colex| colex_to_color_set_id[colex]; // TODO: this actually returns a color set id. Rename here and later.
        let mut sampling_marks = Self::pick_sampled_kmers(sample_distance, &sbwt, lcs, get_colorset_fn, n_threads);

        let color_set_id_bit_width = n_distinct_color_sets.next_power_of_two().trailing_zeros() as usize;
        let mut sampled_color_set_ids = CompactIntVec::new(sampling_marks.count_ones(), color_set_id_bit_width); // In colex order
        let mut n_ids_stored = 0_usize;
        for colex in 0..sbwt.n_sets() {
            if sampling_marks.get(colex) {
                sampled_color_set_ids.set(n_ids_stored, colex_to_color_set_id[colex]);
                n_ids_stored += 1;
            }
        }

        log::info!("Building rank support for sampling marks");
        sampling_marks.enable_rank();

        Self{sbwt, sampling: sampling_marks, color_set_ids: sampled_color_set_ids}
    }

    fn colex_to_color_set_id(&self, mut colex: usize) -> usize {
        if self.sampling.get(colex) {
            // This set is stored
            self.color_set_ids.get(self.sampling.rank(colex)) as usize
        } else {
            // This set is not stored -> walk forward in the de Bruijn graph
            loop {
                for char_idx in 0..self.sbwt.alphabet().len() {
                    if self.sbwt.sbwt().set_contains(colex, char_idx as u8) {
                        // Found the outedge label
                        let new_colex = self.sbwt.lf_step(colex, char_idx);
                        return self.colex_to_color_set_id(new_colex); // Todo: no recursion
                    }
                }

                // No outedges found -> colex is not a suffix group leader position
                assert!(colex > 0);
                colex -= 1;
            }
        }
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        self.sampling.serialize(out).unwrap();
        self.color_set_ids.serialize(out);
    }

    pub fn load(input: &mut impl std::io::Read, sbwt: Arc<SbwtIndex<SubsetMatrix>>) -> Self {
        let sampling = simple_sds_sbwt::bit_vector::BitVector::load(input).unwrap();
        let color_set_ids = CompactIntVec::load(input);

        assert_eq!(sampling.len(), sbwt.n_sets());
        assert_eq!(color_set_ids.len(), sampling.count_ones());

        Self{sbwt: sbwt.clone(), sampling, color_set_ids}
    }

    /// Utility function used in construction
    fn pick_sampled_kmers<'a, F: Fn(usize) -> usize + Sync + Send>(sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>, lcs: Option<&LcsArray>, get_colorset_fn: F, n_threads: usize) -> simple_sds_sbwt::bit_vector::BitVector {
        // Find starts of unitigs. Walk forward to the end of the unitig. Segment by color sets.
        
        let marks = simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
        let marks_mutex = Mutex::new(marks); // Need thread-safe modifications
        let marks_mutex_borrow = &marks_mutex; // Passed into the callback

        let callback = |nodes: &[Node], _: &[u8]| {
            let mut marks = marks_mutex_borrow.lock().unwrap();

            let mut prev_set: Option<usize> = None;
            let mut prev_sample_pos = usize::MAX;
            for (v_pos, v) in nodes.iter().enumerate().rev() {
                let colex = v.id; 
                let cur_set = get_colorset_fn(colex);

                // Sample this node if any of the following are true:
                // - v is the last node of the unitig
                // - v has a different color set than the previous node in iteration order 
                // - v is far enough from the previous sampled node 
                if prev_set.is_none() || cur_set != prev_set.unwrap() || prev_sample_pos - v_pos >= sample_distance {
                    marks.set_bit(colex, true);
                    prev_sample_pos = v_pos;
                }
                prev_set = Some(cur_set);
            }
        };

        log::info!("Initializing the de Bruijn graph");
        let dbg = sbwt::dbg::Dbg::new(sbwt, lcs, n_threads);

        log::info!("Iterating unitigs");
        dbg.iter_unitigs_with_callback(callback, n_threads);

        let marks = marks_mutex.into_inner().unwrap();
        let marks = simple_sds_sbwt::bit_vector::BitVector::from(marks);

        let n_sampled = marks.count_ones();
        log::info!("Sampled {} out of {} k-mers ({:.2}%)", n_sampled, sbwt.n_kmers(), n_sampled as f64 / sbwt.n_kmers() as f64 * 100.0);

        log::info!("Unitig sampling finished");

        marks
    }

}

struct UnitigImportSeqBatch {
    pub concat: Vec<u8>,
    pub starts: Vec<usize>, // Has concat.len() at the end 
    pub color_set_ids: Vec<usize>, // The color set id for each sequence in the concatenation
}

impl UnitigImportSeqBatch {
    fn get_seq(&self, idx: usize) -> &[u8] {
        &self.concat[self.starts[idx]..self.starts[idx+1]]
    }

    fn n_seqs(&self) -> usize {
        self.starts.len() - 1 // Has concat.len() at the end
    }

    fn process(&mut self, results_out: &mut Vec<(usize, usize)>, index: &sbwt::StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>, sample_distance: usize) {
        let k = index.k();

        let mut set_ids_fn = |seq: &[u8]| {
            if seq.len() < k { return }
            let mut distance_from_end = seq.len()-k+1;
            for (start, (match_len, colex_range)) in index.matching_statistics_iter(seq).skip(k-1).enumerate() {
                distance_from_end -= 1;
                if match_len == k {
                    assert_eq!(colex_range.len(), 1);
                    if distance_from_end % (sample_distance-1) == 0 {
                        results_out.push((colex_range.start, color_set_id));
                    }
                } else {
                    panic!( // TODO: return error instead
                        "Error reading unitigs from dump: k-mer {} not found in SBWT", 
                        String::from_utf8_lossy(&seq[start..start+k])
                    );
                }
            }
            assert!(distance_from_end == 0);
        };

        for seq_idx in 0..self.n_seqs() {
            let seq = self.get_seq(seq_idx);
            let color_set_id = self.color_set_ids[seq_idx];

            // Process both forward and reverse complement directions
            set_ids_fn(seq);
            reverse_complement_in_place(&mut seq);
            set_ids_fn(seq);
        }

    }
}

// Todo: the logic associated with this is all over the place. Refactor into one place.
fn unitig_import_parser_thread(unitig_dump: impl std::io::BufRead + Send + 'static, buf_cap: usize, out: Sender<UnitigImportSeqBatch>){
        
    let mut seqs = jseqio::reader::DynamicFastXReader::new(unitig_dump).unwrap();

    let mut cur_concat = Vec::<u8>::with_capacity(buf_cap);
    let mut cur_starts = Vec::<usize>::new();
    let mut cur_color_set_ids = Vec::<usize>::new();
    
    while let Some(rec) = seqs.read_next().unwrap() {
        let color_set_id = index_import::get_color_set_id_from_fasta_header(rec.head);
        cur_color_set_ids.push(color_set_id);

        // Add to concatenation
        cur_starts.push(cur_concat.len());
        cur_concat.extend(rec.seq);

        if cur_concat.len() >= buf_cap {
            cur_starts.push(cur_concat.len()); // End sentinel, as required
            let batch = UnitigImportSeqBatch{concat: cur_concat, starts: cur_starts, color_set_ids: cur_color_set_ids};
            out.send(batch).unwrap();

            // Start a new batch
            cur_concat = Vec::<u8>::with_capacity(buf_cap);
            cur_starts = Vec::<usize>::new();
            cur_color_set_ids.clear();
        }
    }

    if !cur_concat.is_empty() {
        // Send remaining batch
        cur_starts.push(cur_concat.len()); // End sentinel, as required
        let batch = UnitigImportSeqBatch{concat: cur_concat, starts: cur_starts, color_set_ids: cur_color_set_ids};
        out.send(batch).unwrap();
    }

    log::info!("Producer thread: all work pushed to work queue");
    drop(out);
}

impl<CSS: ColorSetStorage> CompactColexKmers<CSS> {

    /// Input: 
    /// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
    ///   color j is present in set i.
    /// TODO: take CompactIntVec for colex_to_color_set_id instead of Vec<usize>
    pub fn new(sbwt: Arc<SbwtIndex<SubsetMatrix>>, lcs: LcsArray, colex_to_color_set_id: Vec<usize>, color_sets: CSS, n_colors: usize, sample_distance: usize, n_threads: usize)
    -> CompactColexKmers<CSS> {
        let colex_map = ColexToColorSetMap::new(sbwt.clone(), Some(&lcs), sample_distance, colex_to_color_set_id, color_sets.n_sets(), n_threads);
        Self {sbwt, lcs, sets: color_sets, map: colex_map}
    }

    pub fn new_from_colored_unitig_dump(
        sbwt: SbwtIndex<SubsetMatrix>, 
        lcs: LcsArray, 
        sample_distance: usize,
        n_threads: usize,
        metadata_dump: impl std::io::BufRead, 
        unitig_dump: impl std::io::BufRead + Send + 'static, 
        color_dump: impl std::io::BufRead) 
        -> Self {

        assert!(sample_distance > 0);

        log::info!("Reading metadata");
        let metadata = index_import::read_index_dump_metadata(metadata_dump);

        let mut reader = jseqio::reader::DynamicFastXReader::new(unitig_dump).unwrap();
        let index = sbwt::StreamingIndex::new(&sbwt, &lcs);

        log::info!("Building (colex, color set id) pairs");
        std::thread::scope(|scope| {
            let (parser_out, worker_in) = unbounded();

            // Create producer
            let producer_handle = scope.spawn(move || {
                unitig_import_parser_thread(unitig_dump, 1 << 20, parser_out);
            });

            // Create workers
            let mut worker_handles = Vec::<std::thread::ScopedJoinHandle::<()>>::new();

            for _ in 0..n_threads {
                let worker_in_clone = worker_in.clone();
                worker_handles.push(scope.spawn(move || {
                    let mut colex_to_color_set_id: Vec<(usize, usize)> = vec![]; // (colex, coled_set_id)
                    while let Ok(mut batch) = worker_in_clone.recv(){
                        batch.process();
                    }
                }));
            }
        });

        let n_colors = metadata.num_colors;

        log::info!("Sorting (colex, color set id) pairs");
        colex_to_color_set_id.sort(); // Sorts by colex
        let bit_width = metadata.num_color_sets.next_power_of_two().trailing_zeros() as usize;
        log::info!("Building compressed representation for color set ids");
        let mut stored_color_set_ids = CompactIntVec::new(colex_to_color_set_id.len(), bit_width);
        let mut sample_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
        for (rank, (colex, id)) in colex_to_color_set_id.into_iter().enumerate() {
            stored_color_set_ids.set(rank, id);
            sample_marks.set_bit(colex, true);
        }
        let mut sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(sample_marks);
        sample_marks.enable_rank();

        log::info!("Reading distinct color sets");
        let color_set_stream = index_import::ColorSetDumpIterGenerator::new(color_dump);
        let distinct_css = CSS::new(color_set_stream, n_colors);
        let distinct_css = *distinct_css; // Unbox

        let sbwt = Arc::new(sbwt);
        let colex_map = ColexToColorSetMap {
            sbwt: sbwt.clone(), // Clones the Arc, not the sbwt
            sampling: sample_marks,
            color_set_ids: stored_color_set_ids,
        };

        Self {sbwt, lcs, sets: distinct_css, map: colex_map}
    }


    pub fn new_single_colored(sbwt: Arc<SbwtIndex<SubsetMatrix>>, lcs: LcsArray, sample_distance: usize, n_threads: usize) -> Self {
        let n_colors = 1;
        let int_bitwidth = 1;

        // The only color set is {0}
        let vv = VecVecUsizeIteratorGenerator{sets: vec![vec![0]], pos: 0};
        let sets = CSS::new(vv, n_colors);

        log::info!("Sampling nodes");
        let mut unitig_samples = ColexToColorSetMap::pick_sampled_kmers(sample_distance, &sbwt, Some(&lcs), |_colex| 0, n_threads);
        unitig_samples.enable_rank();
        log::info!("Storing color set ids for sampled nodes");
        //let color_set_ids = IntVector::with_len(unitig_samples.count_ones(), int_bitwidth, 0).unwrap();
        let color_set_ids = CompactIntVec::new(unitig_samples.count_ones(), int_bitwidth);
        let colex_map = ColexToColorSetMap{
            sbwt: sbwt.clone(),
            sampling: unitig_samples,
            color_set_ids,
        };
        Self {sbwt, lcs, sets: *sets, map: colex_map}
    }

    pub fn colex_to_set_id(&self, colex: usize) -> usize {
        self.map.colex_to_color_set_id(colex)
    }

    pub fn set_id_to_set<'a>(&'a self, id: usize) -> CSS::SetView<'a> {
        self.sets.get_set_view(id)
    }

    pub fn colex_to_set<'a>(&'a self, colex: usize) -> CSS::SetView<'a> {
        self.sets.get_set_view(self.colex_to_set_id(colex))
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();
        self.sets.serialize(out);
        self.map.serialize(out);
    }

    /// If this struct is going to be merged with [crate::coloring::merge_colorings], it will
    /// need select support on the sbwt. We need to build it already during loading because
    /// once the sbwt is put on to the heap into an Arc, it cannot be modified anymore.
    /// Unless we make it an Arc<Refcell<...>>, but that might have overhead because then
    /// it will do run-time borrow checking on every access if I understand correctly.
    pub fn load(input: &mut impl std::io::Read, enable_select: bool) -> Self {
        let mut sbwt = SbwtIndex::<SubsetMatrix>::load(input).unwrap();
        if enable_select {
            log::info!("Building select support");
            sbwt.build_select();
        }
        let sbwt = Arc::new(sbwt);
        let lcs = LcsArray::load(input).unwrap();
        let sets = CSS::load(input);
        let map = ColexToColorSetMap::load(input, sbwt.clone());
        CompactColexKmers{sbwt, lcs, sets, map}
    }

    pub fn lookup_kmer_color_sets(&self, seq: &[u8]) -> Vec<Option<CSS::SetView<'_>>> {
        let k = self.sbwt.k();
        if seq.len() < k {
            return vec![];
        }

        let si = sbwt::StreamingIndex::new(&self.sbwt, &self.lcs);

        let mut set_views = Vec::<Option<CSS::SetView<'_>>>::with_capacity(seq.len()-k+1);
        let mut prev_set_id: Option<usize> = None;
        for (len, range) in si.matching_statistics_iter(seq).skip(k-1) {
            if len == k {
                assert!(range.len() == 1);
                let colex = range.start;
                let set_id = self.colex_to_set_id(colex);
                if prev_set_id.is_some_and(|p| p == set_id) {
                    // Same as previous
                    let prev = set_views.last().unwrap();
                    set_views.push(prev.clone());
                } else {
                    set_views.push(Some(self.set_id_to_set(set_id)));
                }
                prev_set_id = Some(set_id);
            } else {
                prev_set_id = None;
            }
        }

        set_views
    }

    pub fn get_k(&self) -> usize {
        self.sbwt.k()
    }

    pub fn get_set_storage(&self) -> &CSS {
        &self.sets
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct BitKey<'a> { // Bitslice with a custom hash function
    pub bits: &'a BitSlice,
}

impl Hash for BitKey<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash 64 bits at a time
        let len = self.bits.len();
        let n_words = len.div_ceil(64);
        for i in 0..n_words {
            let start = 64*i;
            let end = min(64*(i+1), len);
            let word: u64 = self.bits[start..end].load();
            word.hash(state);
        }
        len.hash(state);  // include length to distinguish e.g. 0b1 from 0b10
    }
}


fn figure_out_if_we_need_to_sample_nonsampled_vs_absent(
    absent_sbwt: &SbwtIndex<SubsetMatrix>, 
    mut absent_colex: usize, // Position in the absent sbwt where k-mer would be inserted
    merged_colex: usize,
    merged_leader_marks: &bitvec::vec::BitVec<u64, Lsb0>,
    absent_merge_marks: &bitvec::vec::BitVec<u64, Lsb0>) -> bool {

    // This node may become the end of a colored unitig in the merged graph. So we may need
    // to sample it. 
    // 
    // This happens if any of the following happen:
    //   (i)   The merged graph has a new outneighbor for this node (unitig ends).
    //   (ii)  The current outneighbor gets a new in-neighbor (unitig ends).
    //   (iii) There will be an edge from the node in the present SBWT to a node
    //         in the absent SBWT. Then the node from the absent SBWT may introduce 
    //         a new color, in which case the colored unitig ends.
    //
    // We assume that all color sets are non-empty, which means that if there is an
    // outedge into the absent sbwt, then this always introduces a new color in case (iii).
    // Under this assumption, if case (i) or (ii) happens, case (iii) also happens, so it's enough
    // to check only for case (iii). If our assumption that all color sets are nonempty
    // does not hold, it only means that we may sample a node unnecessarily, but the
    // color set structure is still correct. 

    assert!(!absent_merge_marks[merged_colex]); // Should be absent
    let mut s = merged_colex;
    while !merged_leader_marks[s] {
        // merged_leader_marks[0] is always set so s > 0 if we are here
        s -= 1;
        if absent_merge_marks[s] {
            absent_colex -= 1;
        }
    }
    let mut e = merged_colex+1;
    while e < merged_leader_marks.len() && !merged_leader_marks[e] {
        e += 1;
    }

    // [s..e) is the suffix group of the present k-mer in the merged sbwt.
    for i in s..e {
        if absent_merge_marks[i] {
            // Suffix group leader in the absent sbwt
            for c_idx in 0..absent_sbwt.alphabet().len() {
                if absent_sbwt.sbwt().set_contains(absent_colex, c_idx as u8) {
                    return true; // Sample x
                }
            }
            return false; // Suffix group leader did not have any edge
        }
    }
    false
}

struct PartitionedIdMap {
    #[allow(clippy::type_complexity)]
    hashmaps: Vec<HashMap::<(Option::<usize>, Option::<usize>), usize>>,
}

struct PartitionedReadOnlyIdMap {
    #[allow(clippy::type_complexity)]
    hashmaps: Vec<HashMap::<(Option::<usize>, Option::<usize>), usize>>,
    cumul_sizes: Vec<usize> // index i contains total length of hash maps [0..i)
}

impl PartitionedIdMap {
    fn hash_pair(x: (Option<usize>, Option<usize>)) -> u64 {
        let mut hasher = FxHasher::default();
        x.hash(&mut hasher);
        hasher.finish()
    }

    fn insert_pair(&mut self, x: (Option<usize>, Option<usize>)) {
        let r = Self::hash_pair(x);
        let hash_map_idx = (r / (u64::MAX / self.hashmaps.len() as u64)) as usize;
        let H = &mut self.hashmaps[hash_map_idx];
        if !H.contains_key(&x) {
            H.insert(x, H.len());
        }
    }

    // There is not method to get a pair. For that, first convert the struct
    // into PartitionedReadOnlyIdMap, which does some precalc to make the
    // lookup faster.
}

impl PartitionedReadOnlyIdMap {
    fn new(old: PartitionedIdMap) -> Self {
        let mut cumul_sizes = Vec::<usize>::with_capacity(old.hashmaps.len() + 1); 
        cumul_sizes.push(0);
        old.hashmaps.iter().for_each(|H| {
            cumul_sizes.push(cumul_sizes.last().unwrap() + H.len());
        });

        PartitionedReadOnlyIdMap{hashmaps: old.hashmaps, cumul_sizes}
    }

    fn total_len(&self) -> usize {
        self.cumul_sizes[self.hashmaps.len()]
    }

    fn get_new_id_of_pair(&self, x: (Option<usize>, Option<usize>)) -> usize {
        let r = PartitionedIdMap::hash_pair(x);
        let hash_map_idx = (r / (u64::MAX / self.hashmaps.len() as u64)) as usize;
        self.cumul_sizes[hash_map_idx] + self.hashmaps[hash_map_idx][&x]
    }

    #[allow(clippy::type_complexity)]
    fn get_old_ids_sorted_by_new_id(&self) -> Vec<(usize, (Option::<usize>, Option::<usize>))> {
        // Collect all elements (new id, old id pair) from the hash maps
        let n_pairs_total = self.total_len();
        let mut id_pairs_in_new_id_order = self.hashmaps.iter().fold(
            Vec::<(usize, (Option::<usize>, Option::<usize>))>::with_capacity(n_pairs_total),
            |mut acc, H| {
                let len_before = acc.len();
                acc.extend(
                    H.iter().map(|(pair, new_id)| (*new_id + len_before, *pair))
                );
                acc
            }
        );
        id_pairs_in_new_id_order.sort();
        id_pairs_in_new_id_order
    }

}


fn compute_color_id_pairs_and_merged_unitig_sampling<CSS: ColorSetStorage>(coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, lcs1: &LcsArray, lcs2: &LcsArray, merge_plan: &MergeInterleaving, n_threads: usize) -> (PartitionedReadOnlyIdMap, simple_sds_sbwt::raw_vector::RawVector) {
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();

    // Distinct color id pairs, inserted into key-disjoint hash maps. The values are
    // color set ids within the hash map.
    let hashmaps = vec![HashMap::<(Option::<usize>, Option::<usize>), usize>::new(); n_threads];
    let mut new_id_map = PartitionedIdMap{hashmaps};

    let mut colex1 = 0_usize;
    let mut colex2 = 0_usize;

    let mut color_set_sample_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(merged_len, false);
    log::info!("Building DBG support");
    let dbg1 = sbwt::dbg::Dbg::new(&(*coloring1.map.sbwt), Some(lcs1), n_threads);
    let dbg2 = sbwt::dbg::Dbg::new(&(*coloring2.map.sbwt), Some(lcs2), n_threads);
    let mut outlabel_buf_1 = Vec::<u8>::new();
    let mut outlabel_buf_2 = Vec::<u8>::new();

    log::info!("Computing new color set id pairs and merged unitig sampling");
    #[derive(Debug)]
    enum Case { // Three cases in a loop below
        Sampled(usize),
        NotSampled,
        Absent,
    }
    for merged_colex in 0..merged_len {
        if !merge_plan.is_dummy[merged_colex] {
            let c1 = if !merge_plan.s1[merged_colex] {
                Case::Absent
            } else if coloring1.map.sampling.get(colex1) {
                Case::Sampled(coloring1.colex_to_set_id(colex1))
            } else {
                Case::NotSampled
            };

            let c2 = if !merge_plan.s2[merged_colex] {
                Case::Absent
            } else if coloring2.map.sampling.get(colex2) {
                Case::Sampled(coloring2.colex_to_set_id(colex2))
            } else {
                Case::NotSampled
            };

            // Ok, this is going to get a bit verbose but bear with me. We have
            // 3 * 3 = 9 cases. There are two symmetric pairs of cases and three unique cases. We could
            // reduce code duplication by making symmetric cases call a common function,
            // but it's so few lines of code anyway so let's just go with this.
            match (c1, c2) {
                (Case::Sampled(id1), Case::Sampled(id2)) => {
                    new_id_map.insert_pair((Some(id1), Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::Sampled(id1), Case::NotSampled) => {
                    let id2 = coloring2.colex_to_set_id(colex2);
                    new_id_map.insert_pair((Some(id1), Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::Sampled(id1), Case::Absent) => {
                    new_id_map.insert_pair((Some(id1), None));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::NotSampled, Case::Sampled(id2)) => {
                    //eprintln!("Case 3");
                    let id1 = coloring1.colex_to_set_id(colex1);
                    new_id_map.insert_pair((Some(id1), Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::NotSampled, Case::NotSampled) => {
                    // K-mer is in both SBWTs but its not sampled in either one.
                    // Since it is not sampled in either SBWT, the outdegree of this k-mer
                    // is 1 in both. But we might still need to sample it in the merged graph.
                    // There are two cases:
                    // 1) The outneighbor k-mers are the same k-mer. Then the outdegree in the merged graph 
                    //    will be 1, and that outneighbor will have the same color set id pair as this
                    //    one -> this node does not need to be sampled
                    // 2) The outneighbor k-mers are different. Now we have a new outgoing branch at this 
                    //    node. Which means this node needs to be sampled.

                    outlabel_buf_1.clear();
                    outlabel_buf_2.clear();
                    dbg1.push_outlabels(Node{id: colex1}, &mut outlabel_buf_1);
                    dbg2.push_outlabels(Node{id: colex2}, &mut outlabel_buf_2);
                    assert_eq!(outlabel_buf_1.len(), 1);
                    assert_eq!(outlabel_buf_2.len(), 1);
                    //eprintln!("{} {}", *outlabel_buf_1.first().unwrap() as char, *outlabel_buf_2.first().unwrap() as char);
                    match (outlabel_buf_1.first(), outlabel_buf_2.first()) {
                        (Some(a), Some(b)) => {
                            if a != b { // Case 2 in the comment above
                                color_set_sample_marks.set_bit(merged_colex, true);
                                let id1 = coloring1.colex_to_set_id(colex1);
                                let id2 = coloring2.colex_to_set_id(colex2);
                                new_id_map.insert_pair((Some(id1), Some(id2)));
                            } else { // The else-branch would be case 1 but then there is nothing to do

                            }
                        }
                        _ => panic!("Bug at computing color set samples bit vector in merge") // Both should have outdegree > 0
                    }
                },
                (Case::NotSampled, Case::Absent) => {
                    let id1 = coloring1.colex_to_set_id(colex1);
                    new_id_map.insert_pair((Some(id1), None));
                    if figure_out_if_we_need_to_sample_nonsampled_vs_absent(&coloring2.map.sbwt, colex2, merged_colex, &merge_plan.is_leader, &merge_plan.s2) {
                        color_set_sample_marks.set_bit(merged_colex, true);
                    }
                },
                (Case::Absent, Case::Sampled(id2)) => {
                    new_id_map.insert_pair((None, Some(id2)));
                    color_set_sample_marks.set_bit(merged_colex, true);
                },
                (Case::Absent, Case::NotSampled) => {
                    let id2 = coloring2.colex_to_set_id(colex2);
                    new_id_map.insert_pair((None, Some(id2)));
                    if figure_out_if_we_need_to_sample_nonsampled_vs_absent(&coloring1.map.sbwt, colex1, merged_colex, &merge_plan.is_leader, &merge_plan.s1) {
                        color_set_sample_marks.set_bit(merged_colex, true);
                    }
                },
                (Case::Absent, Case::Absent) => panic!("Nonexisting merged kmer") // Impossible
            }
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s2[merged_colex] as usize;
    }

    (PartitionedReadOnlyIdMap::new(new_id_map), color_set_sample_marks)

}

struct TwoSetMerger<L: Iterator<Item = usize>, R: Iterator<Item = usize>> {
    left: Option<L>,
    right: Option<R>,
    left_n_colors: usize, // The left set get colors 0..left_n_colors, the right set gets left_n_colors..
}

impl<L: Iterator<Item = usize>, R: Iterator<Item = usize>> Iterator for TwoSetMerger<L,R> {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        // Terrible branch city. TODO: do better.

        // Try to take from left
        if let Some(l) = &mut self.left {
            if let Some(x) = l.next() {
                return Some(x);
            } 
        }

        // Could not take from left -> take from right
        if let Some(r) = &mut self.right {
            r.next().map(|x| self.left_n_colors + x)
        } else {
            None // Finished
        }
    }
}


fn encode_merged_color_sets<CSS: ColorSetStorage>(new_id_map: &PartitionedReadOnlyIdMap, coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>) -> CSS {

    let n_colors_1 = coloring1.sets.get_full_set().iter().count();
    let n_colors_2 = coloring2.sets.get_full_set().iter().count();
    let id_pairs_in_new_id_order = new_id_map.get_old_ids_sorted_by_new_id();

    // Create an iterator of combined sets
    let mut pair_id = 0_usize;
    let n_pairs = id_pairs_in_new_id_order.len();
    let n_colors_1_ref = &n_colors_1; // Reference to move by reference into the closure 
    let iter_of_iters = std::iter::from_fn(move || {
        if pair_id == n_pairs {
            None
        } else {
            let (_, (left, right)) = id_pairs_in_new_id_order[pair_id]; 
            pair_id += 1;

            match (left,right) {
                (Some(x), Some(y)) => {
                    let set1 = coloring1.set_id_to_set(x);
                    let set2 = coloring2.set_id_to_set(y);
                    Some(TwoSetMerger{left: Some(set1.iter()), right: Some(set2.iter()), left_n_colors: *n_colors_1_ref})
                },
                (Some(x), None) => {
                    let set1 = coloring1.set_id_to_set(x);
                    Some(TwoSetMerger{left: Some(set1.iter()), right: None, left_n_colors: *n_colors_1_ref})
                },
                (None, Some(y)) => {
                    let set2 = coloring2.set_id_to_set(y);
                    Some(TwoSetMerger{left: None, right: Some(set2.iter()), left_n_colors: *n_colors_1_ref})
                }
                (None, None) => panic!("Nonexisting color set id pair")
            }
        }
    });

    *CSS::new_from_iter_of_iters(iter_of_iters, n_colors_1 + n_colors_2)

}

fn store_new_sampled_color_ids<CSS: ColorSetStorage>(n_distinct_color_sets: usize, merge_plan: &MergeInterleaving, color_set_sample_marks: &simple_sds_sbwt::bit_vector::BitVector, coloring1: &CompactColexKmers<CSS>, coloring2: &CompactColexKmers<CSS>, pair_to_new_id_maps: &PartitionedReadOnlyIdMap) -> CompactIntVec {
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();

    let bits_per_color_set_id = n_distinct_color_sets.next_power_of_two().trailing_zeros() as usize;
    let mut sampled_ids = CompactIntVec::new(color_set_sample_marks.count_ones(), bits_per_color_set_id);
    let mut n_items_pushed = 0_usize;
    let mut colex1 = 0_usize;
    let mut colex2 = 0_usize;
    for merged_colex in 0..merged_len {
        if color_set_sample_marks.get(merged_colex) {
            let color_set_id_1 = if merge_plan.s1[merged_colex] {
                Some(coloring1.colex_to_set_id(colex1))
            } else {
                None
            };

            let color_set_id_2 = if merge_plan.s2[merged_colex] {
                Some(coloring2.colex_to_set_id(colex2))
            } else {
                None
            };

            // The merge plan should not have a zero-bit at the same position in s1 and s2
            assert!(color_set_id_1.is_some() || color_set_id_2.is_some());
            let id = pair_to_new_id_maps.get_new_id_of_pair((color_set_id_1, color_set_id_2));
            sampled_ids.set(n_items_pushed, id);
            n_items_pushed += 1;
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s2[merged_colex] as usize;
    }

    sampled_ids
}

pub fn merge_compact_colorings<CSS: ColorSetStorage>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, optimize_peak_ram: bool, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::MergeInterleaving::new(&(*coloring1.map.sbwt), &(*coloring2.map.sbwt), optimize_peak_ram, n_threads);

    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();    

    log::info!("Computing color id pairs and merged sampling");
    let (new_id_map, color_set_sample_marks) = compute_color_id_pairs_and_merged_unitig_sampling(&coloring1, &coloring2, &coloring1.lcs, &coloring2.lcs, &merge_plan, n_threads);

    let mut color_set_sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(color_set_sample_marks);
    color_set_sample_marks.enable_rank();
    let n_sampled = color_set_sample_marks.rank(color_set_sample_marks.len());
    log::info!("Sampled {} out of {} SBWT nodes ({:.2}%)", n_sampled, merged_len, n_sampled as f64 / merged_len as f64 * 100.0);

    log::info!("Encoding distinct merged color sets");
    let css = encode_merged_color_sets(&new_id_map, &coloring1, &coloring2);

    log::info!("Storing new sampled color set ids");
    let n_distinct_color_sets = new_id_map.total_len(); 
    let sampled_ids = store_new_sampled_color_ids(n_distinct_color_sets, &merge_plan, &color_set_sample_marks, &coloring1, &coloring2, &new_id_map);

    log::info!("Interleaving SBWTs");
    let precalc_len = max(coloring1.map.sbwt.get_lookup_table().prefix_length, coloring2.map.sbwt.get_lookup_table().prefix_length);

    let sbwt1 = (*coloring1.map.sbwt).clone(); // Todo: avoid clone. Currently unavoidable because we have just a reference to the SBWT, but the merge needs an owned value.
    drop(coloring1);

    let sbwt2 = (*coloring2.map.sbwt).clone(); // Todo: avoid clone
    drop(coloring2);

    log::info!("Interleaving SBWTs");
    let merged_sbwt = Arc::new(sbwt::merge(sbwt1, sbwt2, merge_plan, precalc_len, n_threads));

    log::info!("Computing the merged LCS array"); // Todo: could we do this during the interleave?
    let merged_lcs = LcsArray::from_sbwt(&merged_sbwt, n_threads);

    let new_coloring = CompactColexKmers { 
        sbwt: merged_sbwt.clone(),
        lcs: merged_lcs,
        sets: css, 
        map: ColexToColorSetMap {
            sbwt: merged_sbwt.clone(), 
            sampling: color_set_sample_marks, 
            color_set_ids: sampled_ids 
        }
    };

    log::info!("Color merge finished");
    new_coloring

}




/// Output:
/// - Distinct color sets encoded as ColorSetStorage
/// - HashMap from color set to its index in ColorSets
pub fn hash_and_encode_distinct_sets<'a, CSS: ColorSetStorage>(colex_to_set: &'a CSS, n_colors: usize) -> (CSS, HashMap::<CSS::SetView<'a>, usize, BuildHasherDefault::<FxHasher>>) {
    let n_sets = colex_to_set.n_sets();

    log::info!("Hashing distinct color sets");

    let mut distinct_sets = HashMap::<CSS::SetView<'a>, usize, BuildHasherDefault::<FxHasher>>::default(); // Set -> id
    let mut distinct_set_colex_ranks = Vec::<usize>::new();
    let bar = indicatif::ProgressBar::new(n_sets as u64);
    for colex in 0..n_sets {
        let key = colex_to_set.get_set_view(colex);
        if !distinct_sets.contains_key(&key) {
            distinct_sets.insert(key, distinct_sets.len());
            distinct_set_colex_ranks.push(colex);
        }
        if colex % 1000 == 0 {
            bar.inc(1000);
        }
    }
    bar.finish();

    log::info!("{} distinct color sets found", distinct_sets.len());

    // Create an iterator of iterators, each inner iterator iterating over one color set
    let color_sets_iterator = distinct_set_colex_ranks.into_iter().map(|colex| {
        colex_to_set.get_set_view(colex).iter()
    });

    let colorsets = CSS::new_from_iter_of_iters(color_sets_iterator, n_colors);

    (*colorsets, distinct_sets)

}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use jseqio::seq_db::SeqDB;
    use sbwt::{BitPackedKmerSortingMem, LcsArray, SbwtIndex, SubsetMatrix};
    use simple_sds_sbwt::ops::BitVec;

    use crate::{bitmap_storage::build_from_seq_dbs, build_coloring, colex_colored_kmers::hash_and_encode_distinct_sets, coloring_interface::{ColorSetStorage, ColorSetView}, sparse_dense_storage::SparseDenseStorage};

    use super::{CompactColexKmers, merge_compact_colorings};


    #[cfg(test)]
    pub(crate) fn gen_random_dna_string(len: usize, seed: u64) -> Vec<u8> {
        use rand_chacha::rand_core::{RngCore, SeedableRng};

        let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(seed);
        (0..len).map(|_| { 
            match rng.next_u64() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => panic!("Impossible")
            }
        }).collect()
    }

    fn build_color_sets<CSS: ColorSetStorage>(sbwt1: &SbwtIndex<SubsetMatrix>, lcs1: &LcsArray, dbs1: Vec<SeqDB>, n_threads: usize) 
    -> (Vec<usize>, CSS){
        let n_colors_1 = dbs1.len();
        let bms1 = build_from_seq_dbs(dbs1, &sbwt1, &lcs1, n_threads);

        let iter_of_iters_1 = (0..sbwt1.n_sets()).into_iter().map(|colex| bms1.get_set_view(colex).iter());
        let colex_to_css_1 = *CSS::new_from_iter_of_iters(iter_of_iters_1, n_colors_1);

        let (distinct_css_1, set_to_id_1) = hash_and_encode_distinct_sets(&colex_to_css_1, n_colors_1);
        let colex_to_id: Vec<usize> = (0..sbwt1.n_sets()).into_iter().map(|colex| {
            set_to_id_1[&colex_to_css_1.get_set_view(colex)]
        }).collect(); 

        (colex_to_id, distinct_css_1)
    }

    #[test]
    fn test_merge() {

        if std::env::var("RUST_LOG").is_err() {
            std::env::set_var("RUST_LOG", "info")
        }
        env_logger::init();

        let n_threads = 3;

        for k in 3_usize..10_usize { // k < 3 does not work because construction uses 3-mer binning.

            let input_seqs_1: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (i + k.pow(4)) as u64)).collect();
            let input_seqs_2: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (123456 + i + k.pow(4)) as u64)).collect();

            let mut all_input_seq_slices = Vec::<&[u8]>::new();
            all_input_seq_slices.extend(input_seqs_1.iter().map(|s| s.as_slice()));
            all_input_seq_slices.extend(input_seqs_2.iter().map(|s| s.as_slice()));

            let mut dbs1 = Vec::<SeqDB>::new();
            let mut dbs2 = Vec::<SeqDB>::new();
            let mut dbs_both = Vec::<SeqDB>::new();
            for seq in input_seqs_1.iter() {
                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs1.push(db);

                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs_both.push(db);
            }
            for seq in input_seqs_2.iter() {
                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs2.push(db);

                let mut db = SeqDB::new();
                db.push_seq(seq);
                dbs_both.push(db);
            }

            let (mut sbwt1, lcs1) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(false)
                .k(k)
                .build_lcs(true)
                .n_threads(3)
                .precalc_length(5)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
            .run_from_vecs(&input_seqs_1);

            let (mut sbwt2, lcs2) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(false)
                .k(k)
                .build_lcs(true)
                .n_threads(3)
                .precalc_length(5)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
            .run_from_vecs(&input_seqs_2);

            let (mut sbwt_both, lcs_both) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(false)
                .k(k)
                .build_lcs(true)
                .n_threads(3)
                .precalc_length(5)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
            .run_from_slices(&all_input_seq_slices);

            sbwt1.build_select();
            sbwt2.build_select();
            sbwt_both.build_select();

            let sbwt1 = Arc::new(sbwt1);
            let sbwt2 = Arc::new(sbwt2);
            let sbwt_both = Arc::new(sbwt_both);

            let lcs1 = lcs1.unwrap();
            let lcs2 = lcs2.unwrap();
            let lcs_both = lcs_both.unwrap();


            let sample_distance = 3;

            let (colex_to_id_1, storage_1) = build_color_sets::<SparseDenseStorage>(&sbwt1, &lcs1, dbs1, n_threads); 
            let (colex_to_id_2, storage_2) = build_color_sets::<SparseDenseStorage>(&sbwt2, &lcs2, dbs2, n_threads); 
            let (colex_to_id_both, storage_both)= build_color_sets::<SparseDenseStorage>(&sbwt_both, &lcs_both, dbs_both, n_threads); 
            
            let ccc1 = CompactColexKmers::new(sbwt1, lcs1, colex_to_id_1, storage_1, input_seqs_1.len(), sample_distance, n_threads);
            let ccc2 = CompactColexKmers::new(sbwt2, lcs2, colex_to_id_2, storage_2, input_seqs_2.len(), sample_distance, n_threads);
            let ccc_both = CompactColexKmers::new(sbwt_both, lcs_both, colex_to_id_both, storage_both, input_seqs_1.len() + input_seqs_2.len(), sample_distance, n_threads);

            let ccc_merged = merge_compact_colorings(ccc1, ccc2, true, n_threads);
            let sbwt_merged = &ccc_merged.sbwt;

            for colex in 0..ccc_both.sbwt.n_sets() {
                let kmer = ccc_both.sbwt.access_kmer(colex);

                if kmer.iter().all(|c| *c != b'$') { // Not a dummy k-mer
                    let true_colors: Vec<usize> = ccc_both.colex_to_set(colex).iter().collect();
                    let range = sbwt_merged.search(&kmer).unwrap();
                    assert_eq!(range.len(), 1);
                    let colex_merged = range.start;
                    //let merged_colors = ccc_merged.colex_to_set(colex_merged).as_bitvec(ccc_both.n_colors);
                    let merged_colors: Vec<usize> = ccc_merged.colex_to_set(colex_merged).iter().collect();

                    eprintln!("{} {} {:?} {:?} {} {}", colex, String::from_utf8_lossy(&kmer), true_colors, sbwt_merged.search(&kmer), ccc_merged.map.sampling.get(colex_merged), ccc_merged.colex_to_set_id(colex_merged));
                    assert_eq!(true_colors, merged_colors);
                }

            }
        }
    }
}