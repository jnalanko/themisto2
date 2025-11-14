use bitvec::order::Lsb0;
use bitvec::{field::BitField, slice::BitSlice};
use bitvec::bitvec;
use clap::builder::styling::Color;
use sbwt::MergeInterleaving;
use sbwt::LcsArray;
use sbwt::{dbg::Node, SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::raw_vector::RawVector;
use simple_sds_sbwt::serialize::Serialize;
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Push, Rank, Resize, Vector}, raw_vector::{AccessRaw, PushRaw}};
use rustc_hash::FxHasher;
use std::cmp::max;
use std::sync::Arc;
use std::{cmp::min, collections::HashMap, hash::BuildHasherDefault, sync::Mutex};
use std::hash::{Hash, Hasher};

use crate::coloring_interface::{self, ColorSetView};

/// This is the main data structure in this file: a set of compressed color sets, and a mapping
/// from SBWT colex ranks to color sets such that we can look up the color set of a k-mer by its
/// colex rank in the SBWT. 
pub struct CompactColexColoring {
    // This is on the heap to allow map to refer to it (otherwise assuring lifetime 
    // guarantees becomes problematic). It's reference counted because this struct
    // will have two references to it, the one in self.sbwt, and one in self.map.sbwt.
    // Note that this means that if we replace sbwt here with a new Arc pointing to a new
    // sbwt, then, the map will continue to point to the old sbwt. So don't do that!
    // It's atomic (Arc) because we want to pass this struct to multiple threads.
    sbwt: Arc<SbwtIndex<SubsetMatrix>>, 

    lcs: LcsArray,
    sets: ColorSets, // Distinct color sets
    map: ColexToColorSetMap, // A mapping from the colex rank of a k-mer in the SBWT into a color set id in `sets`
}

/// A data structure for storing arbitary set of sets of integers, such that dense
/// sets are encoded as bitmaps, and sparse sets as lists of integers.
pub struct ColorSets {
    dense_sets: BitMaps,
    sparse_sets: IntVecs,
    n_colors: usize,
    is_dense_marks: simple_sds_sbwt::bit_vector::BitVector, // Has rank support.
}

/// A data structure that stores the color set ids for a subset of sampled k-mers in the SBWT such that
/// the color sets of the rest can be obtained by walking forward in the de Bruijn graph to the
/// closest sampled node.
pub struct ColexToColorSetMap {

    // See the comments inside CompactcolexColoring
    sbwt: Arc<SbwtIndex<SubsetMatrix>>,

    sampling: simple_sds_sbwt::bit_vector::BitVector, // Marks colex ranks that have a color set stored. Has rank support.
    color_set_ids: IntVector, // One color set id for every 1-bit in the sampling
}

// A set of lists of integers, stored in concatenated form.
struct IntVecs {
    intvec_data: IntVector, // Concatenation of IntVecs

    // Ends of individual intvecs, such that ends[0] = 0 and ends[i+1] is the
    // exclusive end of the i-th vector.
    ends: Vec<usize>, 
}

#[derive(Copy, Clone)]
pub struct IntVecSlice<'a> {
    vec: &'a IntVector,
    start: usize,
    end: usize, // Exclusive end
}

// A set of sets encoded as bitmaps.
struct BitMaps {
    bitmap_data: bitvec::vec::BitVec, // Concatenation of bit vectors
    individual_length: usize, // Length of each bitmap in bitmap_data
}

// This enum is only for passing references to individual sets around. The actual
// sets are stored in concatenated form somewhere else in memory. 
#[derive(Copy, Clone)]
pub enum ColorSet<'a> {
    Dense(&'a BitSlice),
    Sparse(IntVecSlice<'a>),
}

impl ColorSet<'_> {

    pub fn extract_and_push_colors_to(&self, buf: &mut Vec<usize>) {
        match self {
            ColorSet::Dense(bv) => {
                for i in bv.iter_ones() {
                    buf.push(i);
                }
            },
            ColorSet::Sparse(iv) => {
                for i in iv.start..iv.end {
                    buf.push(iv.vec.get(i) as usize);
                }
            },
        }
    }

    // Number of elements in the set
    pub fn len(&self) -> usize {
        match self {
            ColorSet::Dense(bv) => {
                bv.count_ones()
            },
            ColorSet::Sparse(iv) => {
                iv.end - iv.start
            },
        }
    }

    pub fn as_bitvec(&self, n_colors: usize) -> bitvec::vec::BitVec {
        match self {
            ColorSet::Dense(bv) => {
                (*bv).into()
            },
            ColorSet::Sparse(iv) => {
                let mut bv = bitvec![0; n_colors];
                for i in iv.start..iv.end {
                    bv.set(iv.vec.get(i) as usize, true);
                }
                bv
            },
        }
    }

    pub fn as_intvec(&self) -> Vec<usize> {
        let mut buf = Vec::<usize>::with_capacity(self.len());
        self.extract_and_push_colors_to(&mut buf);
        buf
    }

}

fn is_dense(bv: &BitSlice) -> bool {
    let n_colors = bv.len();
    let n_elements = bv.count_ones();
    let bits_per_color = n_colors.next_power_of_two().trailing_zeros() as usize;
    let bitmap_size = n_colors;
    let intvec_size = n_elements * bits_per_color;

    bitmap_size <= intvec_size
}

impl IntVecs {
    fn new(bit_width: usize) -> Self {
        IntVecs{intvec_data: IntVector::new(bit_width).unwrap(), ends: vec![0]}
    }

    fn push(&mut self, set: impl IntoIterator<Item = usize>) { // Pushes a new set of integers
        for x in set {
            self.intvec_data.push(x as u64);
        }
        self.ends.push(self.intvec_data.len());
    }

    fn shrink_to_fit(&mut self) {
        self.intvec_data.resize(self.intvec_data.len(), 0);
    }

    fn get(&self, vec_idx: usize) -> IntVecSlice {
        IntVecSlice{vec: &self.intvec_data, start: self.ends[vec_idx], end: self.ends[vec_idx+1]}
    }

    fn n_sets(&self) -> usize {
        self.ends.len() - 1 // Minus 1 because there is a 0 at the start of ends
    }

    fn serialize(&self, out: &mut impl std::io::Write) {
        // Serialize using bincode
        self.intvec_data.serialize(out).unwrap();
        bincode::serialize_into(out, &self.ends).unwrap();
    }

    fn load(input: &mut impl std::io::Read) -> Self {
        // Deserialize using bincode
        let intvec_data = IntVector::load(input).unwrap();
        let ends: Vec<usize> = bincode::deserialize_from(input).unwrap();
        assert!(!ends.is_empty() && ends[0] == 0); // The first end must be 0
        IntVecs{intvec_data, ends}
    }

}


impl BitMaps {
    fn new(individual_length: usize) -> Self {
        BitMaps{bitmap_data: bitvec::vec::BitVec::new(), individual_length}
    }

    fn push(&mut self, bv: &bitvec::slice::BitSlice) {
        assert_eq!(bv.len(), self.individual_length);
        self.bitmap_data.extend_from_bitslice(bv);
    }

    fn shrink_to_fit(&mut self) {
        self.bitmap_data.shrink_to_fit();
    }

    fn get(&self, bitmap_idx: usize) -> &BitSlice {
        &self.bitmap_data[bitmap_idx*self.individual_length .. (bitmap_idx + 1) * self.individual_length]
    }

    #[allow(dead_code)]
    fn n_sets(&self) -> usize {
        self.bitmap_data.len() / self.individual_length
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        // Serialize using bincode
        bincode::serialize_into(out.by_ref(), &self.bitmap_data).unwrap();
        bincode::serialize_into(out.by_ref(), &self.individual_length).unwrap();
    }

    pub fn load(input: &mut impl std::io::Read) -> Self {
        // Deserialize using bincode
        let bitmap_data: bitvec::vec::BitVec = bincode::deserialize_from(input.by_ref()).unwrap();
        let individual_length: usize = bincode::deserialize_from(input.by_ref()).unwrap();
        assert!(individual_length > 0);
        BitMaps{bitmap_data, individual_length}
    }
}

impl ColexToColorSetMap {

    // sets maps from color set to its index in the distinct color sets
    fn new(sbwt: Arc<SbwtIndex<SubsetMatrix>>, sample_distance: usize, color_bitmap: &bitvec::vec::BitVec, sets: &HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>, n_colors: usize, n_threads: usize) -> Self {
        log::info!("Building mapping from colex to color set id");


        let get_colorset_fn = |colex| &color_bitmap[colex*n_colors..(colex+1)*n_colors];
        let mut sampling_marks = Self::pick_sampled_kmers(sample_distance, &sbwt, None, get_colorset_fn, n_threads);

        let color_set_id_bit_width = sets.len().next_power_of_two().trailing_zeros() as usize;
        let mut sampled_color_set_ids = IntVector::new(color_set_id_bit_width).unwrap(); // In colex order
        sampled_color_set_ids.resize(sampling_marks.count_ones(), 0);
        let mut n_ids_stored = 0_usize;
        for colex in 0..sbwt.n_sets() {
            if sampling_marks.get(colex) {
                let set = &color_bitmap[colex*n_colors .. (colex+1)*n_colors];
                let key = BitKey{bits: set};
                let id = sets[&key]; // Should exist in the hash map. Panics if does not exist.
                sampled_color_set_ids.set(n_ids_stored, id as u64);
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
        self.color_set_ids.serialize(out).unwrap();
    }

    pub fn load(input: &mut impl std::io::Read, sbwt: Arc<SbwtIndex<SubsetMatrix>>) -> Self {
        let sampling = simple_sds_sbwt::bit_vector::BitVector::load(input).unwrap();
        let color_set_ids = IntVector::load(input).unwrap();

        assert_eq!(sampling.len(), sbwt.n_sets());
        assert_eq!(color_set_ids.len(), sampling.count_ones());

        Self{sbwt: sbwt.clone(), sampling, color_set_ids}
    }

    /// Utility function used in construction
    fn pick_sampled_kmers<'a, F: Fn(usize) -> &'a BitSlice + Sync + Send>(sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>, lcs: Option<&LcsArray>, get_colorset_fn: F, n_threads: usize) -> simple_sds_sbwt::bit_vector::BitVector {
        // Find starts of unitigs. Walk forward to the end of the unitig. Segment by color sets.
        
        let marks = simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
        let marks_mutex = Mutex::new(marks); // Need thread-safe modifications
        let marks_mutex_borrow = &marks_mutex; // Passed into the callback

        let callback = |nodes: &[Node], _: &[u8]| {
            let mut marks = marks_mutex_borrow.lock().unwrap();

            let mut prev_set: Option<&BitSlice> = None;
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

impl CompactColexColoring {

    /// Input: 
    /// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
    ///   color j is present in set i.
    pub fn new(sbwt: Arc<SbwtIndex<SubsetMatrix>>, lcs: LcsArray, bm: &bitvec::vec::BitVec, n_colors: usize, sample_distance: usize, n_threads: usize) -> Self {
        let (sets, hashmap) = ColorSets::hash_and_encode_distinct_sets(bm, n_colors);
        let colex_map = ColexToColorSetMap::new(sbwt.clone(), sample_distance, bm, &hashmap, n_colors, n_threads);
        Self {sbwt, lcs, sets, map: colex_map}
    }

    pub fn new_single_colored(sbwt: Arc<SbwtIndex<SubsetMatrix>>, lcs: LcsArray, sample_distance: usize, n_threads: usize) -> Self {
        let n_colors = 1;
        let int_bitwidth = 1;

        let mut dense_sets = BitMaps::new(n_colors);
        let sparse_sets = IntVecs::new(int_bitwidth);

        // Let's make the singleton set dense because a bitvector access is
        // probably a bit cheaper than an intvector access.
        let singleton = bitvec![1];
        dense_sets.push(&singleton); // Singleton set

        let mut is_dense_marks = BitVector::from(RawVector::with_len(1, true));
        is_dense_marks.enable_rank();

        let sets = ColorSets {
            dense_sets, sparse_sets, is_dense_marks, n_colors
        };

        log::info!("Sampling nodes");
        let mut unitig_samples = ColexToColorSetMap::pick_sampled_kmers(sample_distance, &sbwt, Some(&lcs), |_colex| &singleton, n_threads);
        unitig_samples.enable_rank();
        log::info!("Storing color set ids for sampled nodes");
        let color_set_ids = IntVector::with_len(unitig_samples.count_ones(), int_bitwidth, 0).unwrap();
        let colex_map = ColexToColorSetMap{
            sbwt: sbwt.clone(),
            sampling: unitig_samples,
            color_set_ids,
        };
        Self {sbwt, lcs, sets, map: colex_map}
    }

    pub fn colex_to_set_id(&self, colex: usize) -> usize {
        self.map.colex_to_color_set_id(colex)
    }

    pub fn set_id_to_set(&self, id: usize) -> ColorSet<'_>  {
        self.sets.get(id)
    }

    pub fn colex_to_set(&self, colex: usize) -> ColorSet<'_> {
        self.set_id_to_set(self.colex_to_set_id(colex))
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
        let sets = ColorSets::load(input);
        let map = ColexToColorSetMap::load(input, sbwt.clone());
        CompactColexColoring{sbwt, lcs, sets, map}
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

impl ColorSets {

    pub fn get(&self, id: usize) -> ColorSet {
        if self.is_dense_marks.get(id) {
            let set_idx = self.is_dense_marks.rank(id);
            ColorSet::Dense(self.dense_sets.get(set_idx))
        } else {
            let set_idx = self.is_dense_marks.rank_zero(id);
            ColorSet::Sparse(self.sparse_sets.get(set_idx))
        }
    }


    /// Input: 
    /// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
    ///   color j is present in set i.
    /// 
    /// Output:
    /// - Distinct color sets encoded as ColorSets
    /// - HashMap from color set to its index in ColorSets
    pub fn hash_and_encode_distinct_sets(bm: &bitvec::vec::BitVec, n_colors: usize) -> (ColorSets, HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>) {
        assert_eq!(bm.len() % n_colors, 0);
        let n_sets = bm.len() / n_colors;

        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;

        let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

        log::info!("Hashing distinct color sets");

        let mut sparse_sets = IntVecs::new(color_id_bit_width);
        let mut dense_sets = BitMaps::new(n_colors);
        let mut distinct_sets = HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>::default(); // Set -> id
        let bar = indicatif::ProgressBar::new(n_sets as u64);
        for colex in 0..n_sets {
            let set = &bm[colex*n_colors .. (colex+1)*n_colors];
            let key = BitKey{bits: set};
            if !distinct_sets.contains_key(&key) {
                distinct_sets.insert(key, distinct_sets.len());
                if is_dense(set) {
                    dense_sets.push(set);
                    is_dense_marks.push_bit(true);
                } else {
                    sparse_sets.push(set.iter_ones());
                    is_dense_marks.push_bit(false);
                }
            }
            if colex % 100 == 0 {
                bar.inc(100);
            }
        }
        bar.finish();

        log::info!("{} distinct color sets found", distinct_sets.len());

        sparse_sets.shrink_to_fit();
        dense_sets.shrink_to_fit();

        log::info!("{}% of the sets are sparse", sparse_sets.n_sets() as f64 / distinct_sets.len() as f64 * 100.0);

        // Add rank support to dense marks
        log::info!("Building rank support for dense marks");
        let mut is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::from(is_dense_marks);
        is_dense_marks.enable_rank();

        let colorsets = ColorSets {
            is_dense_marks, 
            sparse_sets,
            dense_sets,
            n_colors
        };

        (colorsets, distinct_sets)
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        bincode::serialize_into(out.by_ref(), &self.n_colors).unwrap();
        self.is_dense_marks.serialize(out).unwrap();
        self.sparse_sets.serialize(out);
        self.dense_sets.serialize(out);
    }

    pub fn load(input: &mut impl std::io::Read) -> Self {
        let n_colors: usize = bincode::deserialize_from(input.by_ref()).unwrap();
        let is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::load(input).unwrap();
        let sparse_sets = IntVecs::load(input);
        let dense_sets = BitMaps::load(input);

        assert_eq!(is_dense_marks.len(), sparse_sets.n_sets() + dense_sets.n_sets());
        assert_eq!(n_colors, dense_sets.individual_length);
        assert!(sparse_sets.intvec_data.width() >= n_colors.next_power_of_two().trailing_zeros() as usize);

        ColorSets{is_dense_marks, sparse_sets, dense_sets, n_colors}
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


fn compute_color_id_pairs_and_merged_unitig_sampling(coloring1: &CompactColexColoring, coloring2: &CompactColexColoring, lcs1: &LcsArray, lcs2: &LcsArray, merge_plan: &MergeInterleaving, n_threads: usize) -> (PartitionedReadOnlyIdMap, simple_sds_sbwt::raw_vector::RawVector) {
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();

    // Distinct color id pairs, inserted into key-disjoint hash maps. The values are
    // color set ids within the hash map.
    let hashmaps = vec![HashMap::<(Option::<usize>, Option::<usize>), usize>::new(); n_threads];
    let mut new_id_map = PartitionedIdMap{hashmaps};

    let mut colex1 = 0_usize;
    let mut colex2 = 0_usize;

    let mut color_set_sample_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(merged_len, false);
    let dbg1 = sbwt::dbg::Dbg::new(&(*coloring1.map.sbwt), Some(lcs1), n_threads);
    let dbg2 = sbwt::dbg::Dbg::new(&(*coloring2.map.sbwt), Some(lcs2), n_threads);
    let mut outlabel_buf_1 = Vec::<u8>::new();
    let mut outlabel_buf_2 = Vec::<u8>::new();

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

fn is_dense_set(n_elements: usize, bits_per_color: usize, n_colors: usize) -> bool {
    let intvec_size = n_elements * bits_per_color;
    let bitmap_size = n_colors;
    bitmap_size <= intvec_size
}

fn encode_merged_color_sets(new_id_map: &PartitionedReadOnlyIdMap, coloring1: &CompactColexColoring, coloring2: &CompactColexColoring) -> (IntVecs, BitMaps, simple_sds_sbwt::raw_vector::RawVector){

    let n_colors_1 = coloring1.sets.n_colors;
    let n_colors_2 = coloring2.sets.n_colors;
    let n_colors = n_colors_1 + n_colors_2;
    let bits_per_color = n_colors.next_power_of_two().trailing_zeros() as usize;

    let mut sparse_sets = IntVecs::new(bits_per_color);
    let mut dense_sets = BitMaps::new(n_colors);
    let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

    let id_pairs_in_new_id_order = new_id_map.get_old_ids_sorted_by_new_id();

    // Encode sparse and dense sets
    // TODO: avoid small heap allocations here and instead write directly to the final data structure
    for (_, (left, right)) in id_pairs_in_new_id_order.into_iter() {
        match (left,right) {
            (Some(x), Some(y)) => {
                let set1 = coloring1.set_id_to_set(x);
                let set2 = coloring2.set_id_to_set(y);
                let n_elements = set1.len() + set2.len();

                if is_dense_set(n_elements, bits_per_color, n_colors) {
                    // Dense set -> encode as bitmap
                    let mut concat = bitvec::vec::BitVec::with_capacity(n_colors);
                    concat.extend_from_bitslice(&set1.as_bitvec(n_colors_1));
                    concat.extend_from_bitslice(&set2.as_bitvec(n_colors_2));
                    dense_sets.push(&concat);
                    is_dense_marks.push_bit(true);
                } else {
                    // Sparse set -> encode as integers
                    let mut concat = Vec::<usize>::with_capacity(set1.len() + set2.len());
                    concat.extend(set1.as_intvec());

                    // Offset the colors of the second set by the number of colors in the first
                    concat.extend(set2.as_intvec().iter().map(|x| x + n_colors_1));

                    sparse_sets.push(concat);
                    is_dense_marks.push_bit(false);
                }
            },
            (Some(x), None) => {
                let set1 = coloring1.set_id_to_set(x);
                let n_elements = set1.len();
                if is_dense_set(n_elements, bits_per_color, n_colors) {
                    // Dense
                    let mut concat = bitvec::vec::BitVec::with_capacity(n_colors);
                    concat.extend_from_bitslice(&set1.as_bitvec(n_colors_1));
                    concat.extend_from_bitslice(&bitvec![0; n_colors_2]);

                    dense_sets.push(&concat);
                    is_dense_marks.push_bit(true);
                } else {
                    // Sparse
                    sparse_sets.push(set1.as_intvec());
                    is_dense_marks.push_bit(false);
                }

            },
            (None, Some(y)) => {

                let set2 = coloring2.set_id_to_set(y);
                let n_elements = set2.len();

                if is_dense_set(n_elements, bits_per_color, n_colors) {
                    // Dense
                    let mut concat = bitvec::vec::BitVec::with_capacity(n_colors);
                    concat.extend_from_bitslice(&bitvec![0; n_colors_1]);
                    concat.extend_from_bitslice(&set2.as_bitvec(n_colors_2));

                    dense_sets.push(&concat);
                    is_dense_marks.push_bit(true);
                } else {
                    // Sparse
                    sparse_sets.push(set2.as_intvec().iter().map(|x| x + n_colors_1));
                    is_dense_marks.push_bit(false);
                }
            }
            (None, None) => panic!("Nonexisting color set id pair")
        }
    }

    sparse_sets.shrink_to_fit();
    dense_sets.shrink_to_fit();

    (sparse_sets, dense_sets, is_dense_marks)

}

fn store_new_sampled_color_ids(n_distinct_color_sets: usize, merge_plan: &MergeInterleaving, color_set_sample_marks: &simple_sds_sbwt::bit_vector::BitVector, coloring1: &CompactColexColoring, coloring2: &CompactColexColoring, pair_to_new_id_maps: &PartitionedReadOnlyIdMap) -> simple_sds_sbwt::int_vector::IntVector {
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();

    let bits_per_color_set_id = n_distinct_color_sets.next_power_of_two().trailing_zeros() as usize;
    let mut sampled_ids = simple_sds_sbwt::int_vector::IntVector::with_capacity(color_set_sample_marks.count_ones(), bits_per_color_set_id).unwrap();
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
            sampled_ids.push(id as u64);
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s2[merged_colex] as usize;
    }

    sampled_ids
}

pub fn merge_compact_colorings(coloring1: CompactColexColoring, coloring2: CompactColexColoring, optimize_peak_ram: bool, n_threads: usize) -> CompactColexColoring {

    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::MergeInterleaving::new(&(*coloring1.map.sbwt), &(*coloring2.map.sbwt), optimize_peak_ram, n_threads);

    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();    

    let n_colors_1 = coloring1.sets.n_colors;
    let n_colors_2 = coloring2.sets.n_colors;
    let n_colors = n_colors_1 + n_colors_2;

    log::info!("Computing color id pairs and merged sampling");
    let (new_id_map, color_set_sample_marks) = compute_color_id_pairs_and_merged_unitig_sampling(&coloring1, &coloring2, &coloring1.lcs, &coloring2.lcs, &merge_plan, n_threads);

    let mut color_set_sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(color_set_sample_marks);
    color_set_sample_marks.enable_rank();
    let n_sampled = color_set_sample_marks.rank(color_set_sample_marks.len());
    log::info!("Sampled {} out of {} SBWT nodes ({:.2}%)", n_sampled, merged_len, n_sampled as f64 / merged_len as f64 * 100.0);

    log::info!("Encoding distinct merged color sets");
    let (sparse_sets, dense_sets, is_dense_marks) = encode_merged_color_sets(&new_id_map, &coloring1, &coloring2);

    log::info!("{}% of the sets are sparse", sparse_sets.n_sets() as f64 / (sparse_sets.n_sets() + dense_sets.n_sets()) as f64 * 100.0);

    log::info!("Storing new sampled color set ids");
    let n_distinct_color_sets = new_id_map.total_len(); 
    let sampled_ids = store_new_sampled_color_ids(n_distinct_color_sets, &merge_plan, &color_set_sample_marks, &coloring1, &coloring2, &new_id_map);

    // Add rank support to dense marks
    log::info!("Building rank support for dense marks");
    let mut is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::from(is_dense_marks);
    is_dense_marks.enable_rank();

    let colorsets = ColorSets {
        is_dense_marks, 
        sparse_sets,
        dense_sets,
        n_colors
    };

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

    let new_coloring = CompactColexColoring { 
        sbwt: merged_sbwt.clone(),
        lcs: merged_lcs,
        sets: colorsets, 
        map: ColexToColorSetMap {
            sbwt: merged_sbwt.clone(), 
            sampling: color_set_sample_marks, 
            color_set_ids: sampled_ids 
        }
    };

    log::info!("Color merge finished");
    new_coloring

}


#[cfg(test)]
mod tests {
    use std::path::Path;

    use jseqio::seq_db::SeqDB;
    use simple_sds_sbwt::ops::BitVec;

    use crate::colored_kmers::ColoredKmers;

    use super::merge_compact_colorings;


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

    #[test]
    fn test_merge() {

        if std::env::var("RUST_LOG").is_err() {
            std::env::set_var("RUST_LOG", "info")
        }
        env_logger::init();

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

            let cc1 = ColoredKmers::new_from_seq_dbs::<&Path>(dbs1, k, 3, None);
            let cc2 = ColoredKmers::new_from_seq_dbs::<&Path>(dbs2, k, 3, None);
            let mut cc_both = ColoredKmers::new_from_seq_dbs::<&Path>(dbs_both, k, 3, None);

            cc_both.build_sbwt_select_support();

            let sample_distance = 3;
            let ccc1 = cc1.compress_colors(sample_distance, 3);
            let ccc2 = cc2.compress_colors(sample_distance, 3);

            let ccc_merged = merge_compact_colorings(ccc1, ccc2, true, 3);
            let sbwt_merged = &ccc_merged.sbwt;

            for colex in 0..cc_both.sbwt().n_sets() {
                let kmer = cc_both.sbwt().access_kmer(colex);
                let true_colors = cc_both.get_color_set(&kmer);

                if kmer.iter().all(|c| *c != b'$') { // Not a dummy k-mer
                    let range = sbwt_merged.search(&kmer).unwrap();
                    assert_eq!(range.len(), 1);
                    let colex_merged = range.start;
                    let merged_colors = ccc_merged.colex_to_set(colex_merged).as_bitvec(cc_both.n_colors());

                    eprintln!("{} {} {} {:?} {} {}", colex, String::from_utf8_lossy(&kmer), true_colors, sbwt_merged.search(&kmer), ccc_merged.map.sampling.get(colex_merged), ccc_merged.colex_to_set_id(colex_merged));
                    assert_eq!(true_colors, merged_colors);
                }

            }
        }
    }
}

//
//
//
//
//
// NEW CODE HERE BELOW
//
//
//
//
//

pub struct ColorSetViewIterator<'a> {
    set: ColorSet<'a>,
    pos: usize, // Interpreted differently depending of whether this is Sparse or Dense
}

impl<'a> Iterator for ColorSetViewIterator<'a> {
    type Item = usize;

    #[allow(clippy::bool_comparison)]
    fn next(&mut self) -> Option<Self::Item> {
        match &self.set {
            ColorSet::Dense(bit_slice) => {
                // Rewind to the next 1-bit (todo: word parallelism)
                while self.pos < bit_slice.len() && bit_slice[self.pos] == false {
                    self.pos += 1;
                }

                if self.pos < bit_slice.len() {
                    Some(self.pos)
                } else {
                    None
                }
            },
            ColorSet::Sparse(int_vec_slice) => {
                if self.pos == int_vec_slice.end - int_vec_slice.start {
                    None
                } else {
                    let x = int_vec_slice.vec.get(int_vec_slice.start + self.pos);
                    self.pos += 1;
                    Some(x as usize)
                }
            },
        }
    }
}

impl<'a> coloring_interface::ColorSetView<'a> for ColorSet<'a> {
    type Iter = ColorSetViewIterator<'a>;

    fn iter(&self) -> Self::Iter {
        ColorSetViewIterator{
            set: *self,
            pos: 0,
        }
    }
}

impl coloring_interface::ColorSetStorage for CompactColexColoring {
    type SetView<'a> = ColorSet<'a>; 
    type OwnedSet = Vec<usize>; // TODO

    fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow> {
        self.set_id_to_set(id)
    }

    fn get_owned_set(&self, id: usize) -> Self::OwnedSet {
        self.set_id_to_set(id).iter().collect()
    }

    fn get_empty_set(&self) -> Self::OwnedSet {
        vec![]
    }

    fn get_full_set(&self) -> Self::OwnedSet {
        todo!()
    }
}