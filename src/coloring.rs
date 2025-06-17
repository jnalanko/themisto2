use bitvec::{field::BitField, slice::BitSlice};
use bitvec::bitvec;
use rand::seq::index::sample;
use sbwt::{dbg::Node, SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Push, Rank, Resize, Vector}, raw_vector::{AccessRaw, PushRaw}};
use rustc_hash::FxHasher;
use std::cmp::max;
use std::{cmp::min, collections::HashMap, hash::BuildHasherDefault, sync::Mutex};
use std::hash::{Hash, Hasher};

/// This is the main data structure in this file: a set of compressed color sets, and a mapping
/// from SBWT colex ranks to color sets such that we can look up the color set of a k-mer by its
/// colex rank in the SBWT. 
pub struct CompactColexColoring<'a> {
    sets: ColorSets, // Distinct color sets
    map: ColexToColorSetMap<'a>, // A mapping from the colex rank of a k-mer in the SBWT into a color set id in `sets`
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
pub struct ColexToColorSetMap<'a> {
    sbwt: &'a SbwtIndex<SubsetMatrix>,
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

    pub fn as_bitvec(&self) -> bitvec::vec::BitVec {
        match self {
            ColorSet::Dense(bv) => {
                (*bv).into()
            },
            ColorSet::Sparse(iv) => {
                let mut bv = bitvec![0; self.len()];
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
}

impl<'a> ColexToColorSetMap<'a> {

    // sets maps from color set to its index in the distinct color sets
    fn new(sbwt: &'a SbwtIndex<SubsetMatrix>, sample_distance: usize, color_bitmap: &bitvec::vec::BitVec, sets: &HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>, n_colors: usize, n_threads: usize) -> Self {
        log::info!("Building mapping from colex to color set id");

        let mut sampling_marks = Self::pick_sampled_kmers(n_colors, sample_distance, sbwt, color_bitmap, n_threads);

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

    fn colex_to_color_set_id(&self, colex: usize) -> usize {
        if self.sampling.get(colex) {
            // This set is stored
            self.color_set_ids.get(self.sampling.rank(colex)) as usize
        } else {
            // This set is not stored -> walk forward in the de Bruijn graph
            for char_idx in 0..self.sbwt.alphabet().len() {
                if self.sbwt.sbwt().set_contains(colex, char_idx as u8) {
                    let new_colex = self.sbwt.lf_step(colex, char_idx);
                    return self.colex_to_color_set_id(new_colex); // Todo: no recursion
                }
            }
            panic!("Bug in color set sampling: dead end in SBWT graph");
        }
    }

    fn serialize(&self, out: &mut impl std::io::Write) {
        todo!();
    }

    fn load(&self, input: &mut impl std::io::Read, sbwt: &SbwtIndex<SubsetMatrix>) -> Self {
        todo!();
    }

    /// Utility function used in construction
    fn pick_sampled_kmers(n_colors: usize, sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>, bitmaps: &BitSlice, n_threads: usize) -> simple_sds_sbwt::bit_vector::BitVector {
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
                let cur_set = &bitmaps[colex*n_colors..(colex+1)*n_colors];

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
        let dbg = sbwt::dbg::Dbg::new(sbwt, None, n_threads);

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

impl<'a> CompactColexColoring<'a> {

    /// Input: 
    /// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
    ///   color j is present in set i.
    pub fn new(sbwt: &'a SbwtIndex<SubsetMatrix>, bm: &bitvec::vec::BitVec, n_colors: usize, sample_distance: usize, n_threads: usize) -> Self {
        let (sets, hashmap) = ColorSets::hash_and_encode_distinct_sets(bm, n_colors);
        let colex_map = ColexToColorSetMap::new(sbwt, sample_distance, bm, &hashmap, n_colors, n_threads);

        Self {sets, map: colex_map}
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

    pub fn merge<'b>(left: CompactColexColoring<'a>, right: CompactColexColoring<'a>) -> (CompactColexColoring<'b>, SbwtIndex<SubsetMatrix>) {
        todo!();
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
}

fn merge_colorings(coloring1: CompactColexColoring, coloring2: CompactColexColoring, optimize_peak_ram: bool, n_threads: usize) -> CompactColexColoring {
    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::merge::MergeInterleaving::new(coloring1.map.sbwt, coloring2.map.sbwt, optimize_peak_ram, n_threads);

    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());
    let merged_len = merge_plan.s1.len();    

    let n_colors_1 = coloring1.sets.n_colors;
    let n_colors_2 = coloring2.sets.n_colors;
    let n_colors = n_colors_1 + n_colors_2;
    let bits_per_color = n_colors.next_power_of_two().trailing_zeros() as usize;

    log::info!("Hashing distinct color set id pairs");
    let mut distinct_ids = std::collections::HashSet::<(Option<usize>, Option<usize>)>::new();
    let mut colex1 = 0_usize;
    let mut colex2 = 0_usize;

    let mut color_set_sample_marks = simple_sds_sbwt::raw_vector::RawVector::with_len(merged_len, false);
    let dbg1 = sbwt::dbg::Dbg::new(coloring1.map.sbwt, None, n_threads);
    let dbg2 = sbwt::dbg::Dbg::new(coloring2.map.sbwt, None, n_threads);
    let mut outlabel_buf_1 = Vec::<u8>::new();
    let mut outlabel_buf_2 = Vec::<u8>::new();

    for merged_colex in 0..merged_len {
        let color_set_id_1 = if merge_plan.s1[merged_colex] &&  coloring1.map.sampling.get(colex1){
            Some(coloring1.colex_to_set_id(colex1))
        } else {
            None
        };

        let color_set_id_2 = if merge_plan.s2[merged_colex] &&  coloring2.map.sampling.get(colex2){
            Some(coloring2.colex_to_set_id(colex2))
        } else {
            None
        };

        if color_set_id_1.is_some() || color_set_id_2.is_some() {
            distinct_ids.insert((color_set_id_1, color_set_id_2));
            color_set_sample_marks.set_bit(merged_colex, true);
        } else if merge_plan.s1[merged_colex] && merge_plan.s2[merged_colex] {
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
            match (outlabel_buf_1.first(), outlabel_buf_2.first()) {
                (Some(a), Some(b)) => {
                    if a != b { // Case 2 in the comment above
                        color_set_sample_marks.set_bit(merged_colex, true);
                    } // The else-branch would be case 1 but then there is nothing to do
                }
                _ => panic!("Bug at computing color set samples bit vector in merge") // Both should have outdegree > 0
            }
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s1[merged_colex] as usize;
    }

    let color_set_sample_marks = simple_sds_sbwt::bit_vector::BitVector::from(color_set_sample_marks);
    let n_sampled = color_set_sample_marks.count_ones();
    log::info!("Sampled {} out of {} SBWT nodes ({:.2}%)", n_sampled, merged_len, n_sampled as f64 / merged_len as f64 * 100.0);

    log::info!("Constructing distinct merged color sets");
    let mut id_pairs: Vec<(Option<usize>, Option<usize>)> = distinct_ids.into_iter().collect();
    id_pairs.sort_unstable();
    let mut pair_to_new_id = HashMap::<(Option<usize>, Option<usize>), usize>::new();

    let mut sparse_sets = IntVecs::new(bits_per_color);
    let mut dense_sets = BitMaps::new(n_colors);
    let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

    for (new_id, (left, right)) in id_pairs.into_iter().enumerate() {
        pair_to_new_id.insert((left,right), new_id);
        match (left,right) {
            (Some(x), Some(y)) => {
                let set1 = coloring1.set_id_to_set(x);
                let set2 = coloring2.set_id_to_set(y);
                let n_elements = set1.len() + set2.len();

                if n_elements * bits_per_color > n_colors {
                    // Dense set -> encode as bitmap
                    let mut concat = bitvec::vec::BitVec::with_capacity(n_colors);
                    concat.extend_from_bitslice(&set1.as_bitvec());
                    concat.extend_from_bitslice(&set2.as_bitvec());
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
                if n_elements * bits_per_color > n_colors {
                    // Dense
                    let mut concat = bitvec::vec::BitVec::with_capacity(n_colors);
                    concat.extend_from_bitslice(&set1.as_bitvec());
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
                if n_elements * bits_per_color > n_colors {
                    // Dense
                    let mut concat = bitvec::vec::BitVec::with_capacity(n_colors);
                    concat.extend_from_bitslice(&bitvec![0; n_colors_1]);
                    concat.extend_from_bitslice(&set2.as_bitvec());

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

    log::info!("{}% of the sets are sparse", sparse_sets.n_sets() as f64 / (sparse_sets.n_sets() + dense_sets.n_sets()) as f64 * 100.0);

    log::info!("Storing new sampled color set ids");
    let mut sampled_ids = simple_sds_sbwt::int_vector::IntVector::with_capacity(color_set_sample_marks.count_ones(), bits_per_color).unwrap();
    colex1 = 0_usize;
    colex2 = 0_usize;
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
            let id = pair_to_new_id[&(color_set_id_1, color_set_id_2)];
            sampled_ids.push(id as u64);
        }

        colex1 += merge_plan.s1[merged_colex] as usize;
        colex2 += merge_plan.s1[merged_colex] as usize;
    }

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
    let sbwt1 = coloring1.map.sbwt.clone(); // Todo: avoid clone. Currently unavoidable because we have just a reference to the SBWT. 
    drop(coloring1);

    let sbwt2 = coloring2.map.sbwt.clone(); // Todo: avoid clone
    drop(coloring2);

    let merged_sbwt = SbwtIndex::merge(sbwt1, sbwt2, merge_plan, precalc_len, n_threads);
    let new_color_set_ids = IntVector::new(64).unwrap(); // TODO

    CompactColexColoring { sets: colorsets, map: ColexToColorSetMap{sbwt: &merged_sbwt, sampling: color_set_sample_marks, color_set_ids: sampled_ids} }

}