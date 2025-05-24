use bitvec::{field::BitField, order::Lsb0, slice::BitSlice};
use clap::error::KindFormatter;
use sbwt::{dbg::{Dbg, Node}, SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Pack, Push, Rank, Resize, Vector}, raw_vector::{AccessRaw, PushRaw}};
use rustc_hash::FxHasher;
use std::{cmp::min, collections::HashMap, hash::BuildHasherDefault};
use std::hash::{Hash, Hasher};

// This enum is only for passing references to individual sets around.
enum ColorSet<'a> {
    Dense(&'a BitSlice),
    Sparse(IntVecSlice<'a>),
}

impl ColorSet<'_> {
    fn push_colors(&self, buf: &mut Vec<usize>) {
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
}

pub fn pick_sampled_kmers(n_colors: usize, sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>, sets: &HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>) -> simple_sds_sbwt::bit_vector::BitVector {
    // Find starts of unitigs. Walk forward to the end of the unitig. Segment by color sets.
    
    // TODO: for now, just mark every non-dummy node.
    log::info!("WARNING: unitig sampling not implement, marking all nodes instead");

    let dbg = sbwt::dbg::Dbg::new(&sbwt, None, 1); // Todo: n_threads

    let mut marks = simple_sds_sbwt::raw_vector::RawVector::with_len(sbwt.n_sets(), false);
    for node in dbg.node_iterator() {
        marks.set_bit(node.id, true);
    }

    let marks = simple_sds_sbwt::bit_vector::BitVector::from(marks);
    log::info!("Unitig sampling finished");

    marks
}

fn is_dense(bv: &BitSlice) -> bool {
    let n_colors = bv.len();
    let n_elements = bv.count_ones();
    let bits_per_color = n_colors.next_power_of_two().trailing_zeros() as usize;
    let bitmap_size = n_colors;
    let intvec_size = n_elements * bits_per_color;

    bitmap_size <= intvec_size
}

/*
    if bitmap_size <= intvec_size {
        
        ColorSet::Dense(&bv)
    } else {
        let mut iv = IntVector::with_capacity(n_elements, bits_per_color).unwrap();
        for i in bv.iter_ones() {
            iv.push(i as u64);
        }
        ColorSet::Sparse(&iv)
    }

*/

struct IntVecs {
    intvec_data: IntVector, // Concatenation of IntVecs

    // Ends of individual intvecs, such that ends[0] = 0 and ends[i+1] is the
    // exclusive end of the i-th vector.
    ends: Vec<usize>, 
}

struct IntVecSlice<'a> {
    vec: &'a IntVector,
    start: usize,
    end: usize, // Exclusive end
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

struct BitMaps {
    bitmap_data: bitvec::vec::BitVec, // Concatenation of bit vectors
    individual_length: usize, // Length of each bitmap in bitmap_data
}

impl BitMaps {
    fn new(individual_length: usize) -> Self {
        BitMaps{bitmap_data: bitvec::vec::BitVec::new(), individual_length}
    }

    fn push(&mut self, bv: &bitvec::slice::BitSlice) {
        assert_eq!(bv.len(), self.individual_length);
        self.bitmap_data.extend_from_bitslice(&bv);
    }

    fn shrink_to_fit(&mut self) {
        self.bitmap_data.shrink_to_fit();
    }

    fn get(&self, bitmap_idx: usize) -> &BitSlice {
        &self.bitmap_data[bitmap_idx*self.individual_length .. (bitmap_idx + 1) * self.individual_length]
    }

    fn n_sets(&self) -> usize {
        self.bitmap_data.len() / self.individual_length
    }
}

pub struct ColexToColorSetMap<'a> {
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    sampling: simple_sds_sbwt::bit_vector::BitVector, // Marks colex ranks that have a color set stored. Has rank support.
    color_set_ids: IntVector, // One color set id for every 1-bit in the sampling
}

impl<'a> ColexToColorSetMap<'a> {

    // sets maps from color set to its index in the distinct color sets
    fn new(sbwt: &'a SbwtIndex<SubsetMatrix>, sample_distance: usize, color_bitmap: &bitvec::vec::BitVec, sets: &HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>, n_colors: usize) -> Self {
        log::info!("Building mapping from colex to color set id");

        let mut sampling_marks = pick_sampled_kmers(n_colors, sample_distance, sbwt, sets);

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
            self.sampling.rank(colex)
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

    fn sbwt_len(&self) -> usize {
        self.sbwt.n_sets()
    }

    fn serialize(&self, out: &mut impl std::io::Write) {
        todo!();
    }

    fn load(&self, input: &mut impl std::io::Read, sbwt: &SbwtIndex<SubsetMatrix>) -> Self {
        todo!();
    }
}

pub struct ColorSets {
    dense_sets: BitMaps,
    sparse_sets: IntVecs,
    is_dense_marks: simple_sds_sbwt::bit_vector::BitVector, // Has rank support.
}

pub struct CompactColexColoring<'a> {
    sets: ColorSets, // Distinct color sets
    map: ColexToColorSetMap<'a>,
}

impl<'a> CompactColexColoring<'a> {

    /// Input: 
    /// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
    ///   color j is present in set i.
    pub fn new(sbwt: &'a SbwtIndex<SubsetMatrix>, bm: &bitvec::vec::BitVec, n_colors: usize, sample_distance: usize) -> Self {
        let (sets, hashmap) = hash_and_encode_distinct_sets(bm, n_colors);
        let colex_map = ColexToColorSetMap::new(sbwt, sample_distance, bm, &hashmap, n_colors);

        Self {sets, map: colex_map}
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
}

/// Input: 
/// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
///   color j is present in set i.
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
                dense_sets.push(&set);
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
        dense_sets
    };

    (colorsets, distinct_sets)
}

/*
fn is_first_kmer_of_unitig(dbg: &Dbg<SubsetMatrix>, v: Node) -> bool {
    if dbg.indegree(v) > 1 {
        return true;
    }
    if let Some(u) = dbg.follow_inedge(v, 0) {
        dbg.outdegree(u) > 1
    } else {
        true
    }
}

// Returns the sequence of nodes and the label of the unitig
// The out_labels_buf is working space for the function. Don't assume
// anything about its contents when the function returns.
fn walk_unitig_from(dbg: &Dbg<SubsetMatrix>, mut v: Node, out_labels_buf: &mut Vec<u8>) -> (Vec<Node>, Vec<u8>){
    let v0 = v;
    let mut nodes = Vec::<Node>::new();
    nodes.push(v);
    
    let mut label = Vec::<u8>::new();
    dbg.push_node_kmer(v, &mut label); 

    while dbg.outdegree(v) == 1 {
        out_labels_buf.clear();
        dbg.push_outlabels(v, out_labels_buf);
        let c = out_labels_buf[0];
        v = dbg.follow_outedge(v, c).unwrap();
        if v != v0 && dbg.indegree(v) == 1 {
            label.push(c);
            nodes.push(v);
        } else { break; }
    }

    (nodes, label)
}
*/