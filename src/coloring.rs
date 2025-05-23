use bitvec::{field::BitField, order::Lsb0, slice::BitSlice};
use sbwt::{dbg::{Dbg, Node}, SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Pack, Push, Rank, Resize, Vector}};
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;
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

    fn push(&mut self, vec: &IntVector) {
        assert_eq!(vec.width(), self.intvec_data.width());
        self.intvec_data.extend(vec.iter());
        self.ends.push(self.intvec_data.len());
    }

    fn shrink_to_fit(&mut self) {
        self.intvec_data.resize(self.intvec_data.len(), 0);
    }

    fn get(&self, vec_idx: usize) -> IntVecSlice {
        IntVecSlice{vec: &self.intvec_data, start: self.ends[vec_idx], end: self.ends[vec_idx+1]}
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

    fn push(&mut self, bv: bitvec::vec::BitVec) {
        assert_eq!(bv.len(), self.individual_length);
        self.bitmap_data.extend_from_bitslice(&bv);
    }

    fn shrink_to_fit(&mut self) {
        self.bitmap_data.shrink_to_fit();
    }

    fn get(&self, bitmap_idx: usize) -> &BitSlice {
        &self.bitmap_data[bitmap_idx*self.individual_length .. (bitmap_idx + 1) * self.individual_length]
    }
}



pub struct ColorSets<'a> {
    //sets: Vec<ColorSet<'a>>, // Lifetime 'a points to bitmap_data and intvec_data
    bitmaps: BitMaps,// Concatenation of dense sets as bitmaps
    intvecs: IntVecs, // Concatenation of sparse sets as int vecs
    is_dense_marks: simple_sds_sbwt::bit_vector::BitVector, // Has rank support.
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    sampling: simple_sds_sbwt::bit_vector::BitVector // Marks colex ranks that have a color set stored. Has rank support.
}

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

fn pick_sampled_kmers(n_colors: usize, sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> bitvec::vec::BitVec {
    // Find starts of unitigs. Walk forward to the end of the unitig. Segment by color sets.
    
    // TODO: for now, just mark every non-dummy node.
    log::info!("WARNING: unitig sampling not implement, marking all nodes instead");

    let dbg = sbwt::dbg::Dbg::new(&sbwt, None, 1); // Todo: n_threads

    let mut marks = bitvec::vec::BitVec::new();
    marks.resize(sbwt.n_sets(), false);
    for node in dbg.node_iterator() {
        marks.set(node.id, true);
    }

    log::info!("Unitig sampling finished");

    marks
}

#[derive(Debug, Eq, PartialEq)]
pub struct BitKey<'a> {
    pub bits: &'a BitSlice,
}

impl<'a> Hash for BitKey<'a> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let len = self.bits.len();
        assert!(
            len <= usize::BITS as usize,
            "BitSlice too long to load into usize"
        );

        let word: usize = self.bits.load();
        word.hash(state); // hash as an integer
        len.hash(state);  // include length to distinguish e.g. 0b1 from 0b10
    }
}


impl ColorSets<'_> {
    pub fn get(&self, colex: usize) -> ColorSet {
        if self.sampling.get(colex) {
            // This set is stored
            if self.is_dense_marks.get(colex) {
                let set_idx = self.is_dense_marks.rank(colex);
                return ColorSet::Dense(&self.bitmaps.get(set_idx));
            } else {
                let set_idx = self.is_dense_marks.rank_zero(colex);
                return ColorSet::Sparse(self.intvecs.get(set_idx));
            }
        } else {
            // This set is not stored -> walk forward in the de Bruijn graph
            for char_idx in 0..self.sbwt.alphabet().len() {
                if self.sbwt.sbwt().set_contains(colex, char_idx as u8) {
                    let new_colex = self.sbwt.lf_step(colex, char_idx);
                    return self.get(new_colex);
                }
            }
            panic!("Bug in color set sampling: dead end in SBWT graph");
        }
    }



    /// Input: 
    /// - Color sets in bitmap representation: bm[i * n_colors + j] tells whether
    ///   color j is present in set i.
    /// - sample_distance: max walk length to the next sampled color set in a unitig 
    /// Sbwt needs to have select support!
    pub fn new_from_bitmaps(bm: &bitvec::vec::BitVec, n_colors: usize, sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> Self {
        assert_eq!(bm.len() % n_colors, 0);
        let sbwt_len = bm.len() / n_colors;
        assert_eq!(sbwt_len, sbwt.n_sets());

        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;

        let mut is_dense_marks = bitvec::vec::BitVec::<usize, Lsb0>::new();
        is_dense_marks.resize(sbwt_len, false);

        let sampling_marks = pick_sampled_kmers(n_colors, sample_distance, sbwt);
        log::info!("Hashing distinct color sets");

        let mut intvec_data = IntVector::new(color_id_bit_width).unwrap();
        let mut bitvec_data = bitvec::vec::BitVec::<usize, Lsb0>::new();

        let mut intvec_data_ends = vec![0_usize];

        let mut distinct_sets = std::collections::HashMap::<BitKey, usize, BuildHasherDefault::<FxHasher>>::default(); // Set -> id
        let bar = indicatif::ProgressBar::new(sbwt_len as u64);
        for colex in 0..sbwt_len {
            let set = &bm[colex*n_colors .. colex*(n_colors+1)];
            let key = BitKey{bits: set};
            if !distinct_sets.contains_key(&key) {
                distinct_sets.insert(key, distinct_sets.len());
                if is_dense(set) {
                    bitvec_data.extend_from_bitslice(set);
                } else {
                    intvec_data.extend(set.iter_ones());
                    intvec_data_ends.push(intvec_data.len());
                }
            }
            if colex % 100 == 0 {
                bar.inc(100);
            }
        }
        bar.finish();

        log::info!("{} distinct color sets found", distinct_sets.len());
        log::info!("{} of the sets are sparse ({}%)", intvec_data_ends.len() - 1, (intvec_data_ends.len() - 1) as f64 / distinct_sets.len() as f64 * 100.0 );

        log::info!("Storing color set pointers for sampled nodes");

        // Store color set pointers only for the sampled nodes
        let color_set_id_bit_width = distinct_sets.len().next_power_of_two().trailing_zeros() as usize;
        let mut sampled_color_set_ids = IntVector::new(color_set_id_bit_width).unwrap(); // In colex order
        sampled_color_set_ids.resize(sampling_marks.count_ones(), 0);
        let mut n_ids_stored = 0_usize;
        for colex in 0..sbwt_len {
            if sampling_marks[colex] {
                let set = &bm[colex*n_colors .. colex*(n_colors+1)];
                let key = BitKey{bits: set};
                let id = distinct_sets[&key]; // Should exist in the hash map. Panics if does not exist.
                sampled_color_set_ids.set(n_ids_stored, id as u64);
                n_ids_stored += 1;
            }
        }

        log::info!("Color compression finished");
        todo!();
    }
}