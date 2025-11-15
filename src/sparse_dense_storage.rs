use simple_sds_sbwt::int_vector::IntVector;
use simple_sds_sbwt::serialize::Serialize;
use simple_sds_sbwt::{ops::{Access, BitVec, Push, Rank, Resize, Vector}, raw_vector::{AccessRaw, PushRaw}};
use bitvec::order::Lsb0;
use bitvec::{field::BitField, slice::BitSlice};
use bitvec::bitvec;

/// A data structure for storing arbitary set of sets of integers, such that dense
/// sets are encoded as bitmaps, and sparse sets as lists of integers.
pub struct SparseDenseStorage{
    dense_sets: BitMaps,
    sparse_sets: IntVecs,
    n_colors: usize,
    is_dense_marks: simple_sds_sbwt::bit_vector::BitVector, // Has rank support.
}


impl SparseDenseStorage {

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

