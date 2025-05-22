use bitvec::slice::BitSlice;
use sbwt::{SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Pack, Push, Resize, Vector}};

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



struct ColorSets<'a> {
    //sets: Vec<ColorSet<'a>>, // Lifetime 'a points to bitmap_data and intvec_data
    bitmap_data: bitvec::vec::BitVec, // Concatenation of dense sets as bitmaps
    intvec_data: IntVector, // Concatenation of sparse sets as int vecs
    is_dense_marks: simple_sds_sbwt::bit_vector::BitVector, // Has rank support.
    sbwt: &'a SbwtIndex<SubsetMatrix>, // Lifetime 'b can outlive this struct
    sampling: simple_sds_sbwt::bit_vector::BitVector // Marks colex ranks that have a color set stored. Has rank support.
}

impl ColorSets<'_> {
    fn get(&self, colex: usize) -> &ColorSet {
        if self.sampling.get(colex) {
            &self.sets[colex]
        } else {
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
    fn new_from_bitmaps(&self, bm: &bitvec::vec::BitVec, n_colors: usize, sample_distance: usize, sbwt: &SbwtIndex<SubsetMatrix>) -> Self {
        assert_eq!(bm.len() % n_colors, 0);
        let sbwt_len = bm.len() % n_colors;
        assert_eq!(sbwt_len, sbwt.n_sets());

        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;

        let mut intvec_data = IntVector::new(color_id_bit_width).unwrap();
        let mut bitvec_data = bitvec::vec::BitVec::new();

        let mut distinct_sets = std::collections::HashSet::<&BitSlice>::new();
        let mut distinct_sets_encoded = Vec::<ColorSet>::new();
        for colex in 0..sbwt_len {
            let set = &bm[colex*n_colors .. colex*(n_colors+1)];
            if !distinct_sets.contains(set) {
                distinct_sets.insert(set);
                if is_dense(set) {
                    bitvec_data.extend_from_bitslice(set);
                    let new_bits = &bitvec_data[bitvec_data.len() - n_colors .. bitvec_data.len()];
                    distinct_sets_encoded.push(ColorSet::Dense(new_bits));
                } else {
                    let old_end = intvec_data.len();
                    intvec_data.extend(set.iter_ones());
                    distinct_sets_encoded.push
                }
            }
        }

        drop(distinct_sets); // Save memory

        // Encode distinct sets
        for set_idx in 0..n_distinct {
            let slice = distinct_sets[set_idx];
        }

        todo!();
    }
}