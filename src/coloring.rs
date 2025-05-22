use bitvec::slice::BitSlice;
use sbwt::{SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec, Push, Vector}};

struct IntVecSlice<'a> {
    vec: &'a IntVector, // Todo: do not store this
    start: usize,
    end: usize, // Exclusive end
}

// TODO: should not do this: BitSlice store a pointer for each. We only need an offset and a length.
// Actually since the bitslices are of the same known width, we only need an offset.
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


struct ColorSets<'a, 'b> {
    sets: Vec<ColorSet<'a>>, // Lifetime 'a points to bitmap_data and intvec_data
    bitmap_data: bitvec::vec::BitVec, // Concatenation of dense sets as bitmaps
    intvec_data: IntVector, // Concatenation of sparse sets as int vecs
    sbwt: &'b SbwtIndex<SubsetMatrix>, // Lifetime 'b can outlive this struct
    sampling: simple_sds_sbwt::bit_vector::BitVector // Marks colex ranks that have a color set stored. Has rank support.
}

impl ColorSets<'_, '_> {
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