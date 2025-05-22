use sbwt::{SbwtIndex, SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::{int_vector::IntVector, ops::{Access, BitVec}};

enum ColorSet {
    Dense(bitvec::vec::BitVec),
    Sparse(IntVector),
}

impl ColorSet {
    fn push_colors(&self, buf: &mut Vec<usize>) {
        match self {
            ColorSet::Dense(bv) => {
                for i in bv.iter_ones() {
                    buf.push(i);
                }
            },
            ColorSet::Sparse(iv) => {
                for x in iv.iter() {
                    buf.push(x as usize);
                }
            },
        }
    }
}

struct ColorSets<'a> {
    sets: Vec<ColorSet>,
    sbwt: &'a SbwtIndex<SubsetMatrix>,
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
}