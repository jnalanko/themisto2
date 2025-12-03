use std::{cmp::min, io::{Cursor, Read}, ops::Range};

use bitvec::{array::BitArray, order::Lsb0};
use simple_sds_sbwt::serialize::Serialize;

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
pub const IS_DNA: BitArray<[u32; 8]> = bitvec::bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

pub(crate) fn bitvec_to_simple_sds_raw_bitvec(mut bv: bitvec::vec::BitVec) -> simple_sds_sbwt::raw_vector::RawVector {
    // TODO: We really hope that usize equals u64 here, otherwise this this is probably broken.
    // Let's use the deserialization function in simple_sds_sbwt for a raw bitvector.
    // It requires the following header:
    let mut header = [0u64, 0u64]; // bits, words
    header[0] = bv.len() as u64; // Assumes little-endian byte order
    header[1] = bv.len().div_ceil(64) as u64;

    let header_bytes = bytemuck::cast_slice(&header);

    // Make sure the leftover bits in the last word are zeros. Simple-sds
    // depends on this, but the bitvec crate does not guarantee this!
    // The undefined padding bytes have broken my code before, so this is
    // crucial.
    let original_len = bv.len();
    bv.resize(original_len.next_multiple_of(64), false);
    bv.resize(original_len, false);

    let raw_data = bytemuck::cast_slice(bv.as_raw_slice());
    let mut data_with_header = Cursor::new(header_bytes).chain(Cursor::new(raw_data));

    simple_sds_sbwt::raw_vector::RawVector::load(&mut data_with_header).unwrap()
}

pub(crate) fn bitvec_to_simple_sds_bitvec(bv: bitvec::vec::BitVec) -> simple_sds_sbwt::bit_vector::BitVector {
    simple_sds_sbwt::bit_vector::BitVector::from(bitvec_to_simple_sds_raw_bitvec(bv))
}

#[allow(dead_code)]
pub struct VecVecSeqStream{
    vv: Vec<Vec<u8>>,
    pos: usize,
}

impl sbwt::SeqStream for VecVecSeqStream {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.pos == self.vv.len() {
            None
        } else {
            self.pos += 1;
            Some(&self.vv[self.pos-1])
        }
    }
}

impl VecVecSeqStream {
    #[allow(dead_code)]
    pub fn new(seqs: Vec<Vec<u8>>) -> Self {
        Self { vv: seqs, pos: 0 }
    }
}

pub fn segment_range(range: Range<usize>, n_pieces: usize) -> Vec<Range<usize>> {
    let segment_len = range.len().div_ceil(n_pieces);
    let mut pieces: Vec<Range<usize>> = vec![];
    for t in 0..n_pieces{
        let mut s = range.start + t*segment_len;
        let mut e = range.start + min((t+1)*segment_len, range.len());
        if s >= range.end { // Happens e.g. if range.len() == 1 and n_pieces == 10
            s = range.end;
            e = range.end;
        }
        pieces.push(s..e); // Final segments may be empty. Is ok.
    }
    pieces
}

#[cfg(test)]
mod tests {
    use simple_sds_sbwt::ops::BitVec;

    use crate::util::bitvec_to_simple_sds_raw_bitvec;

    #[test]
    fn bitvec_to_simple_sds() {
        let mut bv = bitvec::bitvec![0; 567];
        for i in 0..bv.len() {
            if i % 3 == 0 {
                bv.set(i, true);
            }
        }

        let sds = bitvec_to_simple_sds_raw_bitvec(bv.clone());
        let sds = simple_sds_sbwt::bit_vector::BitVector::from(sds);
        for i in 0..bv.len() {
            assert_eq!(sds.get(i), bv[i]);
        }

        assert_eq!(bv.count_ones(), sds.count_ones());
    }

    #[test]
    fn all_ones() {
        // This exposes the bug in bitvec_to_simple_sds_raw_bitvec
        let bv = bitvec::bitvec![1; 6];
        let bv_clone = bv.clone();
        let sds = bitvec_to_simple_sds_raw_bitvec(bv);
        let sds = simple_sds_sbwt::bit_vector::BitVector::from(sds);

        assert_eq!(bv_clone.count_ones(), sds.count_ones());
    }
}