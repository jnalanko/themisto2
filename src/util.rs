use std::io::{Cursor, Read};

use bitvec::order::Lsb0;
use simple_sds_sbwt::serialize::Serialize;

pub(crate) fn bitvec_to_simple_sds_raw_bitvec(bv: bitvec::vec::BitVec<u64, Lsb0>) -> simple_sds_sbwt::raw_vector::RawVector {
    // Let's use the deserialization function in simple_sds_sbwt for a raw bitvector.
    // It requires the following header:
    let mut header = [0u64, 0u64]; // bits, words
    header[0] = bv.len() as u64; // Assumes little-endian byte order
    header[1] = bv.len().div_ceil(64) as u64;

    let header_bytes = bytemuck::cast_slice(&header);
    let raw_data = bytemuck::cast_slice(bv.as_raw_slice());
    let mut data_with_header = Cursor::new(header_bytes).chain(Cursor::new(raw_data));

    simple_sds_sbwt::raw_vector::RawVector::load(&mut data_with_header).unwrap()
}