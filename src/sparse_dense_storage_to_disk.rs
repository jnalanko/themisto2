use std::{cmp::min, fs::File, io::{Read, Seek, Write}, ops::Range};

use crate::{coloring_interface::{ColorSetStorage, ColorSetView}, int_vec::CompactIntVec, sparse_dense_storage::{SparseDenseColorSetView, SparseDenseStorage}};
use bitvec::prelude::*;
use simple_sds_sbwt::ops::{Select, SelectZero};

#[derive(Clone, Debug)]
struct BitMaps {
    bitmap_data: bitvec::vec::BitVec, // Concatenation of bit vectors
    individual_length: usize, // Length of each bitmap in bitmap_data
}

impl BitMaps {
    fn new_with_zero_init(individual_length: usize, n_sets: usize) -> Self {
        BitMaps{bitmap_data: bitvec![0; n_sets*individual_length], individual_length}
    }

    fn get_mut(&mut self, bitmap_idx: usize) -> &mut BitSlice {
        &mut self.bitmap_data[bitmap_idx*self.individual_length .. (bitmap_idx + 1) * self.individual_length]
    }
}

pub fn write_piece(piece: SparseDenseStorage, color_id_range: Range<usize>, sparse_file: &mut File, dense_file: &mut File, sparse_set_insertion_points: &mut [usize]) {
    write_sparse_sets_piece(&piece, color_id_range.clone(), sparse_file, sparse_set_insertion_points);
    write_dense_sets_piece(&piece, color_id_range, dense_file);
}


fn write_dense_sets_piece(piece: &SparseDenseStorage, color_id_range: Range<usize>, dense_file: &mut File) {
    dense_file.rewind().unwrap();

    // Copy-paste from above:
    // * A file for dense sets: `[data_n_words: usize][total_bitvec_len: usize][n_colors: usize][data]`

    // Read metadata
    let mut metadata = [0u64; 3];
    dense_file.read_exact(bytemuck::cast_slice_mut(&mut metadata)).unwrap();
    let data_n_words = metadata[0] as usize;
    let total_bitvec_len = metadata[1] as usize;
    let total_n_colors = metadata[2] as usize;

    assert!(total_bitvec_len % total_n_colors == 0);
    let total_n_sets = total_bitvec_len / total_n_colors;

    let file_raw_data_start_offset: usize = 8*3;

    let max_n_sets_in_buf = 1000_usize.next_multiple_of(64);
    let max_n_bits_in_buf = max_n_sets_in_buf * total_n_colors;
    let mut buf_bitmap = BitMaps::new_with_zero_init(total_n_colors, min(max_n_sets_in_buf, total_n_sets));

    let mut file_offset = file_raw_data_start_offset;

    // Read raw data from disk and put to a bitmap
    let buf_bytes: &mut [u8] = bytemuck::cast_slice_mut(buf_bitmap.bitmap_data.as_raw_mut_slice());
    dense_file.read_exact(buf_bytes).unwrap();
    let mut real_bytes_in_buf = buf_bytes.len();

    let mut n_bits_in_past_buffers = 0_usize;
    for (dense_id, color_set_id) in piece.get_dense_marks().one_iter() {
        let set_view = piece.get_set_view(color_set_id); // This implies a rank query, which is unnecessary since we're scanning the bit vector
        let piece_bit_slice = match set_view {
            SparseDenseColorSetView::Dense(bit_slice) => bit_slice,
            SparseDenseColorSetView::Sparse(_) => panic!("Expected dense set, got sparse"),
        };
        assert_eq!(piece_bit_slice.len(), color_id_range.len());
        let mut start_bit = dense_id*total_n_colors - n_bits_in_past_buffers;
        while start_bit >= max_n_bits_in_buf {
            // Write the current raw data back to disk, read the next chunk of raw data
            let raw_data = buf_bitmap.bitmap_data.as_raw_mut_slice();
            let raw_bytes: &mut [u8] = bytemuck::cast_slice_mut(raw_data);
            dense_file.seek(std::io::SeekFrom::Start(file_offset as u64)).unwrap();
            dense_file.write_all(&raw_bytes[0..real_bytes_in_buf]).unwrap();
            file_offset += real_bytes_in_buf;
            n_bits_in_past_buffers += raw_bytes.len()*8;

            dense_file.seek(std::io::SeekFrom::Start(file_offset as u64)).unwrap();
            assert!(n_bits_in_past_buffers % 64 == 0);
            let words_remaining = data_n_words - n_bits_in_past_buffers / 64; 
            let bytes_remaining = words_remaining * 8;
            let bytes_to_read = min(bytes_remaining, raw_bytes.len());
            dense_file.read_exact(&mut raw_bytes[0..bytes_to_read]).unwrap();
            real_bytes_in_buf = bytes_to_read;
            start_bit -= max_n_bits_in_buf;
        }
        let target_set = buf_bitmap.get_mut(start_bit / total_n_colors);
        let target_range = &mut target_set[color_id_range.clone()];
        assert!(target_range.count_ones() == 0); // Must not have been set before
        target_range.copy_from_bitslice(piece_bit_slice);
    }

    // Remaining buffer
    let raw_data = buf_bitmap.bitmap_data.as_raw_mut_slice();
    let raw_bytes: &mut [u8] = bytemuck::cast_slice_mut(raw_data);
    dense_file.seek(std::io::SeekFrom::Start(file_offset as u64)).unwrap();
    dense_file.write_all(&raw_bytes[0..real_bytes_in_buf]).unwrap();
}

fn write_sparse_sets_piece(piece: &SparseDenseStorage, color_id_range: Range<usize>, sparse_file: &mut File, sparse_set_insertion_points: &mut [usize]) {

    // Write sparse sets.
    // File format copied from above:
    // * A file for sparse sets: `[data_n_words: usize][data_n_elements: usize][bit_width: usize][n_colors: usize][n_sets: usize][data][starts]`

    sparse_file.rewind().unwrap();

    // Read metadata
    let mut metadata = [0u64; 5];
    sparse_file.read_exact(bytemuck::cast_slice_mut(&mut metadata)).unwrap();
    let data_n_words = metadata[0] as usize;
    let data_n_elements = metadata[1] as usize;
    let bit_width = metadata[2] as usize;
    let _n_colors = metadata[3] as usize;
    let _n_sets = metadata[4] as usize;

    let file_raw_data_start_offset: usize = 8*5;

    // We want a buffer that can hold a number of elements m such that m*bit_width is a multiple of
    // 64. This guarantees that the last word of the buffer is aligned with the
    // end of the bits of the last element. Unless it's the last buffer, but that's ok.
    //log::warn!("USING SMALL BUFFER FOR DEBUG PURPOSES");
    let max_buf_cap_elements = 1_000_000_usize.next_multiple_of(64);
    let buf_cap_elements = min(max_buf_cap_elements, data_n_elements);
    let mut buf_compact_int_vec = CompactIntVec::new(buf_cap_elements, bit_width);
    let mut file_offset = file_raw_data_start_offset;

    // Read raw data from disk and put to buf_compact_int_vec
    let buf_words: &mut [u64] = buf_compact_int_vec.get_mut_raw_data();
    let buf_bytes: &mut [u8] = bytemuck::cast_slice_mut(buf_words);
    sparse_file.read_exact(buf_bytes).unwrap();
    let mut real_bytes_in_buf = buf_bytes.len();
    
    let mut n_elements_in_past_buffers = 0_usize;
    for (sparse_id, color_set_id) in piece.get_dense_marks().zero_iter() {
        for color in piece.get_set_view(color_set_id).iter() { // TODO: this does an unnecessary rank.
            let mut buf_insertion_point = sparse_set_insertion_points[sparse_id] - n_elements_in_past_buffers;
            while buf_insertion_point >= buf_cap_elements {
                let buf_words: &mut [u64] = buf_compact_int_vec.get_mut_raw_data();
                let buf_bytes: &mut [u8] = bytemuck::cast_slice_mut(buf_words);

                sparse_file.seek(std::io::SeekFrom::Start(file_offset as u64)).unwrap();
                sparse_file.write_all(&buf_bytes[0..real_bytes_in_buf]).unwrap();

                file_offset += buf_bytes.len();
                n_elements_in_past_buffers += buf_cap_elements;

                sparse_file.seek(std::io::SeekFrom::Start(file_offset as u64)).unwrap();
                assert!(n_elements_in_past_buffers*bit_width % 64 == 0);
                let words_remaining = data_n_words - n_elements_in_past_buffers*bit_width / 64; 
                let bytes_remaining = words_remaining * 8;
                let bytes_to_read = min(bytes_remaining, buf_bytes.len());
                sparse_file.read_exact(&mut buf_bytes[0..bytes_to_read]).unwrap();
                real_bytes_in_buf = bytes_to_read;

                buf_insertion_point = sparse_set_insertion_points[sparse_id] - n_elements_in_past_buffers;
            }
            assert!(buf_compact_int_vec.get(buf_insertion_point) == 0); // This must not have been written yet
            buf_compact_int_vec.set(buf_insertion_point, color + color_id_range.start);
            sparse_set_insertion_points[sparse_id] += 1;
        }
    }

    // Remaining buffer
    let buf_words: &mut [u64] = buf_compact_int_vec.get_mut_raw_data();
    let buf_bytes: &mut [u8] = bytemuck::cast_slice_mut(buf_words);

    sparse_file.seek(std::io::SeekFrom::Start(file_offset as u64)).unwrap();
    sparse_file.write_all(&buf_bytes[0..real_bytes_in_buf]).unwrap();

}