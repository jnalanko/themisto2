#![allow(clippy::manual_is_multiple_of)]

mod int_vec;

use std::{fs::File, io::{BufReader, Read}, path::Path};
use bitvec::prelude::*;

// A set of sets encoded as bitmaps.
#[derive(Clone, Debug)]
struct BitMaps {
    raw_data: Vec<u64>,
    n_colors: usize,
}

impl BitMaps {
    fn get_set(&self, set_id: usize) -> Vec<usize> {
        let start_bit = set_id * self.n_colors;
        let mut set = vec![];
        for color in 0..self.n_colors {
            let word = (start_bit+color) / 64;
            let off = (start_bit+color) % 64;
            if ((self.raw_data[word] >> off) & 1) == 1 {
                set.push(color);
            }
        }
        set
    }
}

struct SparseSets {
    concat: int_vec::CompactIntVec,
    starts: Vec<usize>,
}

impl SparseSets {
    fn get_set(&self, set_id: usize) -> Vec<usize> {
        let mut set = vec![];
        for i in self.starts[set_id]..self.starts[set_id+1] {
            set.push(self.concat.get(i));
        }
        set
    }
}

fn parse_sparse_file(filename: impl AsRef<Path>) -> SparseSets {
    let mut file = BufReader::new(File::open(filename).unwrap());

    let mut data_n_words = [0_u8; 8];
    let mut data_n_elements = [0_u8; 8];
    let mut bit_width = [0_u8; 8];
    let mut n_colors = [0_u8; 8];
    let mut n_sets = [0_u8; 8];

    file.read_exact(&mut data_n_words).unwrap();
    file.read_exact(&mut data_n_elements).unwrap();
    file.read_exact(&mut bit_width).unwrap();
    file.read_exact(&mut n_colors).unwrap();
    file.read_exact(&mut n_sets).unwrap();

    let data_n_words = usize::from_le_bytes(data_n_words);
    let data_n_elements = usize::from_le_bytes(data_n_elements);
    let bit_width = usize::from_le_bytes(bit_width);
    let n_colors = usize::from_le_bytes(n_colors);
    let n_sets = usize::from_le_bytes(n_sets);

    eprintln!("Sparse metadata: data_n_words={}, data_n_elements={}, bit_width={}, n_colors={}, n_sets={}",
        data_n_words, data_n_elements, bit_width, n_colors, n_sets);

    let mut raw_data: Vec<u64> = vec![0; data_n_words];
    file.read_exact(bytemuck::cast_slice_mut(raw_data.as_mut_slice())).unwrap();

    let mut starts: Vec<usize> = vec![0; n_sets+1];
    file.read_exact(bytemuck::cast_slice_mut(starts.as_mut_slice())).unwrap();

    let concat = int_vec::CompactIntVec::from_parts(raw_data, data_n_elements, bit_width);

    SparseSets {
        concat,
        starts,
    }
}

fn parse_dense_file(filename: impl AsRef<Path>) -> BitMaps {
    let mut file = BufReader::new(File::open(filename).unwrap());

    let mut data_n_words = [0_u8; 8];
    let mut n_bits = [0_u8; 8];
    let mut n_colors = [0_u8; 8];

    file.read_exact(&mut data_n_words).unwrap();
    file.read_exact(&mut n_bits).unwrap();
    file.read_exact(&mut n_colors).unwrap();

    let data_n_words = usize::from_le_bytes(data_n_words);
    let n_bits = usize::from_le_bytes(n_bits);
    let n_colors = usize::from_le_bytes(n_colors);

    eprintln!("Dense metadata: data_n_words={}, n_bits={}, n_colors={}",
        data_n_words, n_bits, n_colors);

    let mut raw_data: Vec<u64> = vec![0; data_n_words];
    file.read_exact(bytemuck::cast_slice_mut(raw_data.as_mut_slice())).unwrap();

    assert!(n_bits % n_colors == 0);

    BitMaps {
        raw_data,
        n_colors,
    }

}

fn parse_marks_file(filename: impl AsRef<Path>) -> BitVec<u64, Lsb0> {
    let mut file = BufReader::new(File::open(filename).unwrap());

    let mut n_words = [0_u8; 8];
    let mut n_bits = [0_u8; 8];

    file.read_exact(&mut n_words).unwrap();
    file.read_exact(&mut n_bits).unwrap();

    let n_words = usize::from_le_bytes(n_words);
    let n_bits = usize::from_le_bytes(n_bits);

    eprintln!("Marks metadata: n_words={}, n_bits={}",
        n_words, n_bits);

    let mut raw_data: Vec<u64> = vec![0; n_words];
    file.read_exact(bytemuck::cast_slice_mut(raw_data.as_mut_slice())).unwrap();

    // Load marks into bitvec
    let mut bits: BitVec<u64, Lsb0> = BitVec::from_vec(raw_data);
    bits.truncate(n_bits);
    bits
}

// Metadata: data_n_words=4790, data_n_elements=51084, bit_width=6, n_colors=43, n_sets=10260


fn main() {
    // CLI that takes an input file prefix
    let cli = clap::Command::new("Themisto debug extractor")
        .arg(clap::Arg::new("input_prefix")
            .help("Input file prefix")
            .required(true)
            .short('i')
            .value_parser(clap::value_parser!(String)));
    
    let matches = cli.get_matches();
    let input_prefix = matches.get_one::<String>("input_prefix").unwrap();

    let mut sparse_filename = input_prefix.clone();
    sparse_filename.push_str(".sparse");

    let mut dense_filename = input_prefix.clone();
    dense_filename.push_str(".dense");

    let mut marks_filename = input_prefix.clone();
    marks_filename.push_str(".marks");

    let sparse_sets = parse_sparse_file(sparse_filename);
    let bitmaps = parse_dense_file(dense_filename);
    let is_dense_marks = parse_marks_file(marks_filename);

    let print = |set_id: usize, set: &[usize]| {
        print!("color_set_id={} size={}", set_id, set.len());
        for c in set {
            print!(" {}", c);
        }
        println!();
    };
    let mut dense_id = 0_usize;
    let mut sparse_id = 0_usize;
    for set_id in 0..is_dense_marks.len() {
        if is_dense_marks[set_id] {
            print(set_id, &bitmaps.get_set(dense_id));
            dense_id += 1;
        } else {
            print(set_id, &sparse_sets.get_set(sparse_id));
            sparse_id += 1;
        }
    }

}