#![allow(clippy::manual_is_multiple_of)]

mod int_vec;

use std::{fs::File, io::{BufReader, Read}, path::Path};

fn parse_sparse_file(filename: impl AsRef<Path>) {
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
    for set_id in 0..n_sets {
        print!("{}:", set_id);
        for i in starts[set_id]..starts[set_id+1] {
            print!(" {}", concat.get(i));
        }
        println!();
    }
}

fn parse_dense_file(filename: impl AsRef<Path>) {
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
    let n_sets = n_bits / n_colors;
    for set_id in 0..n_sets {
        print!("{}:", set_id);
        let start_bit = set_id * n_colors;
        for color in 0..n_colors {
            let word = (start_bit+color) / 64;
            let off = (start_bit+color) % 64;
            if ((raw_data[word] >> off) & 1) == 1 {
                print!(" {}", color); 
            }
        }
        println!();
    }


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

    parse_sparse_file(sparse_filename);
}