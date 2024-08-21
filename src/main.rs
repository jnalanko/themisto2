#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufReader, BufWriter, Read}};
use bitvec::prelude::*;

use sbwt::SubsetMatrix;

use crate::EM::fit_model;

mod EM;
mod colored_kmers;

const FILE_FORMAT_STRING: &[u8] = b"sbwtfile-v1";

struct SbwtFileHeader {
    has_lcs: bool,
}

impl SbwtFileHeader {
    fn read<R: Read>(input: &mut R) -> std::io::Result<SbwtFileHeader> {
        read_and_check_string(input, FILE_FORMAT_STRING, "Invalid or incompatible file format").unwrap();
        let has_lcs: bool = byteorder::ReadBytesExt::read_u8(input).unwrap() != 0;
        Ok(Self{has_lcs})
    }

}

struct SimpleLikelihood {} // Based on compatibility vectors

impl EM::Likelihood for SimpleLikelihood {
    type Observation = BitVec; // Compatibility vector
    fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64 {
        if *x_i.get(k).unwrap() {
            0.99
        }  else {
            0.01
        }
    }
}

// Read a byte string in this format: first a little-endian usize giving the length,
// then the bytes. Check that the bytes match the given slice. Returns an IO error with
// the given error message if the strings do not match.
fn read_and_check_string<R: std::io::Read>(input: &mut R, should_be_this: &[u8], error_message: &str) -> std::io::Result<()> {
    let mut len_buf = [0_u8; 8]; 
    input.read_exact(&mut len_buf)?;

    let len = usize::from_le_bytes(len_buf);
    if len != should_be_this.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            error_message
        ));
    }

    let mut string_buf = vec![0u8; len];
    input.read_exact(&mut string_buf)?;

    if string_buf != should_be_this {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            error_message
        ));
    }

    Ok(())

}

// Removes duplicates and returns the count of each distinct element remaining
// in the vector.
fn reduce_to_classes(compatibility_vectors: &mut Vec<BitVec>) -> Vec<usize> {
    compatibility_vectors.sort_unstable();
    
    let n = compatibility_vectors.len();
    let mut counts = Vec::<usize>::new();
    let mut i = 0_usize;
    let mut j = 0_usize;
    let mut insert_pos = 0_usize;

    let mut borrow_checker_workaround = bitvec![0; compatibility_vectors[0].len()];
    while i < n {
        while j + 1 < n && compatibility_vectors[j+1] == compatibility_vectors[i] {
            j += 1;
        }
        j += 1;
        // vectors [i, j) are the same
        counts.push(j-i);

        borrow_checker_workaround.copy_from_bitslice(compatibility_vectors[i].as_bitslice());
        compatibility_vectors[insert_pos].copy_from_bitslice(&borrow_checker_workaround);
        i = j;
        insert_pos += 1;
    }

    compatibility_vectors.resize(insert_pos, bitvec![]);
    compatibility_vectors.shrink_to_fit();

    counts

}



fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let cli = clap::Command::new("themisto2")
    .arg_required_else_help(true) 
    .subcommand(clap::Command::new("build")
        .arg_required_else_help(true)
        .arg(clap::Arg::new("sbwt")
            .long("sbwt")
            .short('s')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("color-dump")
            .long("color-dump")
            .short('c')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("n-colors")
            .long("n-colors")
            .short('n')
            .value_parser(clap::value_parser!(usize))
            .required(true)
        )
        .arg(clap::Arg::new("out")
            .long("out")
            .short('o')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
    )
    .subcommand(clap::Command::new("pseudoalign")
        .arg_required_else_help(true)
        .arg(clap::Arg::new("index")
            .long("index")
            .short('i')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("query")
            .long("query")
            .short('q')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
    )
    .subcommand(clap::Command::new("pseudoalign-into-EM")
        .arg_required_else_help(true)
        .arg(clap::Arg::new("index")
            .long("index")
            .short('i')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("query")
            .long("query")
            .short('q')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("min-hits")
            .long("min-hits")
            .short('m')
            .value_parser(clap::value_parser!(usize))
            .required(true)
        )
        .arg(clap::Arg::new("n-threads")
            .long("n-threads")
            .short('t')
            .value_parser(clap::value_parser!(usize))
            .default_value("1")
        )
    );

    let matches = cli.get_matches();
    if let Some(sub_matches) = matches.subcommand_matches("import"){
        let sbwt_ascii_dump = sub_matches.get_one::<std::path::PathBuf>("sbwt").unwrap();
        let themisto_dump_prefix = sub_matches.get_one::<std::path::PathBuf>("color-dump-prefix").unwrap();
        let out_path = sub_matches.get_one::<std::path::PathBuf>("out").unwrap();

        let unitig_filename = format!("{}.unitigs.fa", themisto_dump_prefix.to_str().unwrap());
        let color_sets_filename = format!("{}.color_sets.txt", themisto_dump_prefix.to_str().unwrap());
        let metadata_filename = format!("{}.metadata.txt", themisto_dump_prefix.to_str().unwrap());

        let sbwt_in = BufReader::new(File::open(sbwt_ascii_dump).unwrap());
        let unitigs_in = BufReader::new(File::open(unitig_filename).unwrap());
        let color_sets_in = BufReader::new(File::open(color_sets_filename).unwrap());
        let metadata_in = BufReader::new(File::open(metadata_filename).unwrap());

        let mut out = BufWriter::new(File::create(out_path).unwrap());

        let index = colored_kmers::ColoredKmers::new_from_new_themisto_index_dump(sbwt_in, metadata_in, unitigs_in, color_sets_in, 0);
        index.serialize(&mut out);
    } else if let Some(sub_matches) = matches.subcommand_matches("pseudoalign"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        
        log::info!("Loading index");
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();

        log::info!("Pseudoaligning");
        while let Some(rec) = reader.read_next().unwrap(){
            let intersection = index.intersection_pseudoalignment(rec.seq, 1);
            println!("{}", intersection);
        }
    } else if let Some(sub_matches) = matches.subcommand_matches("pseudoalign-into-EM"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        let min_hits = *sub_matches.get_one::<usize>("min-hits").unwrap();
        let n_threads = *sub_matches.get_one::<usize>("n-threads").unwrap();

        log::info!("Loading index");
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();
        log::info!("Computing compatibility vectors");
        let mut compatibility_matrix = Vec::<BitVec>::new();
        while let Some(rec) = reader.read_next().unwrap(){
            compatibility_matrix.push(index.intersection_pseudoalignment(rec.seq, min_hits));
        }
        log::info!("Constructed {} compatibility vectors", compatibility_matrix.len());

        log::info!("Reducing to equivalence classes");
        let counts = reduce_to_classes(&mut compatibility_matrix);
        log::info!("Reduced to {} equivalence classes", counts.len());

        log::info!("Running EM");
        let likelihood = SimpleLikelihood{};
        let n_colors = index.n_colors();
        let mixing_fractions = EM::fit_model(&likelihood, &compatibility_matrix, &counts, &vec![1.0 / n_colors as f64; n_colors], n_threads);
        println!("{:?}", &mixing_fractions);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_to_classes(){
        let mut vec = vec![
            bitvec![1, 1, 1],
            bitvec![1, 0, 1],
            bitvec![0, 1, 0],
            bitvec![1, 0, 1],
            bitvec![1, 0, 0],
            bitvec![0, 1, 0],
            bitvec![0, 1, 0],
            bitvec![1, 0, 1]];
        let counts = reduce_to_classes(&mut vec);
        assert_eq!(vec, vec![
            bitvec![0, 1, 0],
            bitvec![1, 0, 0],
            bitvec![1, 0, 1],
            bitvec![1, 1, 1]]);
        assert_eq!(counts, vec![3, 1, 3, 1]);
    }
}