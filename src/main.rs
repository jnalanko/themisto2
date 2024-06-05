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


fn main() {
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
    );

    let matches = cli.get_matches();
    if let Some(sub_matches) = matches.subcommand_matches("build"){
        let sbwt_path = sub_matches.get_one::<std::path::PathBuf>("sbwt").unwrap();
        let color_dump_path = sub_matches.get_one::<std::path::PathBuf>("color-dump").unwrap();
        let out_path = sub_matches.get_one::<std::path::PathBuf>("out").unwrap();
        let n_colors = *sub_matches.get_one::<usize>("n-colors").unwrap();

        let mut sbwt_reader = BufReader::new(File::open(sbwt_path).unwrap());

        // Read header
        let header = SbwtFileHeader::read(&mut sbwt_reader).unwrap();
        assert!(header.has_lcs);

        // Read sbwt
        let sbwt = sbwt::SbwtIndex::<SubsetMatrix>::load(&mut sbwt_reader).unwrap();
        let lcs = sbwt::LcsArray::load(&mut sbwt_reader).unwrap();
        log::info!("Loaded index with k = {}, precalc length = {}, # kmers = {}, # sbwt sets = {}", sbwt.k(), sbwt.get_lookup_table().prefix_length, sbwt.n_kmers(), sbwt.n_sets());

        let index = colored_kmers::ColoredKmers::new_from_themisto_color_dump(sbwt, lcs, &mut BufReader::new(File::open(color_dump_path).unwrap()), n_colors);
        index.serialize(&mut BufWriter::new(File::create(out_path).unwrap()));
    } else if let Some(sub_matches) = matches.subcommand_matches("pseudoalign"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();
        while let Some(rec) = reader.read_next().unwrap(){
            let intersection = index.intersection_pseudoalignment(rec.seq);
            println!("{}", intersection);
        }
    } else if let Some(sub_matches) = matches.subcommand_matches("pseudoalign-into-EM"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();
        log::info!("Computing compatibility vectors...");
        let mut compatibility_matrix = Vec::<BitVec>::new();
        while let Some(rec) = reader.read_next().unwrap(){
            compatibility_matrix.push(index.intersection_pseudoalignment(rec.seq));
        }
        log::info!("Running EM");
        let likelihood = SimpleLikelihood{};
        let n_colors = index.n_colors();
        let mixing_fractions = EM::fit_model(&likelihood, &compatibility_matrix, &vec![1.0 / n_colors as f64; n_colors]);
        println!("{:?}", &mixing_fractions);
    }
}