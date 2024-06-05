#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufReader, BufWriter, Read}};

use sbwt::SubsetMatrix;

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
    );

    let matches = cli.get_matches();
    let matches = matches.subcommand_matches("build").unwrap();
    let sbwt_path = matches.get_one::<std::path::PathBuf>("sbwt").unwrap();
    let color_dump_path = matches.get_one::<std::path::PathBuf>("color-dump").unwrap();
    let out_path = matches.get_one::<std::path::PathBuf>("out").unwrap();
    let n_colors = *matches.get_one::<usize>("n-colors").unwrap();

    let mut sbwt_reader = BufReader::new(File::open(sbwt_path).unwrap());

    // Read header
    let header = SbwtFileHeader::read(&mut sbwt_reader).unwrap();
    assert!(header.has_lcs);

    // Read sbwt
    let sbwt = sbwt::SbwtIndex::<SubsetMatrix>::load(&mut sbwt_reader).unwrap();
    let lcs = sbwt::LcsArray::load(&mut sbwt_reader).unwrap();
    eprintln!("Loaded index with k = {}, precalc length = {}, # kmers = {}, # sbwt sets = {}", sbwt.k(), sbwt.get_lookup_table().prefix_length, sbwt.n_kmers(), sbwt.n_sets());

    let index = colored_kmers::ColoredKmers::new_from_themisto_color_dump(sbwt, lcs, &mut BufReader::new(File::open(color_dump_path).unwrap()), n_colors);
    index.serialize(&mut BufWriter::new(File::create(out_path).unwrap()));
}