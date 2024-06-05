#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufReader, BufWriter}};

use sbwt::SubsetMatrix;

mod EM;
mod colored_kmers;

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

    let sbwt = sbwt::SbwtIndex::<SubsetMatrix>::load(&mut BufReader::new(File::open(sbwt_path).unwrap())).unwrap();
    let index = colored_kmers::ColoredKmers::new_from_themisto_color_dump(sbwt, &mut BufReader::new(File::open(color_dump_path).unwrap()), n_colors);
    index.serialize(&mut BufWriter::new(File::create(out_path).unwrap()));
}