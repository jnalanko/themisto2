#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufReader, BufWriter}};
use bitvec::prelude::*;

mod EM;
mod colored_kmers;

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

struct LikelihoodMatrix {
    matrix: Vec<Vec<f64>>,
}

impl EM::Likelihood for LikelihoodMatrix {
    type Observation = usize; // Index of the read
    fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64 {
        self.matrix[*x_i][k]
    }
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
    .subcommand(clap::Command::new("import")
        .arg_required_else_help(true)
        .arg(clap::Arg::new("sbwt-ascii-dump")
            .long("sbwt-ascii-dump")
            .short('s')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("color-dump-prefix")
            .long("color-dump-prefix")
            .short('c')
            .value_parser(clap::value_parser!(std::path::PathBuf))
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
        .arg(clap::Arg::new("numerator")
            .long("numerator")
            .value_parser(clap::builder::PossibleValuesParser::new(vec!["hits", "distinguishing"]))
            .default_value("hits")
        )
        .arg(clap::Arg::new("denominator")
            .long("denominator")
            .value_parser(clap::builder::PossibleValuesParser::new(vec!["all", "relevant", "max-distinguishing"]))
            .default_value("all")
        )
        .arg(clap::Arg::new("likelihood")
            .long("likelihood")
            .value_parser(clap::builder::PossibleValuesParser::new(vec!["linear", "softmax", "99p", "beta-binomial"]))
            .default_value("linear")
        )
    );

    let matches = cli.get_matches();
    if let Some(sub_matches) = matches.subcommand_matches("import"){
        let sbwt_ascii_dump = sub_matches.get_one::<std::path::PathBuf>("sbwt-ascii-dump").unwrap();
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
        let numerator = sub_matches.get_one::<String>("numerator").unwrap();
        let denominator = sub_matches.get_one::<String>("denominator").unwrap();

        log::info!("Loading index");
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();

        let mut likelihood_matrix = Vec::<Vec<f64>>::new();
        while let Some(rec) = reader.read_next().unwrap(){
            let data = index.compute_pseudoalignment_data(rec.seq, min_hits);
            
            let mut row: Vec<f64> = vec![0.0; index.n_colors()];
            for color in 0..index.n_colors() {
                let numerator_value = match numerator.as_str() {
                    "hits" => data.hit_counts[color],
                    "distinguishing" => data.distinguishing_hit_counts[color],
                    _ => panic!("Invalid numerator {}", numerator)
                };
                let n_kmers = std::cmp::max(0, rec.seq.len() as isize - index.get_k() as isize + 1) as usize;
                let mut denominator_value = match denominator.as_str() {
                    "all" => n_kmers,
                    "relevant" => data.n_relevant_kmers,
                    "max-distinguishing" => *data.distinguishing_hit_counts.iter().max().unwrap_or(&0),
                    _ => panic!("Invalid denominator {}", denominator)
                };
                denominator_value = std::cmp::max(1, denominator_value); // Avoid division by zero

                row[color] = numerator_value as f64 / denominator_value as f64;
            }

            // Add zero inflation
            row.iter_mut().for_each(|x| *x = x.max(0.01));

            // Normalize
            let rowsum = row.iter().sum::<f64>();
            row.iter_mut().for_each(|x| *x /= rowsum);

            likelihood_matrix.push(row);
        }

        // Observation is now a likelihood matrix row
        let n_rows = likelihood_matrix.len();
        let likelihood = LikelihoodMatrix{matrix: likelihood_matrix};
        let observations: Vec<usize> = (0..n_rows).collect();
        let observation_counts: Vec<usize> = vec![1; n_rows];

        let mixing_fractions = EM::fit_model(&likelihood, &observations, &observation_counts, &vec![1.0 / n_rows as f64; index.n_colors()], n_threads);
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