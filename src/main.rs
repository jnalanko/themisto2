#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufRead, BufReader, BufWriter}, path::PathBuf};
use bitvec::prelude::*;
use colored_kmers::ColoredKmers;

mod EM;
mod colored_kmers;
mod themisto1_compatibility;

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

fn softmax(values: &[f64]) -> Vec<f64> {
    let exp_values: Vec<f64> = values.iter().map(|&x| x.exp()).collect();
    let sum_exp_values: f64 = exp_values.iter().sum();
    exp_values.iter().map(|&x| x / sum_exp_values).collect()
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
        .arg(clap::Arg::new("input")
            .help("A file with one fasta/fastq filename per line")
            .short('i')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("output")
            .help("Outfile filename")
            .short('o')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("temp-dir")
            .help("Directory for temporary files")
            .long("temp-dir")
            .short('d')
            .value_parser(clap::value_parser!(std::path::PathBuf))
            .required(true)
        )
        .arg(clap::Arg::new("k")
            .short('k')
            .value_parser(clap::value_parser!(usize))
            .required(true)
        )
        .arg(clap::Arg::new("n-threads")
            .short('t')
            .value_parser(clap::value_parser!(usize))
            .default_value("4")
        )
    )
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
    .subcommand(clap::Command::new("intersection-pseudoalign")
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
    .subcommand(clap::Command::new("dump-pseudoalignment-data")
        .about("Dumps pseudoalignment data for each read in JSON format")
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
        .arg(clap::Arg::new("max-iterations")
            .long("max-iterations")
            .value_parser(clap::value_parser!(usize))
            .default_value("2000")
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
            .value_parser(clap::builder::PossibleValuesParser::new(vec!["linear", "softmax", "99p"]))
            .default_value("linear")
        )
    )
    .subcommand(clap::Command::new("intersection-pseudoalign-into-EM")
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
        .arg(clap::Arg::new("max-iterations")
            .long("max-iterations")
            .value_parser(clap::value_parser!(usize))
            .default_value("2000")
        )
        .arg(clap::Arg::new("initial-likelihood-ratio")
            .long("initial-likelihood-ratio")
            .short('w')
            .value_parser(clap::value_parser!(f64))
            .default_value("99")
        )
        .arg(clap::Arg::new("static-likelihood")
            .long("static-likelihood")
            .action(clap::ArgAction::SetTrue)
        )
    );

    let matches = cli.get_matches();
    if let Some(sub_matches) = matches.subcommand_matches("build"){
        let input_fof = sub_matches.get_one::<std::path::PathBuf>("input").unwrap();
        let out_path = sub_matches.get_one::<std::path::PathBuf>("output").unwrap();
        let k = *sub_matches.get_one::<usize>("k").unwrap();
        let n_threads = *sub_matches.get_one::<usize>("n-threads").unwrap();
        let temp_dir = sub_matches.get_one::<PathBuf>("temp-dir").unwrap();
        let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
        let index = ColoredKmers::new(input_paths.as_slice(), k, n_threads, temp_dir);
        index.serialize(&mut BufWriter::new(File::create(out_path).unwrap()));

    } else if let Some(sub_matches) = matches.subcommand_matches("import"){
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
    } else if let Some(sub_matches) = matches.subcommand_matches("intersection-pseudoalign"){
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
    } else if let Some(sub_matches) = matches.subcommand_matches("intersection-pseudoalign-into-EM"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        let min_hits = *sub_matches.get_one::<usize>("min-hits").unwrap();
        let n_threads = *sub_matches.get_one::<usize>("n-threads").unwrap();
        let max_iterations = *sub_matches.get_one::<usize>("max-iterations").unwrap();
        let initial_w = *sub_matches.get_one::<f64>("initial-likelihood-ratio").unwrap();
        let optimize_w = !(sub_matches.get_flag("static-likelihood"));
        
        log::info!("Loading index");
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();

        log::info!("Pseudoaligning");
        let mut intersections = Vec::<BitVec>::new(); // Todo: store more compactly with 8 bits per byte?
        while let Some(rec) = reader.read_next().unwrap(){
            let intersection = index.intersection_pseudoalignment(rec.seq, min_hits);
            intersections.push(intersection);
        }

        let class_counts = reduce_to_classes(&mut intersections);

        // Represent with one byte per bit for compatibility with the EM algorithm
        let intersections: Vec<Vec<u8>> = intersections.iter().map(|v| v.iter().map(|b| *b as u8).collect()).collect();

        let (thetas, w) = EM::fit_model_with_intersection_inputs(&intersections, &class_counts, &vec![1.0 / index.n_colors() as f64; index.n_colors()], initial_w, optimize_w, n_threads, max_iterations);
        println!("Final likelihood ratio w: {}", w);
        println!("Mixing fractions: {:?}", thetas);

    } else if let Some(sub_matches) = matches.subcommand_matches("dump-pseudoalignment-data"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        let min_hits = *sub_matches.get_one::<usize>("min-hits").unwrap();
        
        log::info!("Loading index");
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();

        log::info!("Pseudoaligning");
        while let Some(rec) = reader.read_next().unwrap(){
            let data = index.compute_pseudoalignment_data(rec.seq, min_hits);
            let json = serde_json::to_string(&data).unwrap();
            println!("{}", json);
        }
    } else if let Some(sub_matches) = matches.subcommand_matches("pseudoalign-into-EM"){
        let index_path = sub_matches.get_one::<std::path::PathBuf>("index").unwrap();
        let query_path = sub_matches.get_one::<std::path::PathBuf>("query").unwrap();
        let min_hits = *sub_matches.get_one::<usize>("min-hits").unwrap();
        let n_threads = *sub_matches.get_one::<usize>("n-threads").unwrap();
        let numerator = sub_matches.get_one::<String>("numerator").unwrap();
        let denominator = sub_matches.get_one::<String>("denominator").unwrap();
        let likelihood_type = sub_matches.get_one::<String>("likelihood").unwrap();
        let max_iterations = *sub_matches.get_one::<usize>("max-iterations").unwrap();

        log::info!("Loading index");
        let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
        log::info!("Loaded index with {} distinct k-mers and {} colors", index.n_kmers(), index.n_colors());

        let mut reader = jseqio::reader::DynamicFastXReader::from_file(query_path).unwrap();

        let likelihood_function: Box<dyn Fn(&[f64]) -> Vec<f64>> = match likelihood_type.as_str() { // Takes a row of scores, returns a row of likelihoods. That is, f: R^n -> R^n, where n is the number of colors
            "linear" => Box::new(|v: &[f64]| v.to_vec()), // Identity function
            "99p" => Box::new(|v: &[f64]| {
                let (argmax, _max) = v.iter().enumerate().max_by(|(_, a),(_, b)| a.total_cmp(b)).unwrap();
                let mut answer: Vec<f64> = vec![0.01; index.n_colors()];
                answer[argmax] = 0.99;
                answer
            }),
            //"betabinomial" => Box::new(|_: &[f64]| {
            //    todo!(); // Issue: Beta binomial takes in an integer, not a float. But it's almost linear with our hyperparameters, so linear works as a good substitute for this.
            //}),
            "softmax" => Box::new(|v: &[f64]| softmax(v)),
            _ => panic!("Invalid likelihood type: {}", likelihood_type)
        };
        let mut likelihood_matrix = Vec::<Vec<f64>>::new();
        while let Some(rec) = reader.read_next().unwrap(){
            let data = index.compute_pseudoalignment_data(rec.seq, min_hits);
            
            let mut row: Vec<f64> = vec![0.0; index.n_colors()];
            for color in 0..index.n_colors() {
                let numerator_value = match numerator.as_str() {
                    "hits" => data.hit_counts[color],
                    "distinguishing" => data.distinguishing_hit_counts[color],
                    _ => panic!("Invalid numerator type: {}", numerator)
                };
                let n_kmers = std::cmp::max(0, rec.seq.len() as isize - index.get_k() as isize + 1) as usize;
                let mut denominator_value = match denominator.as_str() {
                    "all" => n_kmers,
                    "relevant" => data.n_relevant_kmers,
                    "max-distinguishing" => *data.distinguishing_hit_counts.iter().max().unwrap_or(&0),
                    _ => panic!("Invalid denominator type: {}", denominator)
                };
                denominator_value = std::cmp::max(1, denominator_value); // Avoid division by zero

                row[color] = numerator_value as f64 / denominator_value as f64;
            }

            // Add zero inflation
            row.iter_mut().for_each(|x| *x = x.max(0.01));

            // Apply the likelihood function
            row = likelihood_function(&row);

            // Normalize
            let rowsum = row.iter().sum::<f64>();
            row.iter_mut().for_each(|x| *x /= rowsum);

            likelihood_matrix.push(row);
        }

        likelihood_matrix.shrink_to_fit(); // Saving some memory

        // Observation is now a likelihood matrix row
        let n_rows = likelihood_matrix.len();
        let likelihood = LikelihoodMatrix{matrix: likelihood_matrix};
        let observations: Vec<usize> = (0..n_rows).collect();
        let observation_counts: Vec<usize> = vec![1; n_rows];

        let mixing_fractions = EM::fit_model(&likelihood, &observations, &observation_counts, &vec![1.0 / index.n_colors() as f64; index.n_colors()], n_threads, max_iterations);
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