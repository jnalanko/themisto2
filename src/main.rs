#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufRead, BufReader, BufWriter}, path::PathBuf, sync::Arc};
use bitvec::prelude::*;
use clap::{builder::PossibleValuesParser, Parser, Subcommand};
use colored_kmers::ColoredKmers;
use compatibility_criteria::unique_support_combination_method;
use sbwt::{LcsArray, SbwtIndexVariant};

mod EM;
mod colored_kmers;
mod themisto1_compatibility;
mod compatibility_criteria;
mod compact_colored_kmers;

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

#[derive(Parser)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommands,
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Denominator {
    All,
    Relevant,
    MaxHits,
}

/*
impl FromStr for Denominator {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "all" => Ok(Denominator::All),
            "relevant" => Ok(Denominator::Relevant),
            "maxhits" => Ok(Denominator::MaxHits),
            _ => Err(format!("Invalid denominator value: {}", s)),
        }
    }
}
*/

#[derive(Subcommand)]
pub enum Subcommands {
    #[command(arg_required_else_help = true)]
    Build {
        #[arg(help = "A file with one fasta/fastq filename per line", short, long, required = true)]
        input: PathBuf,

        #[arg(help = "Output filename", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Directory for temporary files", short = 'd', long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(short, required = true)]
        k: usize,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize
    },

    #[command(arg_required_else_help = true)]
    Import {
        #[arg(long = "sbwt-ascii-dump", short = 's', required = true)]
        sbwt_ascii_dump: PathBuf,

        #[arg(long = "color-dump-prefix", short = 'c', required = true)]
        color_dump_prefix: PathBuf,

        #[arg(long = "out", short = 'o', required = true)]
        out: PathBuf,
    },

    #[command(arg_required_else_help = true, name = "intersection-pseudoalign")]
    IntersectionPseudoalign {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,

        #[arg(long = "min-hits", short = 'm', default_value = "1")]
        min_hits: usize,
    },

    #[command(arg_required_else_help = true, name = "threshold-pseudoalign")]
    ThresholdPseudoalign {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,

        #[arg(long = "min-hits", short = 'm', required = true, default_value = "1")]
        min_hits: usize,

        #[arg(long = "threshold", short = 'd', required = true)]
        threshold: f64,

        #[arg(long = "denominator", short = 'n')]
        denominator: Denominator,

        #[arg(long = "unique-weight", short = 'u', default_value = "0", help = "Weight for unique matches (in the range [0,1])")]
        unique_weight: f64,
    },

    #[command(arg_required_else_help = true, name = "segment-consensus-pseudoalign")]
    SegmentConsensusPseudoalign{
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,

        #[arg(long = "min-hits", short = 'm', required = true, help = "This is applied for each segment")]
        min_hits: usize,

        #[arg(long = "segment-length", required = true, default_value = "100")]
        segment_length: usize,

        #[arg(long = "min-unique-segments", required = true, default_value = "1")]
        min_unique_segments: usize,

        #[arg(long = "min-shared-segments", required = true, default_value = "5")]
        min_shared_segments: usize,

        #[arg(long = "consensus-threshold", short = 'd', required = true)]
        consensus_threshold: f64,
    },

    #[command(arg_required_else_help = true, name = "dump-pseudoalignment-data")]
    DumpPseudoalignmentData {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,
    
        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,
    
        #[arg(long = "min-hits", short = 'm', required = true)]
        min_hits: usize,
    },

    #[command(arg_required_else_help = true, name = "pseudoalign-into-EM")]
    PseudoalignIntoEM {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,
    
        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,
    
        #[arg(long = "min-hits", short = 'm', required = true)]
        min_hits: usize,
    
        #[arg(long = "n-threads", short = 't', default_value = "1")]
        n_threads: usize,
    
        #[arg(long = "max-iterations", default_value = "2000")]
        max_iterations: usize,
    
        #[arg(long = "numerator", default_value = "hits", 
            value_parser = PossibleValuesParser::new(["hits", "distinguishing"]))]
        numerator: String,
    
        #[arg(long = "denominator", default_value = "all", 
            value_parser = PossibleValuesParser::new(["all", "relevant", "max-distinguishing"]))]
        denominator: String,
    
        #[arg(long = "likelihood", default_value = "linear", 
            value_parser = PossibleValuesParser::new(["linear", "softmax", "99p"]))]
        likelihood: String,

    },

    #[command(arg_required_else_help = true, name = "intersection-pseudoalign-into-EM")]
    IntersectionPseudoalignIntoEM {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,
    
        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,
    
        #[arg(long = "min-hits", short = 'm', required = true)]
        min_hits: usize,
    
        #[arg(long = "n-threads", short = 't', default_value = "1")]
        n_threads: usize,
    
        #[arg(long = "max-iterations", default_value = "2000")]
        max_iterations: usize,
    
        #[arg(long = "initial-likelihood-ratio", short = 'w', default_value = "99")]
        initial_likelihood_ratio: f64,
    
        #[arg(long = "static-likelihood", action = clap::ArgAction::SetTrue)]
        static_likelihood: bool,
    },

    #[command(arg_required_else_help = true, name = "print-color-sets")]
    PrintColorSets {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,
    
        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,

        #[arg(long = "print-kmers", short = 'p', help = "Also print the k-mers on each line")]
        print_kmers: bool,
    },

    #[command(arg_required_else_help = true, name = "compress-colors")]
    CompressColors {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "sample-distance", short = 'd', default_value = "1")]
        sample_distance: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "validation-queries", help = "For debugging: checks that the compressed colors match the original colors for all k-mers in the input file.")]
        validation_queries: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, name = "merge-compressed-indexes")]
    MergeCompressedIndexes {
        #[arg(long = "index1", required = true)]
        index1_file: PathBuf,

        #[arg(long = "index1-from-sbwt", help = "Load index 1 from format written by sbwt-rs-cli, and color with a single color")]
        index1_from_sbwt: bool,

        #[arg(long = "index2", required = true)]
        index2_file: PathBuf,

        #[arg(long = "index2-from-sbwt", help = "Load index 2 from format written by sbwt-rs-cli, and color with a single color")]
        index2_from_sbwt: bool,

        #[arg(long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,
    },

}


fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();
    match args.command {
        Subcommands::Build { input: input_fof, output: out_path, temp_dir, k, n_threads } => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            let index = ColoredKmers::new_from_files(input_paths.as_slice(), k, n_threads, &temp_dir);
            index.serialize(&mut BufWriter::new(File::create(&out_path).unwrap()));
        },
        Subcommands::Import { sbwt_ascii_dump, color_dump_prefix, out: out_path } => {
            let unitig_filename = format!("{}.unitigs.fa", color_dump_prefix.to_str().unwrap());
            let color_sets_filename = format!("{}.color_sets.txt", color_dump_prefix.to_str().unwrap());
            let metadata_filename = format!("{}.metadata.txt", color_dump_prefix.to_str().unwrap());

            let sbwt_in = BufReader::new(File::open(sbwt_ascii_dump).unwrap());
            let unitigs_in = BufReader::new(File::open(unitig_filename).unwrap());
            let color_sets_in = BufReader::new(File::open(color_sets_filename).unwrap());
            let metadata_in = BufReader::new(File::open(metadata_filename).unwrap());

            let mut out = BufWriter::new(File::create(out_path).unwrap());

            let index = colored_kmers::ColoredKmers::new_from_new_themisto_index_dump(sbwt_in, metadata_in, unitigs_in, color_sets_in, 0);
            index.serialize(&mut out);
        },
        Subcommands::IntersectionPseudoalign { index: index_path, query: query_path, min_hits} => {
            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

            log::info!("Pseudoaligning (intersection method)");
            while let Some(rec) = reader.read_next().unwrap(){
                let intersection = index.intersection_pseudoalignment(rec.seq, min_hits);
                println!("{:?}", intersection.iter_ones().collect::<Vec::<usize>>()); // Todo: print indices of non-zero
            }
        },
        Subcommands::ThresholdPseudoalign { index: index_path, query: query_path, min_hits, threshold, denominator, unique_weight} => {
            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

            log::info!("Pseudoaligning (threshold method, denominator = {:?})", denominator);
            while let Some(rec) = reader.read_next().unwrap(){
                let pa_data = index.compute_pseudoalignment_data(rec.seq, 0);
                let den = match denominator {
                    Denominator::All => pa_data.n_all_kmers as f64,
                    Denominator::Relevant => pa_data.n_unique_kmers as f64 * unique_weight + pa_data.n_relevant_kmers as f64 * (1.0 - unique_weight),
                    Denominator::MaxHits => {
                        (0..index.n_colors()).map(|i| 
                            pa_data.unique_hit_counts[i] as f64 * unique_weight +
                            pa_data.hit_counts[i] as f64 * (1.0 - unique_weight)).
                            max_by(|a,b| a.partial_cmp(b).unwrap())
                            .expect("Programming error: no hit counts found")
                    }
                };
                let compatible_colors = unique_support_combination_method(&pa_data.unique_hit_counts, &pa_data.hit_counts, unique_weight, 0, min_hits, den, threshold);
                println!("{:?}", compatible_colors);
            }

        },
        Subcommands::SegmentConsensusPseudoalign { index: index_path, query: query_path, min_hits, segment_length, min_shared_segments, min_unique_segments, consensus_threshold } => {
            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

            log::info!("Pseudoaligning (segment consensus method)");
            while let Some(rec) = reader.read_next().unwrap(){
                let color_sets: Vec<Vec<usize>> = rec.seq.windows(segment_length).step_by(segment_length).map(|segment| {
                    let bitmap = index.intersection_pseudoalignment(segment, min_hits);
                    bitmap.iter_ones().collect::<Vec::<usize>>()
                }).collect();

                let slices: Vec<&[usize]> = color_sets.iter().map(|v| v.as_slice()).collect();
                let consensus = compatibility_criteria::resolve_consensus_compatibility(&slices, index.n_colors(), min_unique_segments, min_shared_segments, consensus_threshold);
                println!("{:?}", consensus);
            }
        }
        Subcommands::DumpPseudoalignmentData { index: index_path, query: query_path, min_hits } => {
            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

            log::info!("Pseudoaligning");
            while let Some(rec) = reader.read_next().unwrap(){
                let data = index.compute_pseudoalignment_data(rec.seq, min_hits);
                let json = serde_json::to_string(&data).unwrap();
                println!("{}", json);
            }
        },
        Subcommands::PseudoalignIntoEM { 
                index: index_path, 
                query: query_path, 
                min_hits, 
                n_threads, 
                max_iterations, 
                numerator, 
                denominator, 
                likelihood: likelihood_type } => {

            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            log::info!("Loaded index with {} distinct k-mers and {} colors", index.n_kmers(), index.n_colors());

            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

            #[allow(clippy::type_complexity)] // Is fine, stop complaining
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

        },
        Subcommands::IntersectionPseudoalignIntoEM { index: index_path, query: query_path, min_hits, n_threads, max_iterations, initial_likelihood_ratio: initial_w, static_likelihood } => {
            let optimize_w = !static_likelihood;
            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

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
        },
        Subcommands::PrintColorSets { index: index_path, query: query_path, print_kmers } => {
            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();
            while let Some(rec) = reader.read_next().unwrap() {
                println!(">{}", String::from_utf8(rec.head.to_vec()).unwrap());
                let sets = index.get_all_color_sets(rec.seq);
                for (set_idx, set) in sets.iter().enumerate() {
                    if print_kmers {
                        print!("{} ", String::from_utf8(rec.seq[set_idx..set_idx+index.get_k()].to_vec()).unwrap())
                    }
                    set.to_string();
                    let bitstring: String = set.iter().by_vals().map(|b| if b { '1' } else { '0' }).collect();
                    println!("{}", bitstring);
                }
            }
        },
        Subcommands::CompressColors { index: index_path, sample_distance, validation_queries, n_threads, outfile} => {
            let mut out = BufWriter::new(File::create(&outfile).unwrap()); // Open early to fail early if there is a problem

            log::info!("Loading index");
            let index = colored_kmers::ColoredKmers::load(&mut BufReader::new(File::open(index_path).unwrap()));
            log::info!("Compressing colors");

            if let Some(validation_queries) = validation_queries {
                // Clone the original index before construction to be able to compare to it
                let compressed = index.clone().compress_colors(sample_distance, n_threads);

                log::info!("Serializing to {}", outfile.display());
                compressed.serialize(&mut out); // Colors

                log::info!("Validating compressed colors for {}", validation_queries.display());
                let mut reader = jseqio::reader::DynamicFastXReader::from_file(&validation_queries).unwrap();
                while let Some(rec) = reader.read_next().unwrap() {
                    for kmer in rec.seq.windows(index.get_k()) {
                        let old_set: Vec<usize> = index.get_color_set(kmer).iter_ones().collect();

                        let colex_range = index.sbwt().search(kmer);
                        let mut new_set = Vec::<usize>::new();
                        if let Some(colex_range) = colex_range {
                            assert!(colex_range.len() == 1);
                            compressed.colex_to_set(colex_range.start).extract_and_push_colors_to(&mut new_set)
                        }

                        assert_eq!(old_set, new_set);
                    }
                }
                log::info!("All sets match");
            } else {
                log::info!("Copying SBWT to {}", outfile.display());

                let compressed = index.compress_colors(sample_distance, n_threads);
                log::info!("Serializing colors to {}", outfile.display());
                compressed.serialize(&mut out); // Colors
                log::info!("Finished");
            }
        },
        Subcommands::MergeCompressedIndexes{ index1_file, index2_file, index1_from_sbwt, index2_from_sbwt, n_threads, outfile} => {
            let mut out = BufWriter::new(File::create(&outfile).unwrap()); // Open early to fail early if there is a problem

            log::info!("Loading index 1");
            let mut in1 = BufReader::new(File::open(index1_file).unwrap());
            let colors1 = if index1_from_sbwt {
                let SbwtIndexVariant::SubsetMatrix(mut sbwt1) = sbwt::load_sbwt_index_variant(&mut in1).unwrap();
                sbwt1.build_select(); // Required for merge
                let lcs1 = LcsArray::from_sbwt(&sbwt1, n_threads);
                compact_colored_kmers::CompactColexColoring::new_single_colored(Arc::new(sbwt1), lcs1, 10, n_threads) // Todo: 10 to CLI
            } else {
                compact_colored_kmers::CompactColexColoring::load(&mut in1, true) // Select support is required for merge
            };

            log::info!("Loading index 2");
            let mut in2 = BufReader::new(File::open(index2_file).unwrap());
            let colors2 = if index2_from_sbwt {
                let SbwtIndexVariant::SubsetMatrix(mut sbwt2) = sbwt::load_sbwt_index_variant(&mut in2).unwrap();
                sbwt2.build_select(); // Required for merge
                let lcs2 = LcsArray::from_sbwt(&sbwt2, n_threads);
                compact_colored_kmers::CompactColexColoring::new_single_colored(Arc::new(sbwt2), lcs2, 10, n_threads) // Todo: 10 to CLI
            } else {
                compact_colored_kmers::CompactColexColoring::load(&mut in1, true) // Select support is required for merge
            };


            log::info!("Merging");
            let merged_colored_kmers = compact_colored_kmers::merge_compact_colorings(colors1, colors2, true, n_threads);

            log::info!("Serializing");
            merged_colored_kmers.serialize(&mut out);

            log::info!("Finished");

        }
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
