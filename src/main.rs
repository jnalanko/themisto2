#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, path::{Path, PathBuf}, sync::Arc};
use bitmap_storage::BitmapStorage;
use clap::{Parser, Subcommand};
use colex_colored_kmers::CompactColexColoring;
use coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView};
use sbwt::{BitPackedKmerSortingDisk, LcsArray, SubsetMatrix};
use sparse_dense_storage::SparseDenseStorage;

mod EM;
mod bitmap_storage;
//mod themisto1_compatibility; // Out of date after updating sbwt-rs
mod compatibility_criteria;
mod colex_colored_kmers;
mod coloring_interface;
mod sparse_dense_storage;
mod io;

#[derive(Parser)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommands,
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Denominator { // Options for the CLI
    All,
    Relevant,
    MaxHits,
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum ColoringType{ // Options for the CLI
    Bitmaps,
    SparseDense,
}

impl ColoringType {
    pub fn serialization_id(&self) -> [u8; 8] {
        match self {
            ColoringType::Bitmaps => [0, 0, 0, 0, 0, 0, 0, 1],
            ColoringType::SparseDense => [0, 0, 0, 0, 0, 0, 0, 2],
        }
    }
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

        #[arg(help = "Directory for temporary files", long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(short, required = true)]
        k: usize,

        #[arg(long = "sample-distance", short = 'd', default_value = "1")]
        sample_distance: usize,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(long = "index-type")]
        index_type: ColoringType,

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

    #[command(arg_required_else_help = true, name = "merge-compressed-indexes")]
    MergeCompressedIndexes {
        #[arg(long = "index-file-list", required = true)]
        index_file_list: PathBuf,

        #[arg(long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "low-ram-mode", help = "Use more slower but more compact algorithm invert and merge SBWTs")]
        low_ram_mode: bool,
    },

}

fn build_coloring<CSS: ColorSetStorage>(
    sbwt: Arc<sbwt::SbwtIndex<SubsetMatrix>>, lcs: LcsArray, input_paths: &[PathBuf], n_threads: usize, sample_distance: usize) -> CompactColexColoring<CSS> {

    log::info!("Building uncompressed color bitmap");
    let color_storage = bitmap_storage::build_from_files(input_paths, &sbwt, &lcs, n_threads);
    log::info!("Compressing sets with unitig sampling distance {}", sample_distance);
    CompactColexColoring::<CSS>::new(sbwt, lcs, &color_storage.bitmap, color_storage.n_colors, sample_distance, n_threads)
}

#[allow(clippy::large_enum_variant)] // It's saying that it's almost a kilobyte. I don't understand why but ok.
enum IndexVariant {
    BitmapIndex(CompactColexColoring<BitmapStorage>),
    SparseDenseIndex(CompactColexColoring<SparseDenseStorage>),
}

fn load_index_variant(path: &Path, build_select: bool) -> IndexVariant {
    let mut input = BufReader::new(File::open(path).unwrap());
    let mut id_buf = [0u8; 8];
    input.read_exact(&mut id_buf).unwrap();
    if id_buf == ColoringType::Bitmaps.serialization_id() {
        let index = CompactColexColoring::<BitmapStorage>::load(&mut input, build_select);
        IndexVariant::BitmapIndex(index)
    } else if id_buf == ColoringType::SparseDense.serialization_id() {
        let index = CompactColexColoring::<SparseDenseStorage>::load(&mut input, build_select);
        IndexVariant::SparseDenseIndex(index)
    } else {
        panic!("Unrecognized index serialization ID: {:?}", id_buf);
    }
}

fn write_index_variant(index: &IndexVariant, out: &mut impl Write) {
    match index {
        IndexVariant::BitmapIndex(idx) => {
            out.write_all(&ColoringType::Bitmaps.serialization_id()).unwrap();
            idx.serialize(out);
        },
        IndexVariant::SparseDenseIndex(idx) => {
            out.write_all(&ColoringType::SparseDense.serialization_id()).unwrap();
            idx.serialize(out);
        },
    }
}


fn print_color_sets<CSS: ColorSetStorage>(index: &CompactColexColoring<CSS>, query_path: &Path, print_kmers: bool) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();
    // Buffered writing to stdout
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout);

    log::info!("Retrieving color sets for query sequences in {}", query_path.display());

    while let Some(rec) = reader.read_next().unwrap() {
        let sets = index.lookup_kmer_color_sets(rec.seq);
        for (set_idx, set) in sets.iter().enumerate() {
            if print_kmers {
                write!(out, "{} ", String::from_utf8(rec.seq[set_idx..set_idx+index.get_k()].to_vec()).unwrap()).unwrap();
            }
            if let Some(set) = set {
                // Print the set as space-separated list of color IDs
                set.iter().enumerate().map(|(i, id)| (i,id.to_string())).for_each(|(i,s)| {
                    if i == 0 {
                        write!(out,"{}", s).unwrap();
                    } else {
                        write!(out," {}", s).unwrap();
                    }
                }); 
            }
            writeln!(out).unwrap();
        }
    }
}

#[allow(clippy::manual_flatten)]
fn intersection_pseudoalignment<CSS: ColorSetStorage>(index: &CompactColexColoring<CSS>, query_path: &Path, min_hits: usize) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();
    // Buffered writing to stdout
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout);
    let mut query_idx = 0_usize;
    log::info!("Performing intersection pseudoalignment for query sequences in {}", query_path.display());
    while let Some(rec) = reader.read_next().unwrap(){
        let mut intersection = index.get_set_storage().get_full_set();
        let mut n_hits = 0_usize;
        for set in index.lookup_kmer_color_sets(rec.seq) {
            if let Some(set) = set {
                index.get_set_storage().intersect(&mut intersection, &set);
                n_hits += 1;
            }
        }

        // Write output
        write!(out, "{}", query_idx).unwrap();
        if n_hits >= min_hits {
            for color in intersection.iter() {
                write!(out, " {}", color).unwrap();
            }
        }
        writeln!(out).unwrap();

        query_idx += 1;

    }
}

#[allow(clippy::manual_flatten, clippy::len_zero)]
fn threshold_pseudoalignment<CSS: ColorSetStorage>(index: &CompactColexColoring<CSS>, query_path: &Path, min_hits: usize, threshold: f64, denominator: Denominator) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();
    // Buffered writing to stdout
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout);
    let mut query_idx = 0_usize;
    let n_colors = index.get_set_storage().get_full_set().iter().count(); // Todo len() for owned set
    let mut hit_counts = vec![0usize; n_colors];
    let mut nonzero_count_indices = vec![];
    log::info!("Performing threshold pseudoalignment for query sequences in {}", query_path.display());
    while let Some(rec) = reader.read_next().unwrap(){
        let mut n_relevant = 0_usize;
        let mut n_all = 0_usize;
        for set in index.lookup_kmer_color_sets(rec.seq) {
            if let Some(set) = set {
                for color in set.iter() {
                    hit_counts[color] += 1;
                    if hit_counts[color] == 1 {
                        nonzero_count_indices.push(color);
                    }
                }
                n_relevant += 1;
            }
            n_all += 1;
        }
        nonzero_count_indices.sort_unstable(); // Sort to output in sorted order by colors 

        // Write to output all colors that pass the threshold
        write!(out, "{}", query_idx).unwrap();
        if n_relevant >= min_hits && nonzero_count_indices.len() > 0 {
            let den = match denominator {
                Denominator::All => n_all as f64,
                Denominator::Relevant => n_relevant as f64,
                Denominator::MaxHits => {
                    let maxhits = nonzero_count_indices.iter().map(|color| hit_counts[color]).max();
                    maxhits.unwrap() as f64 // Safe because here nonzero_count_indices.len() > 0
                },
            };
            for color in nonzero_count_indices.iter() {
                if hit_counts[color] as f64 / den >= threshold {
                    write!(out, " {}", color).unwrap();
                }
            }
        }
        writeln!(out).unwrap();

        // Clean up
        for &color in &nonzero_count_indices {
            hit_counts[color] = 0;
        }
        nonzero_count_indices.clear();

        query_idx += 1;

    }
}

fn run_merge_tree(infiles: &[PathBuf], temp_dir: &Path, outfile: &Path, n_threads: usize, low_ram_mode: bool) {
    let n_rounds = infiles.len().div_ceil(2);
    let mut current_files: Vec<PathBuf> = infiles.to_vec();
    for round in 0..n_rounds {
        log::info!("Merge round {}", round);
        let mut next_files: Vec<PathBuf> = Vec::new();
        for pair in current_files.chunks(2) {
            if pair.len() == 2 {
                let outpath = if round == n_rounds - 1 {
                    outfile.to_path_buf() // Final output file
                } else {
                    temp_dir.join(format!("merge_round{}_{}.thm2", round, next_files.len()))
                };
                log::info!("Merging {} and {} into {}", pair[0].display(), pair[1].display(), outpath.display());
                let mut out = BufWriter::new(File::create(&outpath).unwrap());
                let colors1 = load_index_variant(&pair[0], true); // Select support is required
                let colors2 = load_index_variant(&pair[1], true); // Select support is required

                match (colors1, colors2) {
                    (IndexVariant::BitmapIndex(c1), IndexVariant::BitmapIndex(c2)) => {
                        log::info!("Merging bitmap indexes");
                        let merged_colored_kmers = colex_colored_kmers::merge_compact_colorings(c1, c2, low_ram_mode, n_threads);
                        log::info!("Serializing merged index to {}", outpath.display());
                        write_index_variant(&IndexVariant::BitmapIndex(merged_colored_kmers), &mut out);
                    },
                    (IndexVariant::SparseDenseIndex(c1), IndexVariant::SparseDenseIndex(c2)) => {
                        log::info!("Merging sparse-dense indexes");
                        let merged_colored_kmers = colex_colored_kmers::merge_compact_colorings(c1, c2, low_ram_mode, n_threads);
                        log::info!("Serializing merged index to {}", outpath.display());
                        write_index_variant(&IndexVariant::SparseDenseIndex(merged_colored_kmers), &mut out);
                    },
                    (IndexVariant::SparseDenseIndex(_), IndexVariant::BitmapIndex(_)) => {
                        panic!("Mismatched index types when merging: {} and {}", pair[0].display(), pair[1].display());
                    }
                    (IndexVariant::BitmapIndex(_), IndexVariant::SparseDenseIndex(_)) => {
                        panic!("Mismatched index types when merging: {} and {}", pair[0].display(), pair[1].display());
                    }
                }
                next_files.push(outpath);
            } else {
                next_files.push(pair[0].clone());
            }
        }
        current_files = next_files;
    }
}

fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { input: input_fof, output, temp_dir, k, n_threads, index_type, sample_distance} => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            let input_stream = io::ChainedInputStream::new(input_paths.clone());
            let mut out = BufWriter::new(File::create(&output).unwrap());

            let (mut sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(true)
                .k(k)
                .build_lcs(true)
                .n_threads(n_threads)
                .precalc_length(8)
                .algorithm(BitPackedKmerSortingDisk::new().dedup_batches(true).temp_dir(&temp_dir))
            .run(input_stream);
            log::info!("Building SBWT select support");
            sbwt.build_select();
            let sbwt = Arc::new(sbwt);
            let lcs = lcs.unwrap(); // Ok because we used .build_lcs(true)

            match index_type {
                ColoringType::Bitmaps => {
                    let index = build_coloring::<BitmapStorage>(sbwt, lcs, &input_paths, n_threads, sample_distance);
                    log::info!("Serializing bitmap index to {}", output.display());
                    write_index_variant(&IndexVariant::BitmapIndex(index), &mut out);
                },
                ColoringType::SparseDense => {
                    let index = build_coloring::<SparseDenseStorage>(sbwt, lcs, &input_paths, n_threads, sample_distance);
                    log::info!("Serializing sparse-dense index to {}", output.display());
                    write_index_variant(&IndexVariant::SparseDenseIndex(index), &mut out);
                },
            }

        },
        Subcommands::IntersectionPseudoalign { index: index_path, query: query_path, min_hits } => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            match index {
                IndexVariant::BitmapIndex(idx) => intersection_pseudoalignment(&idx, &query_path, min_hits),
                IndexVariant::SparseDenseIndex(idx) => intersection_pseudoalignment(&idx, &query_path, min_hits),
            };

        },
        Subcommands::ThresholdPseudoalign { index: index_path, query: query_path, min_hits, threshold, denominator} => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            match index {
                IndexVariant::BitmapIndex(idx) => threshold_pseudoalignment(&idx, &query_path, min_hits, threshold, denominator),
                IndexVariant::SparseDenseIndex(idx) => threshold_pseudoalignment(&idx, &query_path, min_hits, threshold, denominator),
            };
        },
        Subcommands::PrintColorSets { index: index_path, query: query_path, print_kmers } => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            match index {
                IndexVariant::BitmapIndex(idx) => print_color_sets(&idx, &query_path, print_kmers),
                IndexVariant::SparseDenseIndex(idx) => print_color_sets(&idx, &query_path, print_kmers),
            };
        },
        Subcommands::MergeCompressedIndexes { index_file_list, temp_dir, outfile, n_threads, low_ram_mode } => {
            let infiles: Vec<PathBuf> = BufReader::new(File::open(index_file_list).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            run_merge_tree(&infiles, &temp_dir, &outfile, n_threads, low_ram_mode);
        },
    }

/* 
    match args.command {
        Subcommands::Build { input: input_fof, output: out_path, temp_dir, k, n_threads } => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            let index = ColoredKmers::new_from_files(input_paths.as_slice(), k, n_threads, &temp_dir);
            index.serialize(&mut BufWriter::new(File::create(&out_path).unwrap()));
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
                let compressed = index.clone().compress_colors::<ColorSets>(sample_distance, n_threads);

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

                let compressed = index.compress_colors::<ColorSets>(sample_distance, n_threads);
                log::info!("Serializing colors to {}", outfile.display());
                compressed.serialize(&mut out); // Colors
                log::info!("Finished");
            }
        },
        Subcommands::MergeCompressedIndexes{ index_file_list, n_threads, outfile, temp_dir, low_ram_mode} => {
            let infiles: Vec<PathBuf> = BufReader::new(File::open(index_file_list).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            run_merge_tree(&infiles, &temp_dir, &outfile, n_threads, low_ram_mode);
        },
        Subcommands::BuildFromSbwt{ sbwt_file, outfile, n_threads, sample_distance} => {
            let mut out = BufWriter::new(File::create(&outfile).unwrap()); // Open early to fail early if there is a problem
            log::info!("Loading SBWT");
            let mut input = BufReader::new(File::open(sbwt_file).unwrap());
            let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt::load_sbwt_index_variant(&mut input).unwrap();
            log::info!("Building LCS array");
            let lcs = LcsArray::from_sbwt(&sbwt, n_threads);
            let index = compact_colored_kmers::CompactColexColoring::<ColorSets>::new_single_colored(Arc::new(sbwt), lcs, sample_distance, n_threads);
            index.serialize(&mut out);
        }
    } 
*/
}

/*


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

*/