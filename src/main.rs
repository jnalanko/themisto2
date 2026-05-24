#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material
#![allow(clippy::len_zero)] // !is_empty reads as "not is empty" which is not English 
#![allow(clippy::manual_is_multiple_of)] // Oh please

// We assume that usize is 64 bits in many places
#[cfg(not(target_pointer_width = "64"))]
compile_error!("This crate requires a 64-bit usize (target_pointer_width = 64).");

use std::{cmp::min, fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, ops::Range, path::{Path, PathBuf}, process::ExitCode, time::Instant};
use clap::{Parser, Subcommand};
use colex_colored_kmers::CompactColexKmers;
use coloring_interface::{ColorSetStorage, ColorSetView};
use io::{ChainedInputStreamWithRevComp, RewindableSeqStreamGenerator};
use jseqio::record::Record;
use parallel_ms_iteration::{DeduplicatingColorElementGenerator};
use sbwt::{BitPackedKmerSortingDisk, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};
use sparse_dense_storage::SparseDenseStorage;

use crate::{colex_colored_kmers::{ColexToColorSetMap, mark_key_kmers}, io::ChainedInputStream, parallel_ms_iteration::MsElementGenerator, report::report};

mod EM;
mod bitmap_storage;
mod index_import;
mod compatibility_criteria;
mod colex_colored_kmers;
mod coloring_interface;
mod sparse_dense_storage;
mod io;
mod build_from_ggcat;
mod set_of_sets_construction;
mod iterators;
mod parallel_ms_iteration;
mod atomic_bitmap;
mod int_vec;
mod finimizers;
mod util;
mod merge;
mod pseudoalignment;
mod pseudoalignment_metrics;
mod sparse_dense_storage_to_disk;
mod report;

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
pub enum ColoringType{ // Options for the CLI. This is old code: now only SparseDense is supported.
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

#[derive(Subcommand)]
pub enum Subcommands {
    #[command(arg_required_else_help = true)]
    Build {
        #[arg(help = "A file with one fasta/fastq filename per line", long = "file-colors")]
        color_fof: Option<PathBuf>,

        #[arg(help = "A fasta/fastq file, with one color per sequence", long = "seq-colors")]
        sequence_colors_file: Option<PathBuf>,

        #[arg(help = "Precomputed bit matrix SBWT (optional)", long = "sbwt", short = 's')]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Precomputed LCS array for the SBWT (optional)", long = "lcs", short = 'l')]
        lcs_path: Option<PathBuf>,

        #[arg(help = "Output filename. If --to-disk-in-pieces is given, this is interpreted as a filename prefix.", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Optional: Build from unitigs (requires odd k). This makes the construction much faster because now we can exploit the fact that the k-mers have already been deduplicated in the unitigs", short, long)]
        from_unitigs: bool,

        #[arg(help = "Directory for temporary files", long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(short, required = true)]
        k: usize,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(long = "to-disk-in-pieces", help = "Build the final sets in pieces, flushing to disk after each piece. This is useful if there is not enough memory to hold the final data structure in memory at once.")]
        n_pieces: Option<usize>,

        // Hidden because index import and export assume that reverse complements are indexed.
        // The main purpose for this flag currently is to build the index for
        // the running example in the manuscript.
        #[arg(long = "forward-only", hide = true, default_value = "false")]
        forward_only: bool,

    },

    #[command(arg_required_else_help = true, name = "intersection-pseudoalign")]
    IntersectionPseudoalign {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,

        #[arg(long = "min-hits", short = 'm', default_value = "1")]
        min_hits: usize,

        #[arg(long = "n-threads", short = 't', required = true, default_value = "4")]
        n_threads: usize,

        #[arg(long = "output", short = 'o', help = "Output file. If omitted, writes to stdout.")]
        output: Option<PathBuf>,

        #[arg(long = "report-hit-counts", default_value = "false")]
        report_hit_counts: bool,

        #[arg(long = "sort-output", default_value = "false", help = "Write query results in the same order as the input sequences. Without this flag, results may be written in any order.")]
        sort_output: bool,

        #[arg(long = "themisto1-output-format", default_value = "false", help = "Write output in the Themisto 1 format (query index followed by space-separated colors) instead of JSONL. Does not support reporting pseudoalignment metrics.")]
        themisto1_output_format: bool,

        // Hidden option
        #[arg(long = "report-bases-covered", default_value = "false", hide = true)]
        report_bases_covered: bool,
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

        #[arg(long = "n-threads", short = 't', required = true, default_value = "4")]
        n_threads: usize,

        #[arg(long = "output", short = 'o', help = "Output file. If omitted, writes to stdout.")]
        output: Option<PathBuf>,

        #[arg(long = "report-hit-counts", default_value = "false")]
        report_hit_counts: bool,

        #[arg(long = "sort-output", default_value = "false", help = "Write query results in the same order as the input sequences. Without this flag, results may be written in any order.")]
        sort_output: bool,

        #[arg(long = "themisto1-output-format", default_value = "false", help = "Write output in the Themisto 1 format (query index followed by space-separated colors) instead of JSONL. Does not support reporting pseudoalignment metrics.")]
        themisto1_output_format: bool,

        // Hidden option
        #[arg(long = "report-bases-covered", default_value = "false", hide = true)]
        report_bases_covered: bool,
    },

    #[command(arg_required_else_help = true, name = "print-color-sets")]
    PrintColorSets {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', required = true)]
        query: PathBuf,

        #[arg(long = "print-kmers", short = 'p', help = "Also print the k-mers on each line")]
        print_kmers: bool,

        #[arg(long = "output", short = 'o', help = "Output file. If omitted, writes to stdout.")]
        output: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, name = "dump-color-names")]
    DumpColorNames{
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "output", short = 'o', help = "Output file. If omitted, writes to stdout.")]
        output: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, name = "merge")]
    Merge {
        #[arg(long = "index-file-list", required = true)]
        index_file_list: PathBuf,

        #[arg(long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "low-ram-mode", help = "Use more slower but more compact algorithm invert and merge SBWTs")]
        low_ram_mode: bool,
    },

    #[command(arg_required_else_help = true)]
    Import {
        #[arg(help = "Precomputed bit matrix SBWT (optional)", long = "sbwt", short = 's')]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Precomputed LCS array for the SBWT (optional)", long = "lcs", short = 'l')]
        lcs_path: Option<PathBuf>,

        #[arg(help = "Index text dump file prefix, as written by Fulgor 4.0.0", long = "color-dump-prefix", short = 'c', required = true)]
        color_dump_prefix: PathBuf,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(help = "Index output file", long = "out", short = 'o', required = true)]
        out: PathBuf,
    },

    #[command(arg_required_else_help = true)]
    Export {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(help = "Output file prefix", long = "output-prefix", short = 'o', required = true)]
        color_dump_prefix: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,
    },
    #[command(arg_required_else_help = true, name = "stats")]
    Stats {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,
    },

    #[command(name = "report")]
    Report {
        #[arg(help = "Themisto 2 index, used for color names", long = "index", short = 'i', required = true, value_parser = clap::value_parser!(PathBuf))]
        index: PathBuf,

        #[arg(help = "Pseudoalignment JSON file (one record per line). If omitted, reads from stdin.", long = "pseudoalignment-file", short = 'p', value_parser = clap::value_parser!(PathBuf))]
        input: Option<PathBuf>,

        #[arg(help = "Output file. If omitted, writes to stdout.", long = "output", short = 'o', value_parser = clap::value_parser!(PathBuf))]
        output: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, hide = true)]
    FinimizerStats {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,
    },
    #[command(arg_required_else_help = true, hide = true)]
    MinimizerStats {
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "minimizer-length", short = 'm', required = true)]
        m: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,
    },
}

enum Output {
    File(BufWriter<File>),
    Stdout(BufWriter<std::io::Stdout>),
}

fn open_output(path: Option<PathBuf>) -> Output {
    match path {
        Some(p) => Output::File(BufWriter::new(File::create(&p).unwrap())),
        None => Output::Stdout(BufWriter::new(std::io::stdout())),
    }
}

enum BuildMode {
    InMemory,
    ToDisk(PathBuf, usize), // Output path prefix, number of pieces
}

enum ColoredSeqInput {
    FileColors(Vec<PathBuf>), // One file per color
    SequenceColors(PathBuf), // One file, containing one sequence per color.
}

impl ColoredSeqInput {
    fn get_generator(&self) -> Box<dyn RewindableSeqStreamGenerator + Sync + Send> {
        match self {
            ColoredSeqInput::FileColors(path_bufs) => {
                Box::new(io::SeqStreamGeneratorFromFiles::new(path_bufs.clone()))
            },
            ColoredSeqInput::SequenceColors(path_buf) => {
                Box::new(io::SeqStreamGeneratorFromSingleFile::new(path_buf.clone()))
            }
        }
    }
}

// Returns the index if BuildMode is InMemory, otherwise serializes as a set of files to disk.
#[allow(clippy::too_many_arguments)]
fn build_coloring<CSS: ColorSetStorage + Send>(sbwt: sbwt::SbwtIndex<SubsetMatrix>, lcs: LcsArray, input_mode: ColoredSeqInput, color_names: &[String], n_threads: usize, sample_distance: usize, from_unitigs: bool, forward_only: bool, build_mode: BuildMode) -> Option<CompactColexKmers<CSS>>{

    log::info!("Building distinct color set structure");
    let n_colors = color_names.len();
    if let ColoredSeqInput::FileColors(v) = &input_mode {
        assert!(v.len() == n_colors, "Number of color names does not match the number of input files");
    }

    let mut all_input_seq_files = Vec::<PathBuf>::new();
    match &input_mode {
        ColoredSeqInput::FileColors(v) => all_input_seq_files.extend(v.clone()),
        ColoredSeqInput::SequenceColors(file) => all_input_seq_files.push(file.clone()),
    }

    log::info!("=== PHASE 1/3: Marking key k-mers ===");
    let key_kmer_marks = if forward_only {
        let phase1_input_stream = ChainedInputStream::new(all_input_seq_files);
        mark_key_kmers(&sbwt, &lcs, sample_distance, phase1_input_stream, n_threads)
    } else {
        let phase1_input_stream = ChainedInputStreamWithRevComp::new(all_input_seq_files);
        mark_key_kmers(&sbwt, &lcs, sample_distance, phase1_input_stream, n_threads)
    };
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / sbwt.n_kmers() as f64 * 100.0);
    assert_eq!(key_kmer_marks.len(), sbwt.n_sets());

    log::info!("=== PHASE 2/3: Building color set finperprints for key k-mers ===");
    let color_stream_gen = input_mode.get_generator();
    let random_seed = 123123; // Todo: be more random
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id) = if from_unitigs {
        let ms_gen = MsElementGenerator::new(color_stream_gen, StreamingIndex::new(&sbwt, &lcs), !forward_only);
        set_of_sets_construction::find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(ms_gen, key_kmer_marks.clone(), n_colors, n_threads, random_seed)
    } else {
        let ms_gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, color_stream_gen, !forward_only);
        set_of_sets_construction::find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(ms_gen, key_kmer_marks.clone(), n_colors, n_threads, random_seed)
    };

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    match build_mode {
        BuildMode::InMemory => {
            let color_stream_gen = input_mode.get_generator();
            let css = if from_unitigs {
                let gen = MsElementGenerator::new(color_stream_gen, StreamingIndex::new(&sbwt, &lcs), !forward_only);
                set_of_sets_construction::build_color_set_storage(n_colors, repr_kmer_marks, distinct_set_sizes, gen, n_threads)
            } else {
                let gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, color_stream_gen, !forward_only);
                set_of_sets_construction::build_color_set_storage(n_colors, repr_kmer_marks, distinct_set_sizes, gen, n_threads)
            };

            log::info!("Building rank support for key k-mer marks");
            let mut key_kmer_marks = util::bitvec_to_simple_sds_bitvec(key_kmer_marks);
            key_kmer_marks.enable_rank();
            assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));
            let colex_map = ColexToColorSetMap {
                sampling: key_kmer_marks, 
                color_set_ids: key_kmer_idx_to_set_id,
            };

            let cck = CompactColexKmers::<CSS>::new(sbwt, lcs, colex_map, css, Some(&color_names));
            Some(cck) 
        },
        BuildMode::ToDisk(out_prefix, n_pieces) => {
            assert!(n_pieces != 0);
            let input_paths = match &input_mode {
                ColoredSeqInput::FileColors(path_bufs) => path_bufs,
                ColoredSeqInput::SequenceColors(_) => {
                    panic!("ToDisk mode with sequence colors is not yet supported");
                    // The issue is that the current code chunks the input files.
                    // It takes some work to translate this behaviour to a single
                    // input file.
                },
            };
            let chunk_size = n_colors.div_ceil(n_pieces);
            if from_unitigs {
                let mut gens = Vec::<(MsElementGenerator, Range::<usize>)>::new();
                for (chunk_id, chunk) in input_paths.chunks(chunk_size).enumerate() {
                    let color_id_range = chunk_id*chunk_size .. min((chunk_id+1)*chunk_size, n_colors);
                    let gen = Box::new(io::SeqStreamGeneratorFromFiles::new(chunk.to_owned()));
                    let ms_gen = MsElementGenerator::new(gen, StreamingIndex::new(&sbwt, &lcs), !forward_only);
                    gens.push((ms_gen, color_id_range));
                }
                set_of_sets_construction::build_color_set_storage_to_disk::<CSS>(repr_kmer_marks, distinct_set_sizes, gens, &out_prefix, n_threads);
            } else {
                let mut gens = Vec::<(DeduplicatingColorElementGenerator, Range::<usize>)>::new();
                for (chunk_id, chunk) in input_paths.chunks(chunk_size).enumerate() {
                    let color_id_range = chunk_id*chunk_size .. min((chunk_id+1)*chunk_size, n_colors);
                    let gen = Box::new(io::SeqStreamGeneratorFromFiles::new(chunk.to_owned()));
                    let ms_gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, gen, !forward_only);
                    gens.push((ms_gen, color_id_range));
                }
                set_of_sets_construction::build_color_set_storage_to_disk::<CSS>(repr_kmer_marks, distinct_set_sizes, gens, &out_prefix, n_threads);
            };
            None
        },
    }
}

// There used to be two variants. Now there is just one.
// This enum will serialize with an id specifying the variant.
// We'll keep this still in case we want to add support for 
// other variants in the future.
enum IndexVariant {
    SparseDenseIndex(CompactColexKmers<SparseDenseStorage>),
}

fn load_index_color_names_only(path: &Path) -> Vec<String> {
    let mut input = BufReader::new(File::open(path).unwrap());
    let mut id_buf = [0u8; 8];
    input.read_exact(&mut id_buf).unwrap();
    if id_buf == ColoringType::SparseDense.serialization_id() {
        CompactColexKmers::<SparseDenseStorage>::load_color_names_only(&mut input)
    } else {
        panic!("Unrecognized index serialization ID: {:?}", id_buf);
    }
}

fn load_index_variant(path: &Path, build_select: bool) -> IndexVariant {
    let mut input = BufReader::new(File::open(path).unwrap());
    let mut id_buf = [0u8; 8];
    input.read_exact(&mut id_buf).unwrap();
    if id_buf == ColoringType::SparseDense.serialization_id() {
        let index = CompactColexKmers::<SparseDenseStorage>::load(&mut input, build_select);
        IndexVariant::SparseDenseIndex(index)
    } else {
        panic!("Unrecognized index serialization ID: {:?}", id_buf);
    }
}

fn write_index_variant(index: &IndexVariant, out: &mut impl Write) {
    match index {
        IndexVariant::SparseDenseIndex(idx) => {
            out.write_all(&ColoringType::SparseDense.serialization_id()).unwrap();
            idx.serialize(out);
        },
    }
}

fn print_color_names_impl<W: Write>(names: &[String], out: &mut W) {
    for (id, name) in names.iter().enumerate() {
        writeln!(out, "{}\t{}", id, name).unwrap();
    }
}

fn print_color_names(names: &[String], out: Output) {
    match out {
        Output::File(mut w) => print_color_names_impl(names, &mut w),
        Output::Stdout(mut w) => print_color_names_impl(names, &mut w),
    }
}

fn print_stats<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, n_threads: usize) {
    let stats = index.compute_index_stats(n_threads);
    println!("Number of k-mers: {}", stats.n_kmers);
    println!("Number of colors: {}", stats.n_colors);
    println!("SBWT length: {}", stats.n_sbwt_sets);
    println!("Number of distinct color sets: {}", stats.n_distinct_color_sets);
    println!("Mean size of distinct color sets: {:.2}", stats.mean_size_of_distinct_sets());
    println!("Mean k-mer color set size: {:.2}", stats.mean_kmer_color_set_size());
    println!("Fraction of key k-mers: {:.2}%", stats.sample_fraction() * 100.0);
    println!("Number of forward unitigs (not bidirected): {}", stats.n_unitigs);
    println!("Min unitig length: {}", stats.min_unitig_length);
    println!("Max unitig length: {}", stats.max_unitig_length);
    println!("Mean unitig length: {:.2}", stats.mean_unitig_length());
}


fn print_color_sets<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, query_path: &Path, print_kmers: bool, out: Output) {
    match out {
        Output::File(w) => print_color_sets_impl(index, query_path, print_kmers, w),
        Output::Stdout(w) => print_color_sets_impl(index, query_path, print_kmers, w),
    }
}

fn print_color_sets_impl<CSS: ColorSetStorage, W: Write>(index: &CompactColexKmers<CSS>, query_path: &Path, print_kmers: bool, mut out: W) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();

    log::info!("Retrieving color sets for query sequences in {}", query_path.display());

    let mut total_lookup_nanoseconds = 0_u128;
    let mut total_io_nanoseconds = 0_u128;
    let mut total_bases_read = 0_usize; 
    let mut total_reads = 0_usize; 
    while let Some(rec) = reader.read_next().unwrap() {
        total_bases_read += rec.seq.len();
        total_reads += 1;
        let lookup_start = Instant::now();
        let sets = index.lookup_kmer_color_sets(rec.seq);
        total_lookup_nanoseconds += lookup_start.elapsed().as_nanos();

        let io_start = Instant::now();
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
        total_io_nanoseconds += io_start.elapsed().as_nanos();
    }
    log::info!("Processed total of {} bases in {} reads", total_bases_read, total_reads);
    log::info!("Total query time (excluding index loading): {:.2} s", (total_lookup_nanoseconds + total_io_nanoseconds) as f64 / 1_000_000_000.0);
    log::info!("Total time in color set lookup: {:.2} s", total_lookup_nanoseconds as f64 / 1_000_000_000.0);
    log::info!("Total time in output: {:.2} s", total_io_nanoseconds as f64 / 1_000_000_000.0);
    log::info!("Time in lookup per nucleotide: {:.2} ns", total_lookup_nanoseconds as f64 / (total_bases_read as f64));
}

#[allow(clippy::manual_flatten)]
fn intersection_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, query_path: &Path, min_hits: usize, metrics: &[pseudoalignment_metrics::Metric], n_threads: usize, out: Output, sort_output: bool, output_format: pseudoalignment::OutputFormat) {
    log::info!("Running intersection pseudoalignment for query sequences in {}", query_path.display());

    let create_new_aligner = move || {
        let aligner = pseudoalignment::IntersectionPseudoaligner::new(min_hits);
        Box::new(aligner) as Box<dyn pseudoalignment::Pseudoaligner<CSS> + Send>
    };

    match out {
        Output::File(w) => pseudoalignment::run_pseudoalignment(index, query_path, w, create_new_aligner, metrics, n_threads, sort_output, output_format),
        Output::Stdout(w) => pseudoalignment::run_pseudoalignment(index, query_path, w, create_new_aligner, metrics, n_threads, sort_output, output_format),
    }
    log::info!("Finished");
}

#[allow(clippy::manual_flatten, clippy::len_zero)]
fn threshold_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &CompactColexKmers<CSS>, query_path: &Path, min_hits: usize, threshold: f64, denominator: Denominator, metrics: &[pseudoalignment_metrics::Metric], n_threads: usize, out: Output, sort_output: bool, output_format: pseudoalignment::OutputFormat) {
    // Map from the CLI denominator enum to the pseudoalignment denominator enum.
    let denominator = match denominator {
        Denominator::All => pseudoalignment::Denominator::All,
        Denominator::Relevant => pseudoalignment::Denominator::Relevant,
        Denominator::MaxHits => pseudoalignment::Denominator::MaxHits,
    };

    let n_colors = index.get_set_storage().n_colors();

    log::info!("Running threshold pseudoalignment for query sequences in {}", query_path.display());
    let create_new_aligner = move || {
        let aligner = pseudoalignment::ThresholdPseudoaligner::new(
            n_colors,
            threshold,
            min_hits,
            denominator
        );
        Box::new(aligner) as Box<dyn pseudoalignment::Pseudoaligner<CSS> + Send>
    };

    match out {
        Output::File(w) => pseudoalignment::run_pseudoalignment(index, query_path, w, create_new_aligner, metrics, n_threads, sort_output, output_format),
        Output::Stdout(w) => pseudoalignment::run_pseudoalignment(index, query_path, w, create_new_aligner, metrics, n_threads, sort_output, output_format),
    }
    log::info!("Finished");
}

fn run_merge_tree(infiles: &[PathBuf], temp_dir: &Path, outfile: &Path, n_threads: usize, low_ram_mode: bool, sample_distance: usize) {
    let n_rounds = (infiles.len().next_power_of_two()).trailing_zeros() as usize;
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
                    (IndexVariant::SparseDenseIndex(c1), IndexVariant::SparseDenseIndex(c2)) => {
                        log::info!("Merging sparse-dense indexes");
                        let merged_colored_kmers = merge::merge_compact_colex_kmers(c1, c2, low_ram_mode, sample_distance, n_threads);
                        log::info!("Serializing merged index to {}", outpath.display());
                        write_index_variant(&IndexVariant::SparseDenseIndex(merged_colored_kmers), &mut out);
                    },
                }
                next_files.push(outpath);
            } else {
                next_files.push(pair[0].clone());
            }
        }
        current_files = next_files;
    }
}

fn export_index<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, out_prefix: &Path, n_threads: usize) {
    let out_prefix = out_prefix.as_os_str().to_str().unwrap().to_owned();

    let mut metadata_filename = out_prefix.clone();
    metadata_filename.push_str(".metadata.txt");

    let mut unitig_filename = out_prefix.clone();
    unitig_filename.push_str(".unitigs.fa");

    let mut colors_filename = out_prefix.clone();
    colors_filename.push_str(".color_sets.txt");

    let metadata_out = BufWriter::new(File::create(metadata_filename).unwrap());
    let unitigs_out = BufWriter::new(File::create(unitig_filename).unwrap());
    let colors_out = BufWriter::new(File::create(colors_filename).unwrap());

    index.export_colored_unitigs(metadata_out, unitigs_out, colors_out, n_threads); 
}

fn get_sbwt_and_lcs(sbwt_path: &Option<PathBuf>, lcs_path: &Option<PathBuf>, temp_dir: &Path, input_stream: io::ChainedInputStream, k: usize, forward_only: bool, n_threads: usize) -> (SbwtIndex<SubsetMatrix>, LcsArray) {
    let (sbwt, lcs) = if let Some(sbwt_path) = sbwt_path {
        log::info!("Loading SBWT from {}", sbwt_path.display());
        let mut sbwt_in = BufReader::new(File::open(sbwt_path).unwrap());
        let sbwt::SbwtIndexVariant::SubsetMatrix(mut sbwt) = sbwt::load_sbwt_index_variant(&mut sbwt_in).unwrap();

        assert_eq!(sbwt.k(), k);

        log::info!("Building select support for SBWT");
        sbwt.build_select();

        let lcs = if let Some(lcs_path) = lcs_path {
            log::info!("Loading the LCS array from {}", lcs_path.display());
            LcsArray::load(&mut BufReader::new(File::open(lcs_path).unwrap())).unwrap()
        } else {
            log::info!("Building LCS array");
            LcsArray::from_sbwt(&sbwt, n_threads, true)
        };
        (sbwt, lcs)
    } else {
        let (mut sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(!forward_only)
            .k(k)
            .build_lcs(true)
            .n_threads(n_threads)
            .precalc_length(8)
            .algorithm(BitPackedKmerSortingDisk::new().dedup_batches(true).temp_dir(temp_dir))
        .run(input_stream);
        log::info!("Building SBWT select support");
        sbwt.build_select();
        let sbwt = sbwt;
        let lcs = lcs.unwrap(); // Ok because we used .build_lcs(true)
        (sbwt, lcs)
    };

    (sbwt, lcs)

}

fn into_metric_list(report_hit_counts: bool, report_bases_covered: bool) -> Vec<crate::pseudoalignment_metrics::Metric> {
    let mut metrics: Vec<pseudoalignment_metrics::Metric> = vec![];
    if report_hit_counts {
        metrics.push(pseudoalignment_metrics::Metric::KmerHits);
    }
    if report_bases_covered {
        metrics.push(pseudoalignment_metrics::Metric::BasesCovered);
    }
    metrics
}

// Resolve the requested output format, rejecting the unsupported combination of
// Themisto 1 output with per-read metrics before any work starts.
fn resolve_output_format(themisto1_output_format: bool, metrics: &[pseudoalignment_metrics::Metric]) -> pseudoalignment::OutputFormat {
    if themisto1_output_format {
        if !metrics.is_empty() {
            log::error!("The Themisto 1 output format does not support reporting pseudoalignment metrics. Drop --themisto1-output-format or the metric reporting flags.");
            std::process::exit(1);
        }
        pseudoalignment::OutputFormat::Themisto1
    } else {
        pseudoalignment::OutputFormat::Jsonl
    }
}

fn main() -> std::process::ExitCode {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { color_fof, sequence_colors_file, output, temp_dir, k, n_threads, sample_distance, sbwt_path, lcs_path, from_unitigs, forward_only, n_pieces} => {
            if k % 2 == 0 && from_unitigs {
                panic!("--from_unitigs requires odd k");
            }

            let (input_mode, all_input_paths, color_names) = match (color_fof, sequence_colors_file) {
                (None, None) => panic!("Must give one of --file-colors or --seq-colors"),
                (Some(_), Some(_)) => todo!("Must not give both --file-colors and --seq-colors"),
                (Some(color_fof), None) => {
                    let input_paths: Vec<PathBuf> = BufReader::new(File::open(color_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();

                    // Use input paths also as color names
                    let color_names: Vec<String> = input_paths.iter().map(|p| p.clone().into_os_string().into_string().unwrap()).collect();
                    (ColoredSeqInput::FileColors(input_paths.clone()), input_paths, color_names)
                }
                (None, Some(sequence_colors_file)) => {
                    // Read color names from the sequence file
                    log::info!("Reading all sequence names from the input file");
                    let mut color_names = Vec::<String>::new();
                    // TODO: could use a faster parser here
                    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&sequence_colors_file).unwrap();
                    while let Some(rec) = reader.read_next().unwrap() {
                        color_names.push(String::from_utf8(rec.name().to_owned()).unwrap());
                    }
                    log::info!("Read {} sequences", color_names.len());
                    (ColoredSeqInput::SequenceColors(sequence_colors_file.clone()), vec![sequence_colors_file], color_names)
                },
            };

            let input_stream = io::ChainedInputStream::new(all_input_paths.clone());

            // Check that the output file can be created
            let _out_test = File::create(&output).unwrap();
            std::fs::remove_file(&output).unwrap();

            let (sbwt, lcs) = get_sbwt_and_lcs(&sbwt_path, &lcs_path, &temp_dir, input_stream, k, forward_only, n_threads);

            match n_pieces {
                None => {
                    let mut out = BufWriter::new(File::create(&output).unwrap());
                    let index = build_coloring::<SparseDenseStorage>(sbwt, lcs, input_mode, &color_names, n_threads, sample_distance, from_unitigs, forward_only, BuildMode::InMemory);
                    let index = index.unwrap(); // Ok because we passed in InMemory
                    log::info!("Serializing sparse-dense index to {}", output.display());
                    write_index_variant(&IndexVariant::SparseDenseIndex(index), &mut out);
                },
                Some(n_pieces) => {
                    build_coloring::<SparseDenseStorage>(sbwt, lcs, input_mode, &color_names, n_threads, sample_distance, from_unitigs, forward_only, BuildMode::ToDisk(output, n_pieces));
                }
            }
        },

        Subcommands::IntersectionPseudoalign { index: index_path, query: query_path, min_hits, n_threads, report_bases_covered, report_hit_counts, sort_output, themisto1_output_format, output} => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            let metrics = into_metric_list(report_hit_counts, report_bases_covered);
            let output_format = resolve_output_format(themisto1_output_format, &metrics);
            let out = open_output(output);
            match index {
                IndexVariant::SparseDenseIndex(idx) => intersection_pseudoalignment(&idx, &query_path, min_hits, &metrics, n_threads, out, sort_output, output_format),
            };

        },
        Subcommands::ThresholdPseudoalign { index: index_path, query: query_path, min_hits, threshold, denominator, n_threads, report_hit_counts, report_bases_covered, sort_output, themisto1_output_format, output} => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            let metrics = into_metric_list(report_hit_counts, report_bases_covered);
            let output_format = resolve_output_format(themisto1_output_format, &metrics);
            let out = open_output(output);
            match index {
                IndexVariant::SparseDenseIndex(idx) => threshold_pseudoalignment(&idx, &query_path, min_hits, threshold, denominator, &metrics, n_threads, out, sort_output, output_format),
            };
        },
        Subcommands::PrintColorSets { index: index_path, query: query_path, print_kmers, output } => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            let out = open_output(output);
            match index {
                IndexVariant::SparseDenseIndex(idx) => print_color_sets(&idx, &query_path, print_kmers, out),
            };
        },
        Subcommands::DumpColorNames{ index: index_path, output } => {
            log::info!("Loading color names");
            let color_names = load_index_color_names_only(&index_path);
            let out = open_output(output);
            log::info!("Dumping color names");
            print_color_names(&color_names, out);
        },
        Subcommands::Merge { index_file_list, temp_dir, outfile, n_threads, low_ram_mode, sample_distance } => {
            let infiles: Vec<PathBuf> = BufReader::new(File::open(index_file_list).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            run_merge_tree(&infiles, &temp_dir, &outfile, n_threads, low_ram_mode, sample_distance);
        },
        Subcommands::Import { sbwt_path, lcs_path, color_dump_prefix, out: out_path, n_threads, temp_dir, sample_distance} => {
            let unitig_filename = format!("{}.unitigs.fa", color_dump_prefix.to_str().unwrap());
            let color_sets_filename = format!("{}.color_sets.txt", color_dump_prefix.to_str().unwrap());
            let metadata_filename = format!("{}.metadata.txt", color_dump_prefix.to_str().unwrap());

            // Try to open to check that the files are found
            BufReader::new(File::open(&unitig_filename).unwrap());
            BufReader::new(File::open(&color_sets_filename).unwrap());
            BufReader::new(File::open(&metadata_filename).unwrap());

            log::info!("Reading metadata from {}", metadata_filename);
            let metadata = index_import::read_index_dump_metadata(BufReader::new(File::open(&metadata_filename).unwrap()));

            let input_stream = io::ChainedInputStream::new(vec![PathBuf::from(&unitig_filename)]) ;
            let (sbwt,lcs) = get_sbwt_and_lcs(&sbwt_path, &lcs_path, &temp_dir, input_stream, metadata.k, true, n_threads);

            if sbwt.k() != metadata.k {
                log::error!("SBWT k does not match the index dump k ({} vs {})", sbwt.k(), metadata.k);
                return ExitCode::FAILURE;
            }

            let mut out = BufWriter::new(File::create(&out_path).unwrap());

            let unitig_dump = BufReader::new(File::open(&unitig_filename).unwrap());
            let color_dump = BufReader::new(File::open(&color_sets_filename).unwrap());
            let metadata_dump = BufReader::new(File::open(&metadata_filename).unwrap());

            let index = CompactColexKmers::<SparseDenseStorage>::new_from_colored_unitig_dump(
                sbwt, lcs, sample_distance, n_threads, metadata_dump, unitig_dump, color_dump);
            log::info!("Serializing sparse-dense index to {}", out_path.display());
            write_index_variant(&IndexVariant::SparseDenseIndex(index), &mut out);
        },
        Subcommands::Export { index: index_path, color_dump_prefix, n_threads } => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, true); // Select support is required for export
            match index {
                IndexVariant::SparseDenseIndex(idx) => export_index(&idx, &color_dump_prefix, n_threads),
            };
        },
        Subcommands::Stats { index: index_path, n_threads } => {
            let index = load_index_variant(&index_path, true); // Select support required for DBG
            match index {
                IndexVariant::SparseDenseIndex(idx) => print_stats(&idx, n_threads),
            };
        }
        Subcommands::Report { index: index_path, input, output } => {
            log::info!("Loading color names from index");
            let color_names = load_index_color_names_only(&index_path);
            log::info!("Processing pseudoalignment file");
            match (input, output) {
                (Some(in_path), Some(out_path)) => report(BufReader::new(File::open(in_path).unwrap()), BufWriter::new(File::create(out_path).unwrap()), &color_names),
                (Some(in_path), None) => report(BufReader::new(File::open(in_path).unwrap()), BufWriter::new(std::io::stdout()), &color_names),
                (None, Some(out_path)) => report(BufReader::new(std::io::stdin()), BufWriter::new(File::create(out_path).unwrap()), &color_names),
                (None, None) => report(BufReader::new(std::io::stdin()), BufWriter::new(std::io::stdout()), &color_names),
            }
        },
        Subcommands::FinimizerStats { index: index_path, n_threads} => {
            // Hidden subcommand, might be deleted at any moment
            log::info!("Loading index");
            let index = load_index_variant(&index_path, true); // Select support is required for verify
            match index {
                IndexVariant::SparseDenseIndex(idx) => finimizers::minimizer_stats(&idx, n_threads, finimizers::MinimizerType::Finimizer),
            };
        },
        Subcommands::MinimizerStats { index: index_path, n_threads, m} => {
            // Hidden subcommand, might be deleted at any moment
            log::info!("Loading index");
            let index = load_index_variant(&index_path, true); // Select support is required for verify
            match index {
                IndexVariant::SparseDenseIndex(idx) => finimizers::minimizer_stats(&idx, n_threads, finimizers::MinimizerType::Minimizer(m)),
            };
        },
    }

    ExitCode::SUCCESS
}