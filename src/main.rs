#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material
#![allow(clippy::len_zero)] // !is_empty reads as "not is empty" which is not English 
#![allow(clippy::manual_is_multiple_of)] // Oh please

use std::{collections::HashSet, fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, path::{Path, PathBuf}, sync::Arc};
use bitmap_storage::BitmapStorage;
use clap::{Parser, Subcommand};
use colex_colored_kmers::CompactColexKmers;
use coloring_interface::{ColorSetOwned, ColorSetStorage, ColorSetView};
use io::ChainedInputStreamWithRevComp;
use parallel_ms_iteration::{DeduplicatingColorElementGenerator, DistinctColexComputation};
use sbwt::{BitPackedKmerSortingDisk, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix, reverse_complement_in_place};
use simple_sds_sbwt::ops::{BitVec, Rank};
use sparse_dense_storage::SparseDenseStorage;

use crate::{colex_colored_kmers::{ColexToColorSetMap, mark_key_kmers}, int_vec::CompactIntVec, io::ChainedInputStream, iterators::VecIterator, parallel_ms_iteration::MsElementGenerator, set_of_sets_construction::{ParallelElementGenerator, SetElement}};

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
mod old_merge;
mod new_merge;

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

        #[arg(help = "Precomputed SBWT file of k-mers (optional)", short, long)]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Output filename", short, long, required = true)]
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

    #[command(arg_required_else_help = true, name = "dump-color-names")]
    DumpColorNames{
        #[arg(long = "index", short = 'i', required = true)]
        index: PathBuf,
    },

    #[command(arg_required_else_help = true, name = "merge-compressed-indexes")]
    MergeCompressedIndexes {
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
        #[arg(help = "Precomputed bit matrix SBWT", long = "sbwt", short = 's')]
        sbwt: Option<PathBuf>,

        #[arg(help = "Index text dump file prefix, as written by Fulgor 4.0.0", long = "color-dump-prefix", short = 'c', required = true)]
        color_dump_prefix: PathBuf,

        #[arg(long = "index-type")]
        index_type: ColoringType,

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

struct MyBitmapStream {
    bs: BitmapStorage,
    pos: usize,
    buf: Vec<usize>,
}

impl crate::iterators::USizeIteratorGenerator for MyBitmapStream {
    type Iter<'a> = VecIterator<'a> where Self: 'a;
    
    fn next<'b>(&'b mut self) -> Option<Self::Iter<'b>> {
        if self.pos == self.bs.n_sets() {
            None
        } else {
            self.buf.clear(); 
            self.buf.extend(self.bs.get_set_view(self.pos).iter());
            self.pos += 1;
            Some(VecIterator::new(&self.buf))
        }
    }
}

fn build_coloring<CSS: ColorSetStorage + Send>(
    sbwt: Arc<sbwt::SbwtIndex<SubsetMatrix>>, lcs: LcsArray, input_paths: &[PathBuf], n_threads: usize, sample_distance: usize, from_unitigs: bool) -> CompactColexKmers<CSS> {

    let n_colors = input_paths.len();
    log::info!("Building distinct color set structure");

    log::info!("=== PHASE 1/3: Marking key k-mers ===");
    let phase1_input_stream = ChainedInputStreamWithRevComp::new(input_paths.to_owned());
    let key_kmer_marks = mark_key_kmers(&sbwt, &lcs, sample_distance, phase1_input_stream, n_threads);
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / sbwt.n_kmers() as f64 * 100.0);
    assert_eq!(key_kmer_marks.len(), sbwt.n_sets());

    log::info!("=== PHASE 2/3: Building color set finperprints for key k-mers ===");
    let random_seed = 123123; // Todo: be more random
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id) = if from_unitigs {
        let gen = MsElementGenerator::new(input_paths.to_owned(), StreamingIndex::new(&sbwt, &lcs));
        set_of_sets_construction::find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen, key_kmer_marks.clone(), sbwt.n_sets(), n_colors, n_threads, random_seed)
    } else {
        let gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, input_paths.to_owned());
        set_of_sets_construction::find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen, key_kmer_marks.clone(), sbwt.n_sets(), n_colors, n_threads, random_seed)
    };

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    let css = if from_unitigs {
        let gen = MsElementGenerator::new(input_paths.to_owned(), StreamingIndex::new(&sbwt, &lcs));
        set_of_sets_construction::build_color_set_storage(n_colors, repr_kmer_marks, distinct_set_sizes, gen, n_threads)
    } else {
        let gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, input_paths.to_owned());
        set_of_sets_construction::build_color_set_storage(n_colors, repr_kmer_marks, distinct_set_sizes, gen, n_threads)
    };

    log::info!("Building rank support for key k-mer marks");
    let mut key_kmer_marks = util::bitvec_to_simple_sds_bitvec(key_kmer_marks);
    key_kmer_marks.enable_rank();
    assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));
    let colex_map = ColexToColorSetMap {
        sbwt: sbwt.clone(), // Clones just the Arc
        sampling: key_kmer_marks, 
        color_set_ids: key_kmer_idx_to_set_id,
    };

    let color_names: Vec<String> = input_paths.iter().map(|p| p.to_string_lossy().to_string()).collect();
    CompactColexKmers::<CSS>::new(sbwt, lcs, colex_map, css, Some(&color_names))
}

#[allow(clippy::large_enum_variant)] // It's saying that it's almost a kilobyte. I don't understand why but ok.
enum IndexVariant {
    BitmapIndex(CompactColexKmers<BitmapStorage>),
    SparseDenseIndex(CompactColexKmers<SparseDenseStorage>),
}

fn load_index_variant(path: &Path, build_select: bool) -> IndexVariant {
    let mut input = BufReader::new(File::open(path).unwrap());
    let mut id_buf = [0u8; 8];
    input.read_exact(&mut id_buf).unwrap();
    if id_buf == ColoringType::Bitmaps.serialization_id() {
        let index = CompactColexKmers::<BitmapStorage>::load(&mut input, build_select);
        IndexVariant::BitmapIndex(index)
    } else if id_buf == ColoringType::SparseDense.serialization_id() {
        let index = CompactColexKmers::<SparseDenseStorage>::load(&mut input, build_select);
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

fn print_color_names<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>) {
    for (id, name) in index.get_color_names().iter().enumerate() {
        println!("{}\t{}", id, name);
    }
}


fn print_color_sets<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, query_path: &Path, print_kmers: bool) {
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
fn intersection_pseudoalignment<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, query_path: &Path, min_hits: usize) {
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
fn threshold_pseudoalignment<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, query_path: &Path, min_hits: usize, threshold: f64, denominator: Denominator) {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_path).unwrap();
    // Buffered writing to stdout
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout);
    let mut query_idx = 0_usize;
    let n_colors = index.get_set_storage().get_full_set().len();
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
                    (IndexVariant::BitmapIndex(c1), IndexVariant::BitmapIndex(c2)) => {
                        log::info!("Merging bitmap indexes");
                        let merged_colored_kmers = new_merge::new_merge(c1, c2, low_ram_mode, sample_distance, n_threads);
                        log::info!("Serializing merged index to {}", outpath.display());
                        write_index_variant(&IndexVariant::BitmapIndex(merged_colored_kmers), &mut out);
                    },
                    (IndexVariant::SparseDenseIndex(c1), IndexVariant::SparseDenseIndex(c2)) => {
                        log::info!("Merging sparse-dense indexes");
                        let merged_colored_kmers = new_merge::new_merge(c1, c2, low_ram_mode, sample_distance, n_threads);
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

fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { input: input_fof, output, temp_dir, k, n_threads, index_type, sample_distance, sbwt_path, from_unitigs} => {
            if k % 2 == 0 && from_unitigs {
                panic!("--from_unitigs requires odd k");
            }

            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            let input_stream = io::ChainedInputStream::new(input_paths.clone());
            let mut out = BufWriter::new(File::create(&output).unwrap());

            let (sbwt, lcs) = if let Some(sbwt_path) = sbwt_path {
                log::info!("Loading SBWT from {}", sbwt_path.display());
                let mut sbwt_in = BufReader::new(File::open(sbwt_path).unwrap());
                let sbwt::SbwtIndexVariant::SubsetMatrix(mut sbwt) = sbwt::load_sbwt_index_variant(&mut sbwt_in).unwrap();

                log::info!("Building select support for SBWT");
                sbwt.build_select();
                log::info!("Building LCS array");
                let lcs = LcsArray::from_sbwt(&sbwt, n_threads); // TODO: load LCS
                (sbwt, lcs)
            } else {
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
                let sbwt = sbwt;
                let lcs = lcs.unwrap(); // Ok because we used .build_lcs(true)
                (sbwt, lcs)
            };

            match index_type {
                ColoringType::Bitmaps => {
                    let index = build_coloring::<BitmapStorage>(Arc::new(sbwt), lcs, &input_paths, n_threads, sample_distance, from_unitigs);
                    log::info!("Serializing bitmap index to {}", output.display());
                    write_index_variant(&IndexVariant::BitmapIndex(index), &mut out);
                },
                ColoringType::SparseDense => {
                    let index = build_coloring::<SparseDenseStorage>(Arc::new(sbwt), lcs, &input_paths, n_threads, sample_distance, from_unitigs);
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
        Subcommands::DumpColorNames{ index: index_path} => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, false); // No select support required
            match index {
                IndexVariant::BitmapIndex(idx) => print_color_names(&idx),
                IndexVariant::SparseDenseIndex(idx) => print_color_names(&idx),
            };
        },
        Subcommands::MergeCompressedIndexes { index_file_list, temp_dir, outfile, n_threads, low_ram_mode, sample_distance } => {
            let infiles: Vec<PathBuf> = BufReader::new(File::open(index_file_list).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            run_merge_tree(&infiles, &temp_dir, &outfile, n_threads, low_ram_mode, sample_distance);
        },
        Subcommands::Import { sbwt: sbwt_path, color_dump_prefix, out: out_path, n_threads, temp_dir, index_type, sample_distance} => {
            let unitig_filename = format!("{}.unitigs.fa", color_dump_prefix.to_str().unwrap());
            let color_sets_filename = format!("{}.color_sets.txt", color_dump_prefix.to_str().unwrap());
            let metadata_filename = format!("{}.metadata.txt", color_dump_prefix.to_str().unwrap());

            // Try to open to check that the files are found
            BufReader::new(File::open(&unitig_filename).unwrap());
            BufReader::new(File::open(&color_sets_filename).unwrap());
            BufReader::new(File::open(&metadata_filename).unwrap());

            log::info!("Reading metadata from {}", metadata_filename);
            let metadata = index_import::read_index_dump_metadata(BufReader::new(File::open(&metadata_filename).unwrap()));

            let (sbwt, lcs) = if let Some(sbwt_path) = sbwt_path {
                log::info!("Loading SBWT from {}", sbwt_path.display());
                let mut sbwt_in = BufReader::new(File::open(sbwt_path).unwrap());
                let sbwt::SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt::load_sbwt_index_variant(&mut sbwt_in).unwrap();
                if sbwt.k() != metadata.k {
                    log::error!("SBWT k does not match the index dump k ({} vs {})", sbwt.k(), metadata.k);
                    return;
                }

                log::info!("Building LCS array");
                let lcs = LcsArray::from_sbwt(&sbwt, n_threads);
                (sbwt, lcs)
            } else {
                log::info!("No precomputed SBWT given. Building the SBWT and the LCS array");
                let input_stream = io::ChainedInputStream::new(vec![PathBuf::from(&unitig_filename)]);
                let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
                    .add_rev_comp(true)
                    .k(metadata.k)
                    .build_lcs(true)
                    .n_threads(n_threads)
                    .precalc_length(8)
                    .algorithm(BitPackedKmerSortingDisk::new().dedup_batches(false).temp_dir(&temp_dir))
                .run(input_stream); // No batch dedup because unitigs should not have duplicates
                (sbwt, lcs.unwrap())
            };

            let mut out = BufWriter::new(File::create(&out_path).unwrap());

            let unitig_dump = BufReader::new(File::open(&unitig_filename).unwrap());
            let color_dump = BufReader::new(File::open(&color_sets_filename).unwrap());
            let metadata_dump = BufReader::new(File::open(&metadata_filename).unwrap());

            match index_type {
                ColoringType::Bitmaps => {
                    let index = CompactColexKmers::<BitmapStorage>::new_from_colored_unitig_dump(
                        sbwt, lcs, sample_distance, n_threads, metadata_dump, unitig_dump, color_dump);
                    log::info!("Serializing bitmap index to {}", out_path.display());
                    write_index_variant(&IndexVariant::BitmapIndex(index), &mut out);
                },
                ColoringType::SparseDense => {
                    let index = CompactColexKmers::<SparseDenseStorage>::new_from_colored_unitig_dump(
                        sbwt, lcs, sample_distance, n_threads, metadata_dump, unitig_dump, color_dump);
                    log::info!("Serializing sparse-dense index to {}", out_path.display());
                    write_index_variant(&IndexVariant::SparseDenseIndex(index), &mut out);
                },
            }
        },
        Subcommands::Export { index: index_path, color_dump_prefix, n_threads } => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, true); // Select support is required for export
            match index {
                IndexVariant::BitmapIndex(idx) => export_index(&idx, &color_dump_prefix, n_threads),
                IndexVariant::SparseDenseIndex(idx) => export_index(&idx, &color_dump_prefix, n_threads),
            };
        },
        Subcommands::FinimizerStats { index: index_path, n_threads} => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, true); // Select support is required for verify
            match index {
                IndexVariant::BitmapIndex(idx) => finimizers::minimizer_stats(&idx, n_threads, finimizers::MinimizerType::Finimizer),
                IndexVariant::SparseDenseIndex(idx) => finimizers::minimizer_stats(&idx, n_threads, finimizers::MinimizerType::Finimizer),
            };
        },
        Subcommands::MinimizerStats { index: index_path, n_threads, m} => {
            log::info!("Loading index");
            let index = load_index_variant(&index_path, true); // Select support is required for verify
            match index {
                IndexVariant::BitmapIndex(idx) => finimizers::minimizer_stats(&idx, n_threads, finimizers::MinimizerType::Minimizer(m)),
                IndexVariant::SparseDenseIndex(idx) => finimizers::minimizer_stats(&idx, n_threads, finimizers::MinimizerType::Minimizer(m)),
            };
        },
            }
/*
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