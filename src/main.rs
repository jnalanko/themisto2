#![allow(non_snake_case)] // Using upper-case variable names from the source material

use std::{cmp::min, fs::File, io::{BufRead, BufReader, BufWriter, Stdout, Write}, ops::Range, path::{Path, PathBuf}, process::ExitCode, str::FromStr, time::Instant};
use clap::{Parser, Subcommand};
use colex_colored_kmers::{Coloring, ColoredKmers, DbgIndex};
use coloring_interface::{ColorSetStorage, ColorSetView};
use io::RewindableSeqStreamGenerator;
use parallel_ms_iteration::{DeduplicatingColorElementGenerator};
use sbwt::{BitPackedKmerSortingDisk, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix, dbg::Dbg, sbwt_index_variant::SbwtIndexVariant};
use simple_sds_sbwt::ops::{BitVec, Rank};
use sparse_dense_storage::SparseDenseStorage;

use themisto2::{index_import, colex_colored_kmers, coloring_interface, sparse_dense_storage, io, set_of_sets_construction, parallel_ms_iteration, finimizers, pseudoalignment, pseudoalignment_metrics, filter_reads, unitig_export};
use themisto2::{colex_colored_kmers::{ColexToColorSetMap, mark_key_kmers}, parallel_ms_iteration::MsElementGenerator, report::report};
use themisto2::set_operations::{difference, intersect, merge, subset};
use themisto2::{Component, IndexDir, IndexWriter, PseudoalignmentIndex};


#[derive(Parser)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommands,
}

#[derive(clap::ValueEnum, Clone, Copy, Debug)]
pub enum Denominator { // Options for the CLI
    All,
    Relevant,
    MaxHits,
}

#[derive(Subcommand)]
pub enum Subcommands {
    #[command(arg_required_else_help = true)]
    Build {

        // The maximum number of colors is 2^32 because we use 32-bit color ids in phase 2 of the construction.
        #[arg(help = "A file with one fasta/fastq filename per line, optionally followed by a tab and a color name. A line without a color name uses the filename as the color name. Maximum supported number of colors: 2^32.", long = "file-colors")]
        color_fof: Option<PathBuf>,

        // The maximum number of colors is 2^32 because we use 32-bit color ids in phase 2 of the construction.
        #[arg(help = "A fasta/fastq file, with one color per sequence. The color name of a sequence is its header up to the first space. Maximum supported number of colors: 2^32.", long = "seq-colors")]
        sequence_colors_file: Option<PathBuf>,

        #[arg(help = "With --seq-colors: sequences with the same color name form a single color, instead of each sequence having its own color. IMPORTANT: all sequences with the same color name must be CONSECUTIVE in the file; otherwise the build stops with an error before it starts. All sequences of a color are held in memory at once while the color is processed (several colors at a time, up to a few per thread), so for colors too large for that, use --file-colors with one or more files per color instead, which reads them one sequence at a time.", long = "seq-colors-by-name", requires = "sequence_colors_file")]
        seq_colors_by_name: bool,

        #[arg(help = "Precomputed bit matrix SBWT (optional). Maximum supported SBWT length: 2^40.", long = "sbwt", short = 's')]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Precomputed LCS array for the SBWT (optional)", long = "lcs", short = 'l')]
        lcs_path: Option<PathBuf>,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty. If --to-disk-in-pieces is given, this is interpreted as a filename prefix instead.", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Optional: Build from per-color unitigs (requires odd k). This makes the construction much faster, because the k-mers of a color are not deduplicated. IMPORTANT: the user is responsible for making sure that no k-mer occurs more than once within a color, counting both orientations, over all sequences of that color (for example, the unitigs of each color from GGCAT). This is not checked: if it does not hold, the index is silently wrong. A color may consist of several files (--file-colors) or several sequences (--seq-colors-by-name).", long)]
        from_unitigs: bool,

        #[arg(help = "Directory for temporary files. Required if --sbwt is not given. Should point to a location with fast IO and plenty of free space.", long = "temp-dir", required_unless_present="sbwt_path")]
        temp_dir: Option<PathBuf>,

        #[arg(short, required = true)]
        k: usize,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(help = "Number of parallel worker threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Number of parallel input parsing threads (only matters with --file-colors). This might speed up construction if the file system is slow or the number of input sequence files is massive. If you observe that the program is using less CPU time than what is available during construction phases 2 or 3, increasing the number of parsing threads might help with keeping the cores busy. Only use this if the storage hardware actually supports concurrent reads (i.e. not a single spinning HDD). Otherwise this will probably just make things slower.", long = "n-parser-threads", default_value = "1")]
        n_parser_threads: usize,
        
        // Hidden because the output is three files, not a single themisto2 index. TODO: Should put the files together into a single index.
        #[arg(long = "to-disk-in-pieces", help = "Build the final sets in pieces, flushing to disk after each piece. This is useful if there is not enough memory to hold the final data structure in memory at once.", hide = true)]
        n_pieces: Option<usize>,

    },

    #[command(arg_required_else_help = true, name = "intersection-pseudoalign")]
    IntersectionPseudoalign {
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', conflicts_with_all = ["query_list", "query_output_list"], help = "Query file (FASTA/FASTQ, optionally gzipped). If omitted and --query-list is also omitted, reads queries from stdin. Mutually exclusive with --query-list and --query-output-list.")]
        query: Option<PathBuf>,

        #[arg(long = "output", short = 'o', conflicts_with_all = ["query_list", "query_output_list"], help = "Output file. If omitted, writes to stdout. Mutually exclusive with --query-list and --query-output-list.")]
        output: Option<PathBuf>,

        #[arg(long = "query-list", help = "A list of query files, one per line. All query results are printed to stdout, unless --query-output-list is given.")]
        query_list: Option<PathBuf>,

        #[arg(long = "query-output-list", requires = "query_list", help = "A list of filenames, one per line. The file on the i-th line is the output file for the i-th input file in the file for --query-list. Requires --query-list.")]
        query_output_list: Option<PathBuf>,

        #[arg(long = "min-hits", short = 'm', default_value = "1")]
        min_hits: usize,

        #[arg(long = "n-threads", short = 't', required = true, default_value = "4")]
        n_threads: usize,

        #[arg(long = "report-hit-counts", default_value = "false")]
        report_hit_counts: bool,

        #[arg(long = "positive-reads-only", default_value = "false", help = "Output formatting: only write reads with at least one positive color.")]
        positive_reads_only: bool,

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
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "query", short = 'q', conflicts_with_all = ["query_list", "query_output_list"], help = "Query file (FASTA/FASTQ, optionally gzipped). If omitted and --query-list is also omitted, reads queries from stdin. Mutually exclusive with --query-list and --query-output-list.")]
        query: Option<PathBuf>,

        #[arg(long = "output", short = 'o', conflicts_with_all = ["query_list", "query_output_list"], help = "Output file. If omitted, writes to stdout. Mutually exclusive with --query-list and --query-output-list.")]
        output: Option<PathBuf>,

        #[arg(long = "query-list", help = "A list of query files, one per line. All query results are printed to stdout, unless --query-output-list is given.")]
        query_list: Option<PathBuf>,

        #[arg(long = "query-output-list", requires = "query_list", help = "A list of filenames, one per line. The file on the i-th line is the output file for the i-th input file in the file for --query-list. Requires --query-list.")]
        query_output_list: Option<PathBuf>,

        #[arg(long = "min-hits", short = 'm', default_value = "1")]
        min_hits: usize,

        #[arg(long = "threshold", short = 'd', required = true)]
        threshold: f64,

        #[arg(long = "denominator", short = 'n')]
        denominator: Denominator,

        #[arg(long = "n-threads", short = 't', required = true, default_value = "4")]
        n_threads: usize,

        #[arg(long = "report-hit-counts", default_value = "false")]
        report_hit_counts: bool,

        #[arg(long = "positive-reads-only", default_value = "false", help = "Output formatting: only write reads with at least one positive color.")]
        positive_reads_only: bool,

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
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
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
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "output", short = 'o', help = "Output file. If omitted, writes to stdout.")]
        output: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, name = "merge")]
    Merge {
        #[arg(help = "A file with one Themisto 2 index directory per line", long = "index-file-list", required = true)]
        index_file_list: PathBuf,

        #[arg(long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty.", long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "low-ram-mode", help = "Use slower algorithms that need less memory to combine the SBWTs. The result is the same.")]
        low_ram_mode: bool,

        #[arg(long = "merge-shared-colors", help = "Treat colors with identical names in different input indexes as the same color. Color names must be unique within each input index.")]
        merge_shared_colors: bool,

        #[arg(long = "keep-redundant-dummies", help = "Do not remove the dummy nodes that become redundant when merging SBWTs. The result has the same k-mers and colors, but may have more dummy nodes. With more than two input indexes, the intermediate merges always keep them, and only the last merge removes them (unless this option is given).")]
        keep_redundant_dummies: bool,
    },

    #[command(arg_required_else_help = true, name = "intersect", about = "Intersect two indexes: keep the k-mers that are in both")]
    Intersect {
        #[arg(help = "Themisto 2 index directory", long = "index1", required = true)]
        index1: PathBuf,

        #[arg(help = "Themisto 2 index directory", long = "index2", required = true)]
        index2: PathBuf,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty.", long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "colors", value_enum, default_value = "union", help = "How the color set of each k-mer is computed from its color sets in the two indexes. With 'intersect', colors are matched by name, only the colors in both indexes are kept, and k-mers whose color sets do not intersect get an empty color set. Color names must be unique within each input index.")]
        colors: intersect::IntersectColors,

        #[arg(long = "merge-shared-colors", help = "With '--colors union': treat colors with identical names in the two indexes as the same color. Color names must be unique within each input index.")]
        merge_shared_colors: bool,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "low-ram-mode", help = "Use slower algorithms that need less memory to combine the SBWTs. The result is the same.")]
        low_ram_mode: bool,
    },

    #[command(arg_required_else_help = true, name = "subset", about = "Keep the k-mers of an index whose color set contains at least one of the given colors, with only those colors: the same as an index built from the sequences of these colors")]
    Subset {
        #[arg(help = "Themisto 2 index directory", long = "index", required = true)]
        index: PathBuf,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty.", long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "colors", value_delimiter = ',', help = "Comma-separated names of the colors to keep (to remove, with --complement)")]
        colors: Vec<String>,

        #[arg(long = "colors-file", help = "A file with the names of the colors to keep (to remove, with --complement), one per line")]
        colors_file: Option<PathBuf>,

        #[arg(long = "complement", help = "Remove the given colors and keep all others")]
        complement: bool,

        #[arg(long = "keep-all-colors", help = "Keep the whole color sets of the kept k-mers, with all colors of the index, instead of only the given colors")]
        keep_all_colors: bool,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "low-ram-mode", help = "Use slower algorithms that need less memory to combine the SBWTs. The result is the same.")]
        low_ram_mode: bool,
    },

    #[command(arg_required_else_help = true, name = "difference", about = "Subtract two indexes: keep the k-mers of index1 that are not in index2, with their colors in index1")]
    Difference {
        #[arg(help = "Themisto 2 index directory", long = "index1", required = true)]
        index1: PathBuf,

        #[arg(help = "Themisto 2 index directory", long = "index2", required = true)]
        index2: PathBuf,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty.", long = "output", short = 'o', required = true)]
        outfile: PathBuf,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "low-ram-mode", help = "Use slower algorithms that need less memory to combine the SBWTs. The result is the same.")]
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

        #[arg(help = "Directory for temporary files. Required if --sbwt is not given. Should point to a location with fast IO and plenty of free space.", long = "temp-dir", required_unless_present="sbwt_path")]
        temp_dir: Option<PathBuf>,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty.", long = "out", short = 'o', required = true)]
        out: PathBuf,
    },

    #[command(arg_required_else_help = true, name = "import-sbwt")]
    ImportSbwt {
        #[arg(help = "Precomputed bit matrix SBWT", long = "sbwt", short = 's', required = true)]
        sbwt_path: PathBuf,

        #[arg(help = "Precomputed LCS array for the SBWT (optional)", long = "lcs", short = 'l')]
        lcs_path: Option<PathBuf>,

        #[arg(help = "Output directory for the Themisto 2 index. It must not exist or must be empty.", long = "output", short = 'o', required = true)]
        output: PathBuf,

        #[arg(long = "sample-distance", short = 'd', default_value = "30")]
        sample_distance: usize,

        #[arg(help = "Number of parallel threads", long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(help = "Name for the single color in the resulting Themisto 2 index", long = "color-name", required = true)]
        color_name: String,
    },

    #[command(arg_required_else_help = true)]
    Export {
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(help = "Output file prefix", long = "output-prefix", short = 'o', required = true)]
        color_dump_prefix: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,

        #[arg(long = "gfa", action = clap::ArgAction::SetTrue, help = "Export as GFA 1.0 with links (written to <prefix>.unitigs.gfa)")]
        gfa: bool,

        #[arg(long = "colors-to-stdout", action = clap::ArgAction::SetTrue, help = "Write the color dump to stdout instead of a file.")]
        colors_to_stdout: bool,

        #[arg(long = "no-unitigs", action = clap::ArgAction::SetTrue, help = "Do not export unitig sequences.")]
        no_unitigs: bool,

        #[arg(long = "no-color-sets", action = clap::ArgAction::SetTrue, help = "Do not export color sets.")]
        no_color_sets: bool,
    },
    #[command(arg_required_else_help = true, name = "stats")]
    Stats {
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
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

    #[command(arg_required_else_help = true, name = "filter-reads")]
    FilterReads {
        #[arg(help = "Pseudoalignment JSONL file produced by a pseudoalign subcommand", long = "pseudoalignment-file", short = 'p', required = true)]
        pseudoalignment: PathBuf,

        #[arg(help = "Reads file (FASTA/FASTQ, optionally gzipped) that was the input to the pseudoalignment", long = "reads", short = 'r', required = true)]
        reads: PathBuf,

        #[arg(help = "Output reads file. Extension determines format (.fa/.fq, optionally .gz)", long = "output", short = 'o', required = true)]
        output: PathBuf,

        #[arg(help = "File with one color id (integer) per line. A read is kept if its pseudoalignment matches at least one of these colors.", long = "color-ids", short = 'c', required = true)]
        color_ids: PathBuf,
    },

    #[command(arg_required_else_help = true, hide = true)]
    FinimizerStats {
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
        index: PathBuf,

        #[arg(long = "n-threads", short = 't', default_value = "4")]
        n_threads: usize,
    },
    #[command(arg_required_else_help = true, hide = true)]
    MinimizerStats {
        #[arg(help = "Themisto 2 index directory", long = "index", short = 'i', required = true)]
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
    FileColors(Vec<Vec<PathBuf>>), // One group of files per color
    SequenceColors(PathBuf, Vec<usize>), // One file, and the color of each of its sequences
}

impl ColoredSeqInput {
    fn get_generator(&self) -> Box<dyn RewindableSeqStreamGenerator + Sync + Send> {
        match self {
            ColoredSeqInput::FileColors(file_groups) => {
                Box::new(io::SeqStreamGeneratorFromFiles::new(file_groups.clone()))
            },
            ColoredSeqInput::SequenceColors(path_buf, record_colors) => {
                Box::new(io::SeqStreamGeneratorFromSingleFile::new(path_buf.clone(), record_colors.clone()))
            }
        }
    }
}

/// Parses a --file-colors file. Each non-empty line is a path, optionally followed by a tab and
/// a color name. A line without a name uses the path as its name. Lines with the same name are
/// one color. Colors are numbered in the order their names first appear. Returns the files of
/// each color and the color names.
fn parse_color_fof(reader: impl BufRead) -> Result<(Vec<Vec<PathBuf>>, Vec<String>), String> {
    let mut file_groups = Vec::<Vec<PathBuf>>::new();
    let mut color_names = Vec::<String>::new();
    let mut name_to_color = std::collections::HashMap::<String, usize>::new();
    for (line_idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading the --file-colors file: {}", e))?;
        if line.is_empty() { continue; }

        let (path, name) = match line.split_once('\t') {
            Some((path, name)) => (path, name),
            None => (line.as_str(), line.as_str()),
        };
        if path.is_empty() || name.is_empty() {
            return Err(format!("Line {} of the --file-colors file has an empty path or color name: {:?}", line_idx + 1, line));
        }

        let color = *name_to_color.entry(name.to_string()).or_insert_with(|| {
            file_groups.push(Vec::new());
            color_names.push(name.to_string());
            color_names.len() - 1
        });
        file_groups[color].push(PathBuf::from(path));
    }
    Ok((file_groups, color_names))
}

/// Assigns colors to the sequences of a --seq-colors file, given the color name of each
/// sequence. Without by_name, each sequence is its own color. With by_name, sequences with the
/// same name are one color, and colors are numbered in file order. The sequences of a color
/// must then be consecutive in the file, which is checked here so that the build fails before
/// it starts. Returns the color of each sequence and the color names.
fn assign_seq_colors(seq_names: Vec<String>, by_name: bool) -> Result<(Vec<usize>, Vec<String>), String> {
    if !by_name {
        return Ok(((0..seq_names.len()).collect(), seq_names));
    }

    let mut record_colors = Vec::<usize>::with_capacity(seq_names.len());
    let mut color_names = Vec::<String>::new();
    let mut seen_names = std::collections::HashSet::<String>::new();
    for (seq_idx, name) in seq_names.into_iter().enumerate() {
        if name.is_empty() {
            return Err(format!("Sequence {} of the --seq-colors file has an empty color name, which --seq-colors-by-name can not group", seq_idx + 1));
        }
        if color_names.last() != Some(&name) {
            if !seen_names.insert(name.clone()) {
                return Err(format!("Sequence {} of the --seq-colors file has color name {:?}, but the sequences with that name are not consecutive in the file. With --seq-colors-by-name, all sequences with the same color name must be consecutive.", seq_idx + 1, name));
            }
            color_names.push(name);
        }
        record_colors.push(color_names.len() - 1);
    }
    Ok((record_colors, color_names))
}

// Returns the index and the LCS array of its SBWT if BuildMode is InMemory, otherwise serializes
// the color sets as a set of files to disk.
#[allow(clippy::too_many_arguments)]
fn build_coloring<CSS: ColorSetStorage + Send>(sbwt: sbwt::SbwtIndex<SubsetMatrix>, lcs: LcsArray, input_mode: ColoredSeqInput, color_names: &[String], n_threads: usize, sample_distance: usize, from_unitigs: bool, build_mode: BuildMode, n_parser_threads: usize) -> Option<(DbgIndex<CSS>, LcsArray)>{

    log::info!("Building distinct color set structure");
    let n_colors = u32::try_from(color_names.len()).unwrap_or_else( |_| {
        log::error!("Maximum number of colors 2^32 exceeded");
        panic!();
    });

    if sbwt.n_sets() > set_of_sets_construction::SET_OF_SETS_CONSTRUCTION_MAX_SBWT_LEN {
        log::error!("Maximum SBWT length exceeded. Found: {}, maximum supported: {}", sbwt.n_sets(), set_of_sets_construction::SET_OF_SETS_CONSTRUCTION_MAX_SBWT_LEN);
        panic!();
    }

    if let ColoredSeqInput::FileColors(v) = &input_mode {
        assert!(v.len() == n_colors as usize, "Number of color names does not match the number of input file groups");
    }

    log::info!("Initializing DBG");
    let dbg = Dbg::new(&sbwt, Some(&lcs), n_threads);

    log::info!("=== PHASE 1/3: Marking key k-mers ===");
    let key_kmer_marks = {
        let mut phase1_input_stream = input_mode.get_generator();
        mark_key_kmers(&sbwt, &dbg, sample_distance, &mut phase1_input_stream, n_threads, n_parser_threads, true)
    };
    drop(dbg); // Only needed in phase 1
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / sbwt.n_kmers() as f64 * 100.0);
    assert_eq!(key_kmer_marks.len(), sbwt.n_sets());

    log::info!("=== PHASE 2/3: Building color set finperprints for key k-mers ===");
    let color_stream_gen = input_mode.get_generator();
    let random_seed = 123123; // Todo: be more random
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id, mut key_kmer_marks) = if from_unitigs {
        let ms_gen = MsElementGenerator::new(color_stream_gen, StreamingIndex::new(&sbwt, &lcs), true, n_parser_threads);
        set_of_sets_construction::find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(ms_gen, key_kmer_marks, n_colors, n_threads, random_seed)
    } else {
        let ms_gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, color_stream_gen, true);
        set_of_sets_construction::find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(ms_gen, key_kmer_marks, n_colors, n_threads, random_seed)
    };

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    match build_mode {
        BuildMode::InMemory => {
            let color_stream_gen = input_mode.get_generator();
            let css = if from_unitigs {
                let gen = MsElementGenerator::new(color_stream_gen, StreamingIndex::new(&sbwt, &lcs), true, n_parser_threads);
                set_of_sets_construction::build_color_set_storage(n_colors as usize, repr_kmer_marks, distinct_set_sizes, gen, n_threads)
            } else {
                let gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, color_stream_gen, true);
                set_of_sets_construction::build_color_set_storage(n_colors as usize, repr_kmer_marks, distinct_set_sizes, gen, n_threads)
            };

            log::info!("Building rank support for key k-mer marks");
            key_kmer_marks.enable_rank(); // Should already have but does not hurt
            assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));
            let colex_map = ColexToColorSetMap {
                sampling: key_kmer_marks, 
                color_set_ids: key_kmer_idx_to_set_id,
            };

            let coloring = Coloring::new(colex_map, css, Some(color_names));
            log::info!("Marking the suffix groups and dummy k-mers");
            let index = DbgIndex::from_lcs(sbwt, &lcs, coloring, n_threads);
            Some((index, lcs))
        },
        BuildMode::ToDisk(out_prefix, n_pieces) => {
            assert!(n_pieces != 0);
            let file_groups = match &input_mode {
                ColoredSeqInput::FileColors(file_groups) => file_groups,
                ColoredSeqInput::SequenceColors(..) => {
                    panic!("ToDisk mode with sequence colors is not yet supported");
                    // The issue is that the current code chunks the input files.
                    // It takes some work to translate this behaviour to a single
                    // input file.
                },
            };
            let chunk_size = (n_colors as usize).div_ceil(n_pieces);
            if from_unitigs {
                let mut gens = Vec::<(MsElementGenerator, Range::<usize>)>::new();
                for (chunk_id, chunk) in file_groups.chunks(chunk_size).enumerate() {
                    let color_id_range = chunk_id*chunk_size .. min((chunk_id+1)*chunk_size, n_colors as usize);
                    let gen = Box::new(io::SeqStreamGeneratorFromFiles::new(chunk.to_owned()));
                    let ms_gen = MsElementGenerator::new(gen, StreamingIndex::new(&sbwt, &lcs), true, n_parser_threads);
                    gens.push((ms_gen, color_id_range));
                }
                set_of_sets_construction::build_color_set_storage_to_disk::<CSS>(repr_kmer_marks, distinct_set_sizes, gens, &out_prefix, n_threads);
            } else {
                let mut gens = Vec::<(DeduplicatingColorElementGenerator, Range::<usize>)>::new();
                for (chunk_id, chunk) in file_groups.chunks(chunk_size).enumerate() {
                    let color_id_range = chunk_id*chunk_size .. min((chunk_id+1)*chunk_size, n_colors as usize);
                    let gen = Box::new(io::SeqStreamGeneratorFromFiles::new(chunk.to_owned()));
                    let ms_gen = DeduplicatingColorElementGenerator::new(&sbwt, &lcs, gen, true);
                    gens.push((ms_gen, color_id_range));
                }
                set_of_sets_construction::build_color_set_storage_to_disk::<CSS>(repr_kmer_marks, distinct_set_sizes, gens, &out_prefix, n_threads);
            };
            None
        },
    }
}

/// Returns the value of `result`, or logs the error and exits.
fn or_exit<T>(result: std::io::Result<T>) -> T {
    result.unwrap_or_else(|e| {
        log::error!("{}", e);
        std::process::exit(1);
    })
}

/// Opens the Themisto 2 index at `path` and checks that it has all of `components`, so that a
/// missing component is reported before anything is loaded.
fn open_index(path: &Path, components: &[Component]) -> IndexDir {
    let index = or_exit(IndexDir::open(path));
    or_exit(index.require(components));
    index
}

const PSEUDOALIGNMENT_COMPONENTS: [Component; 5] = [Component::Sbwt, Component::Sufgroups, Component::ColorNames, Component::ColorSampling, Component::ColorStorage];
const DBG_COMPONENTS: [Component; 6] = [Component::Sbwt, Component::Sufgroups, Component::DummyMarks, Component::ColorNames, Component::ColorSampling, Component::ColorStorage];

/// The components of an index for pseudoalignment: the SBWT (no select support needed), the
/// suffix group marks and the coloring.
fn load_for_pseudoalignment(path: &Path) -> (SbwtIndex<SubsetMatrix>, bitvec::vec::BitVec, Coloring<SparseDenseStorage>) {
    let index = open_index(path, &PSEUDOALIGNMENT_COMPONENTS);
    (or_exit(index.load_sbwt(false)), or_exit(index.load_sufgroups()), or_exit(index.load_coloring()))
}

fn load_dbg_index(path: &Path) -> DbgIndex<SparseDenseStorage> {
    or_exit(open_index(path, &DBG_COMPONENTS).load_dbg_index())
}

/// Creates the output directory of an index. Called before the index is built, so that an
/// unusable output path is reported before the work is done.
fn create_index_writer(path: &Path) -> IndexWriter {
    or_exit(IndexWriter::create(path))
}

fn write_index(writer: IndexWriter, index: &DbgIndex<SparseDenseStorage>, lcs: &LcsArray, both_orientations: Option<bool>) {
    log::info!("Writing the Themisto 2 index to {}", writer.path().display());
    or_exit(writer.write_all(index, lcs, both_orientations));
}

/// The orientation of the result of a set operation on indexes with orientations `a` and `b`:
/// known only if it is known and the same for both.
fn combined_orientation(a: Option<bool>, b: Option<bool>) -> Option<bool> {
    match (a, b) {
        (Some(a), Some(b)) if a == b => Some(a),
        _ => None,
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

fn print_stats(index: &DbgIndex<SparseDenseStorage>, index_dir: &IndexDir, n_threads: usize) {
    // Compute quick-to-access statistics first

    let stdout = std::io::stdout();

    // Macro to flush stdout after every print
    macro_rules! printstat {
        ($($arg:tt)*) => {
            println!($($arg)*);
            stdout.lock().flush().unwrap();
        };
    }

    printstat!("n_kmers\t{}", index.sbwt().n_kmers());
    printstat!("n_colors\t{}", index.get_set_storage().n_colors());
    printstat!("sbwt_length\t{}", index.sbwt().n_sets());
    printstat!("n_sampled_positions\t{}", index.get_map().sampling.count_ones());

    let n_color_sets = index.get_set_storage().n_sets();
    printstat!("n_distinct_color_sets\t{}", n_color_sets);
    log::info!("Counting sparse and dense color sets");
    let n_dense_color_sets = index.get_set_storage().get_dense_marks().count_ones();
    let n_sparse_color_sets = n_color_sets - n_dense_color_sets;
    printstat!("n_sparse_color_sets\t{}", n_sparse_color_sets);
    printstat!("n_dense_color_sets\t{}", n_dense_color_sets);

    log::info!("Computing mean distinct color set size");
    let total_distinct_size: usize = (0..n_color_sets).map(|i| index.get_set_storage().get_set_view(i).len()).sum();
    printstat!("mean_distinct_color_set_size\t{:.6}", total_distinct_size as f64 / n_color_sets as f64);

    // The sizes of the files of the components. NA for components that the index does not have.
    let mut full_size_bytes = 0;
    for (component, size) in index_dir.component_sizes() {
        let name = component.file_name().replace('.', "_");
        match size {
            Some(size) => {
                printstat!("{}_size_bytes\t{}", name, size);
                full_size_bytes += size;
            },
            None => { printstat!("{}_size_bytes\tNA", name); },
        }
    }
    printstat!("full_index_size_bytes\t{}", full_size_bytes);

    // DBG traversal. Slow: computes unitig and k-mer color set stats
    log::info!("Computing unitig and k-mer color set stats (DBG traversal)");
    let stats = index.compute_index_stats(n_threads);
    printstat!("n_unitigs\t{}", stats.n_unitigs);
    printstat!("mean_unitig_length\t{:.6}", stats.mean_unitig_length());
    printstat!("mean_kmer_color_set_size\t{:.6}", stats.mean_kmer_color_set_size());

}


fn print_color_sets<CSS: ColorSetStorage>(index: &PseudoalignmentIndex<'_, CSS>, query_path: &Path, print_kmers: bool, out: Output) {
    match out {
        Output::File(w) => print_color_sets_impl(index, query_path, print_kmers, w),
        Output::Stdout(w) => print_color_sets_impl(index, query_path, print_kmers, w),
    }
}

fn print_color_sets_impl<CSS: ColorSetStorage, W: Write>(index: &PseudoalignmentIndex<'_, CSS>, query_path: &Path, print_kmers: bool, mut out: W) {
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

fn open_query_reader(query_path: Option<&Path>) -> jseqio::reader::DynamicFastXReader {
    match query_path {
        Some(p) => jseqio::reader::DynamicFastXReader::from_file(&p).unwrap(),
        None => jseqio::reader::DynamicFastXReader::from_stdin().unwrap(),
    }
}

#[allow(clippy::too_many_arguments)]
fn intersection_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &PseudoalignmentIndex<'_, CSS>, query_path: Option<&Path>, min_hits: usize, metrics: &[pseudoalignment_metrics::Metric], n_threads: usize, out: Output, sort_output: bool, output_settings: pseudoalignment::OutputSettings) {
    let create_new_aligner = move || {
        let aligner = pseudoalignment::IntersectionPseudoaligner::new(min_hits);
        Box::new(aligner) as Box<dyn pseudoalignment::Pseudoaligner<CSS> + Send>
    };

    let reader = open_query_reader(query_path);
    match out {
        Output::File(w) => pseudoalignment::run_pseudoalignment(index, reader, w, create_new_aligner, metrics, n_threads, sort_output, output_settings),
        Output::Stdout(w) => pseudoalignment::run_pseudoalignment(index, reader, w, create_new_aligner, metrics, n_threads, sort_output, output_settings),
    }
    log::info!("Finished");
}

#[allow(clippy::too_many_arguments)]
fn threshold_pseudoalignment<CSS: ColorSetStorage + Send + Sync>(index: &PseudoalignmentIndex<'_, CSS>, query_path: Option<&Path>, min_hits: usize, threshold: f64, denominator: Denominator, metrics: &[pseudoalignment_metrics::Metric], n_threads: usize, out: Output, sort_output: bool, output_settings: pseudoalignment::OutputSettings) {
    // Map from the CLI denominator enum to the pseudoalignment denominator enum.
    let denominator = match denominator {
        Denominator::All => pseudoalignment::Denominator::All,
        Denominator::Relevant => pseudoalignment::Denominator::Relevant,
        Denominator::MaxHits => pseudoalignment::Denominator::MaxHits,
    };

    let n_colors = index.get_set_storage().n_colors();

    let query_log_string = query_path.map(|p| p.display().to_string()).unwrap_or_else(|| "stdin".to_string());
    log::info!("Running threshold pseudoalignment for query sequences in {}", query_log_string);
    let create_new_aligner = move || {
        let aligner = pseudoalignment::ThresholdPseudoaligner::new(
            n_colors,
            threshold,
            min_hits,
            denominator
        );
        Box::new(aligner) as Box<dyn pseudoalignment::Pseudoaligner<CSS> + Send>
    };

    let reader = open_query_reader(query_path);
    match out {
        Output::File(w) => pseudoalignment::run_pseudoalignment(index, reader, w, create_new_aligner, metrics, n_threads, sort_output, output_settings),
        Output::Stdout(w) => pseudoalignment::run_pseudoalignment(index, reader, w, create_new_aligner, metrics, n_threads, sort_output, output_settings),
    }
    log::info!("Finished");
}

#[allow(clippy::too_many_arguments)]
fn run_merge_tree(infiles: &[PathBuf], temp_dir: &Path, outfile: &Path, n_threads: usize, merge_shared_colors: bool, low_ram_mode: bool, keep_redundant_dummies: bool, sample_distance: usize) {
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
                    temp_dir.join(format!("merge_round{}_{}", round, next_files.len()))
                };
                // Only the last round removes the redundant dummy nodes: it also removes those kept
                // in the inputs, so the result is the same, without a cleanup pass in every round
                let keep_dummies = keep_redundant_dummies || round != n_rounds - 1;
                log::info!("Merging {} and {} into {}", pair[0].display(), pair[1].display(), outpath.display());
                let writer = create_index_writer(&outpath);
                let orientation1 = open_index(&pair[0], &DBG_COMPONENTS).metadata().both_orientations;
                let orientation2 = open_index(&pair[1], &DBG_COMPONENTS).metadata().both_orientations;
                let c1 = load_dbg_index(&pair[0]);
                let c2 = load_dbg_index(&pair[1]);
                let (merged, merged_lcs) = merge::merge_compact_colex_kmers(c1, c2, merge_shared_colors, low_ram_mode, keep_dummies, sample_distance, n_threads);
                write_index(writer, &merged, &merged_lcs, combined_orientation(orientation1, orientation2));
                next_files.push(outpath);
            } else {
                next_files.push(pair[0].clone());
            }
        }
        current_files = next_files;
    }
}

// Dispatch with monomorphized output streams
fn export_index_dispatch<CSS: ColorSetStorage + Sync + Send>(
    index: DbgIndex<CSS>, 
    unitigs_out: Option<impl Write + Sync + Send>,
    metadata_out: impl Write + Sync + Send,
    colors_out: Option<impl Write + Sync + Send>,
    unitigs_as_gfa: bool,
    n_threads: usize
) {

    // Destructure the index, so that the DBG can take over the marks without copying them
    let (sbwt, sufgroups, dummy_marks, coloring) = index.into_parts();
    let (map, sets, _color_names) = coloring.into_parts();

    if let Some(mut colors_out) = colors_out {
        log::info!("Exporting color sets");
        unitig_export::write_color_sets(&mut colors_out, &sets);
    } else {
        log::info!("Skipped exporting color sets");
    }

    let dbg = Dbg::from_marks(&sbwt, sufgroups, dummy_marks);

    let n_unitigs: Option<usize> = if let Some(unitigs_out) = unitigs_out {
        if unitigs_as_gfa {
            log::info!("Exporting unitigs (format: GFA)");
            Some(unitig_export::export_gfa(&sbwt, &dbg, &map, unitigs_out, n_threads))
        } else {
            log::info!("Exporting unitigs (format: fasta)");
            Some(unitig_export::export_colored_unitigs(&sbwt, &dbg, &map, unitigs_out, n_threads))
        }
    } else { 
        log::info!("Skipped exporting unitigs");
        None 
    };

    log::info!("Exporting metadata");
    unitig_export::write_metadata(metadata_out, n_unitigs, &sets, sbwt.k());

}

fn export_index<CSS: ColorSetStorage + Sync + Send>(
    index: DbgIndex<CSS>, 
    out_prefix: &Path,
    colors_to_stdout: bool,
    unitigs_as_gfa: bool,
    no_unitigs: bool,
    no_color_sets: bool,
    n_threads: usize
) {

    let out_prefix = out_prefix.as_os_str().to_str().unwrap().to_owned();

    let mut metadata_filename = out_prefix.clone();
    metadata_filename.push_str(".metadata.txt");

    let mut unitig_filename = out_prefix.clone();
    if unitigs_as_gfa {
        unitig_filename.push_str(".unitigs.gfa");
    } else {
        unitig_filename.push_str(".unitigs.fa");
    }

    let mut colors_filename = out_prefix.clone();
    colors_filename.push_str(".color_sets.txt");

    let metadata_out = BufWriter::new(File::create(metadata_filename).unwrap());

    let unitigs_out: Option<BufWriter<File>> = if no_unitigs {
        None    
    } else {
        Some(BufWriter::new(File::create(unitig_filename).unwrap()))
    };

    if colors_to_stdout {
        let colors_out: Option<BufWriter<Stdout>> = if no_color_sets {
            None
        } else {
            Some(BufWriter::new(std::io::stdout()))
        };

        export_index_dispatch(index, unitigs_out, metadata_out, colors_out, unitigs_as_gfa, n_threads);
    } else {

        let colors_out: Option<BufWriter<File>> = if no_color_sets {
            None
        } else {
            Some(BufWriter::new(File::create(colors_filename).unwrap()))
        };

        export_index_dispatch(index, unitigs_out, metadata_out, colors_out, unitigs_as_gfa, n_threads);
    }
}

fn get_sbwt_and_lcs(sbwt_path: &Option<PathBuf>, lcs_path: &Option<PathBuf>, temp_dir: &Option<PathBuf>, input_stream: io::ChainedInputStream, k: usize, n_threads: usize) -> (SbwtIndex<SubsetMatrix>, LcsArray) {
    let (sbwt, lcs) = if let Some(sbwt_path) = sbwt_path {
        log::info!("Loading SBWT from {}", sbwt_path.display());
        let mut sbwt_in = BufReader::new(File::open(sbwt_path).unwrap());
        // TODO: support other SBWT variants (e.g. SubsetCorrectionSets) instead of only SubsetMatrix.
        // Themisto is hardcoded to SbwtIndex<SubsetMatrix> throughout, so this needs generics or a conversion.
        let SbwtIndexVariant::SubsetMatrix(mut sbwt) = SbwtIndexVariant::load(&mut sbwt_in).unwrap() else {
            log::error!("Only SBWT indexes with the SubsetMatrix subset rank structure are supported");
            std::process::exit(1);
        };

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
        log::info!("SBWT not provided -> building the SBWT.");
        let temp_dir = temp_dir.as_ref().expect("Tempory directory not specified (must be specified for SBWT construction)");
        let (mut sbwt, lcs) = BitPackedKmerSortingDisk::new(input_stream, k)
            .add_rev_comp(true)
            .build_lcs(true)
            .n_threads(n_threads)
            .dedup_batches(true)
            .temp_dir(temp_dir)
            .run();
        log::info!("Building SBWT select support");
        sbwt.build_select();
        let sbwt = sbwt;
        let lcs = lcs.unwrap();
        (sbwt, lcs)
    };

    (sbwt, lcs)

}

fn into_metric_list(report_hit_counts: bool, report_bases_covered: bool) -> Vec<pseudoalignment_metrics::Metric> {
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
fn resolve_output_settings(themisto1_output_format: bool, positive_reads_only: bool, metrics: &[pseudoalignment_metrics::Metric]) -> pseudoalignment::OutputSettings {
    let format = if themisto1_output_format {
        if !metrics.is_empty() {
            log::error!("The Themisto 1 output format does not support reporting pseudoalignment metrics. Drop --themisto1-output-format or the metric reporting flags.");
            std::process::exit(1);
        }
        pseudoalignment::OutputFormat::Themisto1
    } else {
        pseudoalignment::OutputFormat::Jsonl
    };
    pseudoalignment::OutputSettings { format, report_only_reads_that_hit_something: positive_reads_only}
}

fn read_file_of_paths(file: &PathBuf) -> Vec<PathBuf> {
    BufReader::new(File::open(file).unwrap())
        .lines()
        .map(|s| PathBuf::from_str(&s.unwrap()).unwrap())
        .collect()
}

// Returns pairs (input file, output file). The input file can be None, in which
// case the input should be read from stdin. The output file can be None, in which
// case the output should be written to stdout.
fn get_input_and_output_files(query_path: Option<PathBuf>, output_path: Option<PathBuf>, query_listfile: Option<PathBuf>, output_listfile: Option<PathBuf>)
-> Vec<(Option<PathBuf>, Option<PathBuf>)>{
    let mut io_pairs: Vec<(Option<PathBuf>, Option<PathBuf>)> = vec![];

    let no_input_source_given = query_path.is_none() && query_listfile.is_none();

    if let Some(f) = query_path {
        io_pairs.push((Some(f), output_path.clone()))
    }

    if let Some(query_listfile) = query_listfile {
        let query_files = read_file_of_paths(&query_listfile);
        if query_files.is_empty() {
            eprintln!("Error: --query-list file {} is empty", query_listfile.display());
            std::process::exit(1);
        }
        let output_files: Vec<Option<PathBuf>> = if let Some(output_listfile) = output_listfile {
            read_file_of_paths(&output_listfile).into_iter().map(Some).collect()
        } else {
            vec![None; query_files.len()]
        };

        assert!(query_files.len() == output_files.len(), "Error: number of input files ({}) does not match the number of output files ({})", query_files.len(), output_files.len());

        for i in 0..query_files.len() {
            io_pairs.push((Some(query_files[i].clone()), output_files[i].clone()));
        }
    }

    if no_input_source_given {
        // Default: read query sequences from stdin
        io_pairs.push((None, output_path));
    }

    io_pairs
}

fn main() -> std::process::ExitCode {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { color_fof, sequence_colors_file, seq_colors_by_name, output, temp_dir, k, n_threads, sample_distance, sbwt_path, lcs_path, from_unitigs, n_pieces, n_parser_threads} => {
            if k % 2 == 0 && from_unitigs {
                panic!("--from_unitigs requires odd k");
            }

            let (input_mode, all_input_paths, color_names) = match (color_fof, sequence_colors_file) {
                (None, None) => panic!("Must give one of --file-colors or --seq-colors"),
                (Some(_), Some(_)) => todo!("Must not give both --file-colors and --seq-colors"),
                (Some(color_fof), None) => {
                    let (file_groups, color_names) = parse_color_fof(BufReader::new(File::open(color_fof).unwrap())).unwrap_or_else(|e| {
                        log::error!("{}", e);
                        panic!("{}", e);
                    });
                    if from_unitigs {
                        let n_multi = file_groups.iter().filter(|group| group.len() > 1).count();
                        if n_multi > 0 {
                            log::warn!("{} colors have more than one unitig file. With --from-unitigs, no k-mer may occur more than once within a color (counting both orientations), which is not checked. Files from separate unitig runs usually share k-mers, which would make the index silently wrong.", n_multi);
                        }
                    }
                    let input_paths: Vec<PathBuf> = file_groups.concat();
                    (ColoredSeqInput::FileColors(file_groups), input_paths, color_names)
                }
                (None, Some(sequence_colors_file)) => {
                    // Read color names from the sequence file
                    log::info!("Reading all sequence names from the input file");
                    let mut seq_names = Vec::<String>::new();
                    let mut reader = needletail::parse_fastx_file(&sequence_colors_file).unwrap();
                    while let Some(rec) = reader.next() {
                        let rec = rec.unwrap();
                        let id = rec.id();
                        let name = id.split(|&b| b == b' ').next().unwrap_or(id);
                        seq_names.push(String::from_utf8(name.to_owned()).unwrap());
                    }
                    log::info!("Read {} sequences", seq_names.len());
                    let n_seqs = seq_names.len();
                    let (record_colors, color_names) = assign_seq_colors(seq_names, seq_colors_by_name).unwrap_or_else(|e| {
                        log::error!("{}", e);
                        panic!("{}", e);
                    });
                    if seq_colors_by_name {
                        log::info!("Grouped {} sequences into {} colors by name", n_seqs, color_names.len());
                    }
                    (ColoredSeqInput::SequenceColors(sequence_colors_file.clone(), record_colors), vec![sequence_colors_file], color_names)
                },
            };

            let input_stream = io::ChainedInputStream::new(all_input_paths.clone());

            match n_pieces {
                None => {
                    // Create the output directory first, to fail before the work if it can't be
                    let writer = create_index_writer(&output);
                    let (sbwt, lcs) = get_sbwt_and_lcs(&sbwt_path, &lcs_path, &temp_dir, input_stream, k, n_threads);
                    let index = build_coloring::<SparseDenseStorage>(sbwt, lcs, input_mode, &color_names, n_threads, sample_distance, from_unitigs, BuildMode::InMemory, n_parser_threads);
                    let (index, lcs) = index.unwrap(); // Ok because we passed in InMemory
                    // The build adds the reverse complements of the input and colors both orientations
                    write_index(writer, &index, &lcs, Some(true));
                },
                Some(n_pieces) => {
                    // Check that the output file can be created
                    let _out_test = File::create(&output).unwrap();
                    std::fs::remove_file(&output).unwrap();
                    let (sbwt, lcs) = get_sbwt_and_lcs(&sbwt_path, &lcs_path, &temp_dir, input_stream, k, n_threads);
                    build_coloring::<SparseDenseStorage>(sbwt, lcs, input_mode, &color_names, n_threads, sample_distance, from_unitigs, BuildMode::ToDisk(output, n_pieces), n_parser_threads);
                }
            }
        },

        Subcommands::IntersectionPseudoalign { index: index_path, query: query_path, query_list, query_output_list, min_hits, n_threads, report_bases_covered, report_hit_counts, sort_output, themisto1_output_format, output, positive_reads_only} => {
            let (sbwt, sufgroups, coloring) = load_for_pseudoalignment(&index_path);
            let index = PseudoalignmentIndex::new(&sbwt, &sufgroups, &coloring);
            let metrics = into_metric_list(report_hit_counts, report_bases_covered);
            let output_settings = resolve_output_settings(themisto1_output_format, positive_reads_only, &metrics);

            let io_file_pairs = get_input_and_output_files(query_path, output, query_list, query_output_list);

            for (infile, outfile) in io_file_pairs {
                let input_log_string = infile.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdin".to_string());
                let output_log_string = outfile.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdout".to_string());
                log::info!("Running intersection pseudoalignment for query sequences in {}, writing output to {}", input_log_string, output_log_string);

                let out = open_output(outfile);
                intersection_pseudoalignment(&index, infile.as_deref(), min_hits, &metrics, n_threads, out, sort_output, output_settings);

            }

        },
        Subcommands::ThresholdPseudoalign { index: index_path, query: query_path, query_list, query_output_list, min_hits, threshold, denominator, n_threads, report_hit_counts, report_bases_covered, sort_output, themisto1_output_format, output, positive_reads_only} => {
            let (sbwt, sufgroups, coloring) = load_for_pseudoalignment(&index_path);
            let index = PseudoalignmentIndex::new(&sbwt, &sufgroups, &coloring);
            let metrics = into_metric_list(report_hit_counts, report_bases_covered);
            let output_settings = resolve_output_settings(themisto1_output_format, positive_reads_only, &metrics);
            let io_file_pairs = get_input_and_output_files(query_path, output, query_list, query_output_list);

            for (infile, outfile) in io_file_pairs {
                let input_log_string = infile.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdin".to_string());
                let output_log_string = outfile.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdout".to_string());
                log::info!("Running threshold pseudoalignment for query sequences in {}, writing output to {}", input_log_string, output_log_string);

                let out = open_output(outfile);
                threshold_pseudoalignment(&index, infile.as_deref(), min_hits, threshold, denominator, &metrics, n_threads, out, sort_output, output_settings);
            }
        },
        Subcommands::PrintColorSets { index: index_path, query: query_path, print_kmers, output } => {
            let (sbwt, sufgroups, coloring) = load_for_pseudoalignment(&index_path);
            let index = PseudoalignmentIndex::new(&sbwt, &sufgroups, &coloring);
            let out = open_output(output);
            print_color_sets(&index, &query_path, print_kmers, out);
        },
        Subcommands::DumpColorNames{ index: index_path, output } => {
            let color_names = or_exit(open_index(&index_path, &[Component::ColorNames]).load_color_names());
            let out = open_output(output);
            log::info!("Dumping color names");
            print_color_names(&color_names, out);
        },
        Subcommands::Merge { index_file_list, temp_dir, outfile, n_threads, low_ram_mode, merge_shared_colors, keep_redundant_dummies, sample_distance } => {
            let infiles: Vec<PathBuf> = BufReader::new(File::open(index_file_list).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            run_merge_tree(&infiles, &temp_dir, &outfile, n_threads, merge_shared_colors, low_ram_mode, keep_redundant_dummies, sample_distance);
        },
        Subcommands::Intersect { index1, index2, outfile, colors, merge_shared_colors, sample_distance, n_threads, low_ram_mode } => {
            if merge_shared_colors && colors != intersect::IntersectColors::Union {
                log::warn!("--merge-shared-colors has no effect with --colors {:?}: colors are always matched by name", colors);
            }
            let writer = create_index_writer(&outfile);
            let orientation = combined_orientation(open_index(&index1, &DBG_COMPONENTS).metadata().both_orientations, open_index(&index2, &DBG_COMPONENTS).metadata().both_orientations);
            let c1 = load_dbg_index(&index1);
            let c2 = load_dbg_index(&index2);
            log::info!("Intersecting {} and {}", index1.display(), index2.display());
            let (result, result_lcs) = intersect::intersect_compact_colex_kmers(c1, c2, colors, merge_shared_colors, low_ram_mode, sample_distance, n_threads);
            write_index(writer, &result, &result_lcs, orientation);
        },
        Subcommands::Subset { index, outfile, colors, colors_file, complement, keep_all_colors, sample_distance, n_threads, low_ram_mode } => {
            let mut names = colors;
            if let Some(colors_file) = colors_file {
                let content = std::fs::read_to_string(&colors_file).unwrap_or_else(|e| {
                    log::error!("Could not read {}: {}", colors_file.display(), e);
                    std::process::exit(1);
                });
                names.extend(content.lines().map(|line| line.trim().to_string()).filter(|line| !line.is_empty()));
            }
            if names.is_empty() {
                log::error!("No colors given: use --colors or --colors-file");
                std::process::exit(1);
            }
            {
                {
                    let writer = create_index_writer(&outfile);
                    let index_dir = open_index(&index, &DBG_COMPONENTS);
                    let orientation = index_dir.metadata().both_orientations;
                    let coloring = or_exit(index_dir.load_dbg_index());
                    let color_names = coloring.get_color_names();
                    // Sets of the names, so that matching them takes time linear in the number of
                    // names: with lists, keeping 108,930 of the 121,034 colors of an index took 46 s
                    let index_names: std::collections::HashSet<&str> = color_names.iter().map(|name| name.as_str()).collect();
                    if let Some(unknown) = names.iter().find(|name| !index_names.contains(name.as_str())) {
                        log::error!("The index has no color named {:?}", unknown);
                        std::process::exit(1);
                    }
                    let kept_names: std::collections::HashSet<&str> = names.iter().map(|name| name.as_str()).collect();
                    let kept_colors: Vec<bool> = color_names.iter().map(|name| kept_names.contains(name.as_str()) != complement).collect();
                    log::info!("Taking the subset of {} for {} of its {} colors{}", index.display(), kept_colors.iter().filter(|&&kept| kept).count(), kept_colors.len(), if keep_all_colors { ", keeping all colors of its k-mers" } else { "" });
                    let (result, result_lcs) = subset::subset_compact_colex_kmers(coloring, &kept_colors, keep_all_colors, low_ram_mode, sample_distance, n_threads);
                    write_index(writer, &result, &result_lcs, orientation);
                }
            }
        },
        Subcommands::Difference { index1, index2, outfile, sample_distance, n_threads, low_ram_mode } => {
            let writer = create_index_writer(&outfile);
            let index_dir1 = open_index(&index1, &DBG_COMPONENTS);
            // Only the SBWT of index2 is needed: the result has the colors of index1
            let index_dir2 = open_index(&index2, &[Component::Sbwt]);
            let orientation = combined_orientation(index_dir1.metadata().both_orientations, index_dir2.metadata().both_orientations);
            let c1 = or_exit(index_dir1.load_dbg_index());
            let sbwt2 = or_exit(index_dir2.load_sbwt(false));
            log::info!("Subtracting {} from {}", index2.display(), index1.display());
            let (result, result_lcs) = difference::difference_compact_colex_kmers(c1, sbwt2, low_ram_mode, sample_distance, n_threads);
            write_index(writer, &result, &result_lcs, orientation);
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
            let (sbwt,lcs) = get_sbwt_and_lcs(&sbwt_path, &lcs_path, &temp_dir, input_stream, metadata.k, n_threads);

            if sbwt.k() != metadata.k {
                log::error!("SBWT k does not match the index dump k ({} vs {})", sbwt.k(), metadata.k);
                return ExitCode::FAILURE;
            }

            let writer = create_index_writer(&out_path);

            let unitig_dump = BufReader::new(File::open(&unitig_filename).unwrap());
            let color_dump = BufReader::new(File::open(&color_sets_filename).unwrap());
            let metadata_dump = BufReader::new(File::open(&metadata_filename).unwrap());

            let coloring = Coloring::<SparseDenseStorage>::new_from_colored_unitig_dump(
                &sbwt, &lcs, sample_distance, n_threads, metadata_dump, unitig_dump, color_dump);
            let index = DbgIndex::from_lcs(sbwt, &lcs, coloring, n_threads);
            // The import colors both orientations of every unitig
            write_index(writer, &index, &lcs, Some(true));
        },
        Subcommands::ImportSbwt { sbwt_path, lcs_path, output, sample_distance, n_threads, color_name } => {
            let writer = create_index_writer(&output);
            log::info!("Loading SBWT from {}", sbwt_path.display());
            let mut sbwt_in = BufReader::new(File::open(&sbwt_path).unwrap());
            // TODO: support other SBWT variants, see get_sbwt_and_lcs
            let SbwtIndexVariant::SubsetMatrix(mut sbwt) = SbwtIndexVariant::load(&mut sbwt_in).unwrap() else {
                log::error!("Only SBWT indexes with the SubsetMatrix subset rank structure are supported");
                std::process::exit(1);
            };

            log::info!("Building select support for SBWT");
            sbwt.build_select();

            let lcs = if let Some(lcs_path) = lcs_path {
                log::info!("Loading the LCS array from {}", lcs_path.display());
                LcsArray::load(&mut BufReader::new(File::open(&lcs_path).unwrap())).unwrap()
            } else {
                log::info!("Building LCS array");
                LcsArray::from_sbwt(&sbwt, n_threads, true)
            };

            log::info!("Building single-colored Themisto 2 index");
            let coloring = Coloring::<SparseDenseStorage>::new_single_colored(&sbwt, &lcs, sample_distance, n_threads, color_name);
            let index = DbgIndex::from_lcs(sbwt, &lcs, coloring, n_threads);
            // Whether the SBWT has both orientations of its k-mers is not known
            write_index(writer, &index, &lcs, None);
        },
        Subcommands::Export { index: index_path, color_dump_prefix, n_threads, gfa, colors_to_stdout, no_unitigs, no_color_sets } => {
            let index = load_dbg_index(&index_path); // Select support is required for export
            export_index(index, &color_dump_prefix, colors_to_stdout, gfa, no_unitigs, no_color_sets, n_threads);
        },
        Subcommands::Stats { index: index_path, n_threads } => {
            let index_dir = open_index(&index_path, &DBG_COMPONENTS);
            let index = or_exit(index_dir.load_dbg_index());
            print_stats(&index, &index_dir, n_threads);
        }
        Subcommands::Report { index: index_path, input, output } => {
            let color_names = or_exit(open_index(&index_path, &[Component::ColorNames]).load_color_names());
            log::info!("Processing pseudoalignment file");
            match (input, output) {
                (Some(in_path), Some(out_path)) => report(BufReader::new(File::open(in_path).unwrap()), BufWriter::new(File::create(out_path).unwrap()), &color_names),
                (Some(in_path), None) => report(BufReader::new(File::open(in_path).unwrap()), BufWriter::new(std::io::stdout()), &color_names),
                (None, Some(out_path)) => report(BufReader::new(std::io::stdin()), BufWriter::new(File::create(out_path).unwrap()), &color_names),
                (None, None) => report(BufReader::new(std::io::stdin()), BufWriter::new(std::io::stdout()), &color_names),
            }
        },
        Subcommands::FilterReads { pseudoalignment, reads, output, color_ids } => {
            return filter_reads::filter_reads(&pseudoalignment, &reads, &output, &color_ids);
        },
        Subcommands::FinimizerStats { index: index_path, n_threads} => {
            // Hidden subcommand, might be deleted at any moment
            let index_dir = open_index(&index_path, &[Component::Sbwt, Component::Lcs, Component::ColorNames, Component::ColorSampling, Component::ColorStorage]);
            let sbwt = or_exit(index_dir.load_sbwt(true)); // Select support is required for verify
            let lcs = or_exit(index_dir.load_lcs());
            let coloring = or_exit(index_dir.load_coloring());
            finimizers::minimizer_stats(&ColoredKmers::new(&sbwt, &coloring), &lcs, n_threads, finimizers::MinimizerType::Finimizer);
        },
        Subcommands::MinimizerStats { index: index_path, n_threads, m} => {
            // Hidden subcommand, might be deleted at any moment
            let index_dir = open_index(&index_path, &[Component::Sbwt, Component::Lcs, Component::ColorNames, Component::ColorSampling, Component::ColorStorage]);
            let sbwt = or_exit(index_dir.load_sbwt(true)); // Select support is required for verify
            let lcs = or_exit(index_dir.load_lcs());
            let coloring = or_exit(index_dir.load_coloring());
            finimizers::minimizer_stats(&ColoredKmers::new(&sbwt, &coloring), &lcs, n_threads, finimizers::MinimizerType::Minimizer(m));
        },
    }

    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    #[test]
    fn test_parse_color_fof() {
        let fof = "a.fna\nb.fna\tX\n\nc.fna\tY\nd.fna\tX\na.fna\n";
        let (groups, names) = super::parse_color_fof(fof.as_bytes()).unwrap();
        assert_eq!(names, vec!["a.fna", "X", "Y"]);
        let paths = |v: &[&str]| -> Vec<PathBuf> { v.iter().map(PathBuf::from).collect() };
        assert_eq!(groups, vec![paths(&["a.fna", "a.fna"]), paths(&["b.fna", "d.fna"]), paths(&["c.fna"])]);
    }

    fn names(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_assign_seq_colors_duplicate_names_collapse() {
        let (colors, color_names) = super::assign_seq_colors(names(&["chr1", "chr1", "chr1", "chr2", "chr2"]), true).unwrap();
        assert_eq!(colors, vec![0, 0, 0, 1, 1]);
        assert_eq!(color_names, names(&["chr1", "chr2"]));
    }

    #[test]
    fn test_assign_seq_colors_distinct_names() {
        // With distinct names, grouping by name gives the same colors as one color per sequence
        for by_name in [false, true] {
            let (colors, color_names) = super::assign_seq_colors(names(&["a", "b", "c"]), by_name).unwrap();
            assert_eq!(colors, vec![0, 1, 2]);
            assert_eq!(color_names, names(&["a", "b", "c"]));
        }
    }

    #[test]
    fn test_assign_seq_colors_mixed() {
        let (colors, color_names) = super::assign_seq_colors(names(&["x", "x", "y", "z", "z", "z", "w"]), true).unwrap();
        assert_eq!(colors, vec![0, 0, 1, 2, 2, 2, 3]);
        assert_eq!(color_names, names(&["x", "y", "z", "w"]));
    }

    #[test]
    fn test_assign_seq_colors_rejects_non_consecutive_names() {
        let e = super::assign_seq_colors(names(&["x", "y", "x"]), true).unwrap_err();
        assert!(e.contains("Sequence 3") && e.contains("not consecutive"), "{}", e);
        // Without grouping, the order of the names does not matter
        assert!(super::assign_seq_colors(names(&["x", "y", "x"]), false).is_ok());
    }

    #[test]
    fn test_assign_seq_colors_without_by_name_keeps_one_color_per_sequence() {
        let (colors, color_names) = super::assign_seq_colors(names(&["chr1", "chr1", "chr2"]), false).unwrap();
        assert_eq!(colors, vec![0, 1, 2]);
        assert_eq!(color_names, names(&["chr1", "chr1", "chr2"]));
    }

    #[test]
    fn test_assign_seq_colors_single_sequence() {
        for by_name in [false, true] {
            let (colors, color_names) = super::assign_seq_colors(names(&["only"]), by_name).unwrap();
            assert_eq!(colors, vec![0]);
            assert_eq!(color_names, names(&["only"]));
        }
    }

    #[test]
    fn test_assign_seq_colors_all_same_name() {
        let (colors, color_names) = super::assign_seq_colors(names(&["s", "s", "s", "s"]), true).unwrap();
        assert_eq!(colors, vec![0, 0, 0, 0]);
        assert_eq!(color_names, names(&["s"]));
    }

    #[test]
    fn test_assign_seq_colors_empty_name() {
        // Grouping would silently merge all unnamed sequences into one color, so it is an error
        assert!(super::assign_seq_colors(names(&["a", "", "a"]), true).is_err());
        // Without grouping, unnamed sequences are still accepted, one color each
        let (colors, color_names) = super::assign_seq_colors(names(&["", ""]), false).unwrap();
        assert_eq!(colors, vec![0, 1]);
        assert_eq!(color_names, names(&["", ""]));
    }

    #[test]
    fn test_parse_color_fof_rejects_empty_fields() {
        assert!(super::parse_color_fof("a.fna\t\n".as_bytes()).is_err());
        assert!(super::parse_color_fof("\tX\n".as_bytes()).is_err());
    }
}
