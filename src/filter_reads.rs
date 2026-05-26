use std::{fs::File, io::{BufRead, BufReader}, path::Path, process::ExitCode};

use jseqio::record::Record;
use jseqio::writer::SeqRecordWriter;
use rustc_hash::{FxHashMap, FxHashSet};

#[derive(serde::Deserialize)]
struct PseudoalignmentRecord<'a> {
    #[serde(borrow)]
    name: &'a str,
    colors: Vec<usize>,
}

fn read_color_ids(path: &Path) -> FxHashSet<usize> {
    let mut set = FxHashSet::<usize>::default();
    let reader = BufReader::new(File::open(path).unwrap_or_else(|e| {
        log::error!("Could not open color id file {}: {}", path.display(), e);
        std::process::exit(1);
    }));
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let id: usize = trimmed.parse().unwrap_or_else(|_| {
            log::error!("Color id file {} line {}: not a non-negative integer: {:?}",
                        path.display(), line_no + 1, trimmed);
            std::process::exit(1);
        });
        set.insert(id);
    }
    set
}

// Returns name -> keep flag. Keep flag is true iff the read's color set intersects target_colors.
fn build_keep_map(pseudoalignment_path: &Path, target_colors: &FxHashSet<usize>) -> FxHashMap<Vec<u8>, bool> {
    let reader = BufReader::new(File::open(pseudoalignment_path).unwrap_or_else(|e| {
        log::error!("Could not open pseudoalignment file {}: {}", pseudoalignment_path.display(), e);
        std::process::exit(1);
    }));

    let mut keep: FxHashMap<Vec<u8>, bool> = FxHashMap::default();
    let mut n_keep = 0_usize;
    let mut line = String::new();
    let mut line_no = 0_usize;
    let mut input = reader;
    loop {
        line.clear();
        let n = input.read_line(&mut line).unwrap();
        if n == 0 { break; }
        if line.ends_with('\n') { line.pop(); }
        if line.is_empty() {
            line_no += 1;
            continue;
        }
        let rec: PseudoalignmentRecord = serde_json::from_str(&line).unwrap_or_else(|e| {
            log::error!("Failed to parse JSON in {} on line {}: {}", pseudoalignment_path.display(), line_no + 1, e);
            std::process::exit(1);
        });
        let should_keep = rec.colors.iter().any(|c| target_colors.contains(c));
        if should_keep { n_keep += 1; }
        let prev = keep.insert(rec.name.as_bytes().to_vec(), should_keep);
        if prev.is_some() {
            log::error!("Duplicate read name in pseudoalignment file {}: {}", pseudoalignment_path.display(), rec.name);
            std::process::exit(1);
        }
        line_no += 1;
    }
    log::info!("Read {} pseudoalignment records ({} matching the target colors)", keep.len(), n_keep);
    keep
}

pub fn filter_reads(pseudoalignment_path: &Path, reads_path: &Path, output_path: &Path, color_ids_path: &Path) -> ExitCode {
    log::info!("Loading target color ids from {}", color_ids_path.display());
    let target_colors = read_color_ids(color_ids_path);
    log::info!("Loaded {} target color id(s)", target_colors.len());

    log::info!("Scanning pseudoalignment file {}", pseudoalignment_path.display());
    let mut keep = build_keep_map(pseudoalignment_path, &target_colors);

    log::info!("Streaming reads from {} and writing filtered output to {}", reads_path.display(), output_path.display());
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&reads_path).unwrap_or_else(|e| {
        log::error!("Could not open reads file {}: {}", reads_path.display(), e);
        std::process::exit(1);
    });
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_file(&output_path).unwrap_or_else(|e| {
        log::error!("Could not create output file {}: {}", output_path.display(), e);
        std::process::exit(1);
    });

    let mut n_reads = 0_usize;
    let mut n_written = 0_usize;
    while let Some(rec) = reader.read_next().unwrap() {
        n_reads += 1;
        let name = rec.name();
        match keep.remove(name) {
            Some(true) => {
                writer.write(&rec).unwrap();
                n_written += 1;
            }
            Some(false) => {}
            None => {
                log::error!("Read name {:?} in {} has no record in pseudoalignment file {}",
                            String::from_utf8_lossy(name), reads_path.display(), pseudoalignment_path.display());
                std::process::exit(1);
            }
        }
    }
    writer.flush().unwrap();

    if !keep.is_empty() {
        let example: &[u8] = keep.keys().next().unwrap();
        log::error!("{} pseudoalignment record(s) in {} have no matching read in {} (e.g. {:?})",
                    keep.len(), pseudoalignment_path.display(), reads_path.display(), String::from_utf8_lossy(example));
        std::process::exit(1);
    }

    log::info!("Processed {} reads, wrote {} to {}", n_reads, n_written, output_path.display());
    ExitCode::SUCCESS
}
