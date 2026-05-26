use std::{fs::File, io::{BufRead, BufReader}, path::Path, process::ExitCode};

use jseqio::record::Record;
use jseqio::writer::SeqRecordWriter;
use rustc_hash::FxHashMap;

#[derive(serde::Deserialize)]
struct PseudoalignmentRecord<'a> {
    #[serde(borrow)]
    name: &'a str,
    colors: Vec<usize>,
}

fn read_color_ids(path: &Path) -> Vec<usize> {
    let mut ids = Vec::<usize>::new();
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
        ids.push(id);
    }
    ids
}

// A `Sorted` value is a slice that is statically guaranteed to be sorted ascending,
// because the only way to construct one is via `Sorted::sort`, which sorts the underlying slice.
#[derive(Clone, Copy)]
struct Sorted<'a>(&'a [usize]);

impl<'a> Sorted<'a> {
    fn sort(v: &'a mut [usize]) -> Self {
        v.sort_unstable();
        Self(v)
    }
}

// Linear merge-style intersection check. Both inputs are statically guaranteed to be sorted ascending.
fn sorted_lists_intersect(a: Sorted, b: Sorted) -> bool {
    let (a, b) = (a.0, b.0);
    let (mut i, mut j) = (0_usize, 0_usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => return true,
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
    false
}

struct KeepEntry {
    keep: bool,
    seen: bool,
}

// Returns name -> KeepEntry. `keep` is true iff the read's color set intersects target_colors.
// `seen` starts false and is flipped to true during the read-streaming pass.
fn build_keep_map(pseudoalignment_path: &Path, target_colors: Sorted) -> FxHashMap<Vec<u8>, KeepEntry> {
    let reader = BufReader::new(File::open(pseudoalignment_path).unwrap_or_else(|e| {
        log::error!("Could not open pseudoalignment file {}: {}", pseudoalignment_path.display(), e);
        std::process::exit(1);
    }));

    let mut keep: FxHashMap<Vec<u8>, KeepEntry> = FxHashMap::default();
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
        let mut rec: PseudoalignmentRecord = serde_json::from_str(&line).unwrap_or_else(|e| {
            log::error!("Failed to parse JSON in {} on line {}: {}", pseudoalignment_path.display(), line_no + 1, e);
            std::process::exit(1);
        });
        let rec_colors = Sorted::sort(&mut rec.colors);
        let should_keep = sorted_lists_intersect(rec_colors, target_colors);
        if should_keep { n_keep += 1; }
        let prev = keep.insert(rec.name.as_bytes().to_vec(), KeepEntry { keep: should_keep, seen: false });
        if prev.is_some() {
            log::warn!("Duplicate read name in pseudoalignment file {}: {}", pseudoalignment_path.display(), rec.name);
        }
        line_no += 1;
    }
    log::info!("Read {} pseudoalignment records ({} matching the target colors)", keep.len(), n_keep);
    keep
}

pub fn filter_reads(pseudoalignment_path: &Path, reads_path: &Path, output_path: &Path, color_ids_path: &Path) -> ExitCode {
    log::info!("Loading target color ids from {}", color_ids_path.display());
    let mut target_colors_buf = read_color_ids(color_ids_path);
    log::info!("Loaded {} target color id(s)", target_colors_buf.len());
    let target_colors = Sorted::sort(&mut target_colors_buf);

    log::info!("Scanning pseudoalignment file {}", pseudoalignment_path.display());
    let mut keep = build_keep_map(pseudoalignment_path, target_colors);

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
        match keep.get_mut(name) {
            Some(entry) => {
                if entry.seen {
                    log::warn!("Duplicate read name {:?} in {}", String::from_utf8_lossy(name), reads_path.display());
                }
                entry.seen = true;
                if entry.keep {
                    writer.write(&rec).unwrap();
                    n_written += 1;
                }
            }
            None => {
                log::error!("Read name {:?} in {} has no record in pseudoalignment file {}",
                            String::from_utf8_lossy(name), reads_path.display(), pseudoalignment_path.display());
                std::process::exit(1);
            }
        }
    }
    writer.flush().unwrap();

    let n_unseen = keep.values().filter(|e| !e.seen).count();
    if n_unseen > 0 {
        let example: &[u8] = keep.iter().find(|(_, e)| !e.seen).map(|(k, _)| k.as_slice()).unwrap();
        log::error!("{} pseudoalignment record(s) in {} have no matching read in {} (e.g. {:?})",
                    n_unseen, pseudoalignment_path.display(), reads_path.display(), String::from_utf8_lossy(example));
        std::process::exit(1);
    }

    log::info!("Processed {} reads, wrote {} to {}", n_reads, n_written, output_path.display());
    ExitCode::SUCCESS
}
