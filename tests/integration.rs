use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(serde::Deserialize, Debug)]
struct PseudoalignmentRecord {
    #[allow(dead_code)]
    name: String,
    colors: Vec<usize>,
}

const BIN: &str = env!("CARGO_BIN_EXE_themisto2");
const PROJECT_DIR: &str = env!("CARGO_MANIFEST_DIR");

static COUNTER: AtomicU64 = AtomicU64::new(0);

fn tmp_dir() -> PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::SeqCst);
    let dir = PathBuf::from(PROJECT_DIR)
        .join("temp")
        .join(format!("integration-test-{}", id));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn themisto2() -> Command {
    let mut cmd = Command::new(BIN);
    cmd.current_dir(PROJECT_DIR);
    cmd
}

#[test]
fn file_colors_and_seq_colors_produce_same_index() {
    let dir = tmp_dir();
    let tmp1 = dir.join("tmp1");
    let tmp2 = dir.join("tmp2");
    std::fs::create_dir_all(&tmp1).unwrap();
    std::fs::create_dir_all(&tmp2).unwrap();

    let index1 = dir.join("file_colors.thm2");
    let index2 = dir.join("seq_colors.thm2");

    let status = themisto2()
        .args(["build", "--file-colors", "example/fof.txt", "--temp-dir"])
        .arg(&tmp1)
        .args(["-k", "3", "-t", "1", "-o"])
        .arg(&index1)
        .status()
        .unwrap();
    assert!(status.success(), "--file-colors build failed");

    let status = themisto2()
        .args(["build", "--seq-colors", "example/seq-colors.fna", "--temp-dir"])
        .arg(&tmp2)
        .args(["-k", "3", "-t", "1", "-o"])
        .arg(&index2)
        .status()
        .unwrap();
    assert!(status.success(), "--seq-colors build failed");

    let bytes1 = std::fs::read(&index1).unwrap();
    let bytes2 = std::fs::read(&index2).unwrap();
    assert_eq!(bytes1, bytes2, "Indexes built with --file-colors and --seq-colors differ");
}

fn read_jsonl(path: &PathBuf) -> Vec<PseudoalignmentRecord> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l)
            .unwrap_or_else(|e| panic!("Failed to parse JSONL line `{}`: {}", l, e)))
        .collect()
}

#[test]
fn pseudoalign_example_data() {
    let dir = tmp_dir();
    let tmp = dir.join("tmp");
    std::fs::create_dir_all(&tmp).unwrap();

    let index = dir.join("index.thm2");

    let status = themisto2()
        .args(["build", "--file-colors", "example/fof.txt", "--temp-dir"])
        .arg(&tmp)
        .args(["-k", "3", "-t", "1", "-o"])
        .arg(&index)
        .status()
        .unwrap();
    assert!(status.success(), "build failed");

    // Query each input file against the index using both pseudoalignment modes.
    // With --file-colors, color i corresponds to the i-th file in fof.txt
    // (C1.fna -> 0, C2.fna -> 1, C3.fna -> 2). Every k-mer extracted from
    // Ci.fna has color i in its color set, so the intersection of color
    // sets across all k-mers of any sequence in Ci.fna must contain i.
    for (color, query_file) in ["example/C1.fna", "example/C2.fna", "example/C3.fna"]
        .iter()
        .enumerate()
    {
        let intersection_out = dir.join(format!("intersection-c{}.jsonl", color + 1));
        let intersection_args = || {
            let mut cmd = themisto2();
            cmd.args(["intersection-pseudoalign", "-i"])
                .arg(&index)
                .args(["-q", query_file, "-t", "1", "--sort-output"]);
            cmd
        };

        let status = intersection_args()
            .arg("-o")
            .arg(&intersection_out)
            .status()
            .unwrap();
        assert!(status.success(), "intersection-pseudoalign failed for {}", query_file);

        let records = read_jsonl(&intersection_out);
        assert!(!records.is_empty(), "no records for intersection on {}", query_file);
        for rec in &records {
            assert!(
                rec.colors.contains(&color),
                "intersection on {}: expected color {} in {:?}",
                query_file, color, rec.colors,
            );
        }

        // Without -o, the same flags should produce identical output on stdout.
        let stdout = intersection_args().output().unwrap();
        assert!(stdout.status.success(), "intersection-pseudoalign (stdout) failed for {}", query_file);
        assert_eq!(
            stdout.stdout,
            std::fs::read(&intersection_out).unwrap(),
            "intersection-pseudoalign stdout differs from -o output for {}",
            query_file,
        );

        // Without -q, queries are read from stdin; output should match.
        let stdin_out = run_with_stdin(
            themisto2()
                .args(["intersection-pseudoalign", "-i"])
                .arg(&index)
                .args(["-t", "1", "--sort-output"]),
            query_file,
        );
        assert_eq!(
            stdin_out,
            std::fs::read(&intersection_out).unwrap(),
            "intersection-pseudoalign stdin output differs from --query output for {}",
            query_file,
        );

        let threshold_out = dir.join(format!("threshold-c{}.jsonl", color + 1));
        let threshold_args = || {
            let mut cmd = themisto2();
            cmd.args(["threshold-pseudoalign", "-i"])
                .arg(&index)
                .args([
                    "-q", query_file,
                    "-m", "1",
                    "-d", "0.7",
                    "-n", "relevant",
                    "-t", "1",
                    "--sort-output",
                ]);
            cmd
        };

        let status = threshold_args()
            .arg("-o")
            .arg(&threshold_out)
            .status()
            .unwrap();
        assert!(status.success(), "threshold-pseudoalign failed for {}", query_file);

        let records = read_jsonl(&threshold_out);
        assert!(!records.is_empty(), "no records for threshold on {}", query_file);
        for rec in &records {
            assert!(
                rec.colors.contains(&color),
                "threshold on {}: expected color {} in {:?}",
                query_file, color, rec.colors,
            );
        }

        let stdout = threshold_args().output().unwrap();
        assert!(stdout.status.success(), "threshold-pseudoalign (stdout) failed for {}", query_file);
        assert_eq!(
            stdout.stdout,
            std::fs::read(&threshold_out).unwrap(),
            "threshold-pseudoalign stdout differs from -o output for {}",
            query_file,
        );

        let stdin_out = run_with_stdin(
            themisto2()
                .args(["threshold-pseudoalign", "-i"])
                .arg(&index)
                .args([
                    "-m", "1",
                    "-d", "0.7",
                    "-n", "relevant",
                    "-t", "1",
                    "--sort-output",
                ]),
            query_file,
        );
        assert_eq!(
            stdin_out,
            std::fs::read(&threshold_out).unwrap(),
            "threshold-pseudoalign stdin output differs from --query output for {}",
            query_file,
        );
    }
}

fn write_file(path: &PathBuf, content: &str) {
    std::fs::write(path, content).unwrap();
}

fn fastq_record_names(content: &str) -> Vec<&str> {
    // A FASTQ record is four lines (@head, seq, +, qual), so the header line
    // is at indices 0, 4, 8, ... We rely on that structure instead of scanning
    // for `@` because quality strings can start with `@`.
    content.lines().step_by(4).map(|l| l.strip_prefix('@').expect("header missing @")).collect()
}

#[test]
fn filter_reads_happy_path() {
    let dir = tmp_dir();
    let pa = dir.join("pa.jsonl");
    let reads = dir.join("reads.fastq");
    let colors = dir.join("colors.txt");
    let out = dir.join("filtered.fastq");

    write_file(&pa, "{\"name\": \"r1\", \"colors\": [0]}\n\
                     {\"name\": \"r2\", \"colors\": [1]}\n\
                     {\"name\": \"r3\", \"colors\": [2]}\n\
                     {\"name\": \"r4\", \"colors\": []}\n");
    write_file(&reads, "@r1\nACGT\n+\nIIII\n\
                        @r2\nACGT\n+\nIIII\n\
                        @r3\nACGT\n+\nIIII\n\
                        @r4\nACGT\n+\nIIII\n");
    write_file(&colors, "0\n2\n");

    let result = themisto2()
        .args(["filter-reads", "-p"]).arg(&pa)
        .arg("-r").arg(&reads)
        .arg("-o").arg(&out)
        .arg("-c").arg(&colors)
        .output().unwrap();
    assert!(result.status.success(),
            "filter-reads failed: {}", String::from_utf8_lossy(&result.stderr));

    let filtered = std::fs::read_to_string(&out).unwrap();
    assert_eq!(fastq_record_names(&filtered), vec!["r1", "r3"]);
}

#[test]
fn filter_reads_duplicate_read_names() {
    let dir = tmp_dir();
    let pa = dir.join("pa.jsonl");
    let reads = dir.join("reads.fastq");
    let colors = dir.join("colors.txt");
    let out = dir.join("filtered.fastq");

    write_file(&pa, "{\"name\": \"r1\", \"colors\": [0]}\n\
                     {\"name\": \"r2\", \"colors\": [1]}\n");
    // r1 appears twice in the read file; both occurrences should be written and a warning logged.
    write_file(&reads, "@r1\nACGT\n+\nIIII\n\
                        @r2\nACGT\n+\nIIII\n\
                        @r1\nACGT\n+\nIIII\n");
    write_file(&colors, "0\n");

    let result = themisto2()
        .args(["filter-reads", "-p"]).arg(&pa)
        .arg("-r").arg(&reads)
        .arg("-o").arg(&out)
        .arg("-c").arg(&colors)
        .output().unwrap();
    assert!(result.status.success(),
            "filter-reads failed: {}", String::from_utf8_lossy(&result.stderr));

    let filtered = std::fs::read_to_string(&out).unwrap();
    assert_eq!(fastq_record_names(&filtered), vec!["r1", "r1"]);

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("Duplicate read name \"r1\""),
            "expected duplicate-name warning in stderr, got: {}", stderr);
}

#[test]
fn filter_reads_unknown_read_errors() {
    let dir = tmp_dir();
    let pa = dir.join("pa.jsonl");
    let reads = dir.join("reads.fastq");
    let colors = dir.join("colors.txt");
    let out = dir.join("filtered.fastq");

    write_file(&pa, "{\"name\": \"r1\", \"colors\": [0]}\n");
    // r2 has no record in the pseudoalignment file.
    write_file(&reads, "@r1\nACGT\n+\nIIII\n\
                        @r2\nACGT\n+\nIIII\n");
    write_file(&colors, "0\n");

    let result = themisto2()
        .args(["filter-reads", "-p"]).arg(&pa)
        .arg("-r").arg(&reads)
        .arg("-o").arg(&out)
        .arg("-c").arg(&colors)
        .output().unwrap();
    assert!(!result.status.success(), "expected failure but command succeeded");

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("has no record in pseudoalignment file"),
            "expected error about missing pseudoalignment record, got: {}", stderr);
}

#[test]
fn filter_reads_missing_read_errors() {
    let dir = tmp_dir();
    let pa = dir.join("pa.jsonl");
    let reads = dir.join("reads.fastq");
    let colors = dir.join("colors.txt");
    let out = dir.join("filtered.fastq");

    // pseudoalignment has 2 records but the read file only has r1.
    write_file(&pa, "{\"name\": \"r1\", \"colors\": [0]}\n\
                     {\"name\": \"r2\", \"colors\": [0]}\n");
    write_file(&reads, "@r1\nACGT\n+\nIIII\n");
    write_file(&colors, "0\n");

    let result = themisto2()
        .args(["filter-reads", "-p"]).arg(&pa)
        .arg("-r").arg(&reads)
        .arg("-o").arg(&out)
        .arg("-c").arg(&colors)
        .output().unwrap();
    assert!(!result.status.success(), "expected failure but command succeeded");

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("have no matching read"),
            "expected error about unmatched pseudoalignment records, got: {}", stderr);
}

#[test]
fn import_sbwt_builds_single_colored_index() {
    use sbwt::{BitPackedKmerSortingMem, SbwtIndexBuilder, SbwtIndexVariant, write_sbwt_index_variant};

    let dir = tmp_dir();
    let sbwt_path = dir.join("input.sbwt");
    let index_path = dir.join("single.thm2");
    let query_file = "example/C1.fna";
    let k = 3;
    let color_name = "my_color";

    let (sbwt, _lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
        .k(k)
        .n_threads(1)
        .add_rev_comp(true)
        .build_lcs(false)
        .algorithm(BitPackedKmerSortingMem::default())
        .run_from_fasta(std::fs::File::open(PathBuf::from(PROJECT_DIR).join(query_file)).unwrap());

    let mut sbwt_out = std::io::BufWriter::new(std::fs::File::create(&sbwt_path).unwrap());
    write_sbwt_index_variant(&SbwtIndexVariant::SubsetMatrix(sbwt), &mut sbwt_out).unwrap();
    drop(sbwt_out);

    let status = themisto2()
        .args(["import-sbwt", "-s"]).arg(&sbwt_path)
        .args(["-o"]).arg(&index_path)
        .args(["--color-name", color_name, "-t", "1"])
        .status()
        .unwrap();
    assert!(status.success(), "import-sbwt failed");

    let dump = themisto2()
        .args(["dump-color-names", "-i"]).arg(&index_path)
        .output()
        .unwrap();
    assert!(dump.status.success(), "dump-color-names failed: {}", String::from_utf8_lossy(&dump.stderr));
    assert_eq!(String::from_utf8(dump.stdout).unwrap(), format!("0\t{}\n", color_name));

    let pa_out = dir.join("pa.jsonl");
    let status = themisto2()
        .args(["intersection-pseudoalign", "-i"]).arg(&index_path)
        .args(["-q", query_file, "-t", "1", "--sort-output", "-o"]).arg(&pa_out)
        .status()
        .unwrap();
    assert!(status.success(), "intersection-pseudoalign on imported single-colored index failed");

    let records = read_jsonl(&pa_out);
    assert!(!records.is_empty(), "no pseudoalignment records produced");
    for rec in &records {
        assert_eq!(rec.colors, vec![0], "expected single color 0 for every read, got {:?}", rec.colors);
    }
}

/// Run `cmd`, piping the contents of `query_file` to its stdin, and return the
/// captured stdout bytes. Asserts the process exited successfully.
fn run_with_stdin(cmd: &mut Command, query_file: &str) -> Vec<u8> {
    let query_bytes = std::fs::read(PathBuf::from(PROJECT_DIR).join(query_file)).unwrap();
    let mut child = cmd
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .unwrap();
    child.stdin.as_mut().unwrap().write_all(&query_bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "stdin-fed run failed for {}", query_file);
    output.stdout
}
