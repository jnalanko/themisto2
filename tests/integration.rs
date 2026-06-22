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
    // Pipe stderr so subprocess log output doesn't leak into the test terminal.
    // Tests that need to inspect stderr can still read it via .output().stderr.
    cmd.stderr(Stdio::piped());
    cmd
}


#[test]
fn file_colors_and_seq_colors_single_seq_per_color_produce_same_index() {
    let dir = tmp_dir();
    let tmp1 = dir.join("tmp1");
    let tmp2 = dir.join("tmp2");
    std::fs::create_dir_all(&tmp1).unwrap();
    std::fs::create_dir_all(&tmp2).unwrap();

    let index1 = dir.join("file_colors.thm2");
    let index2 = dir.join("seq_colors.thm2");

    let status = themisto2()
        .args(["build", "--file-colors", "tests/data/fof.txt", "--temp-dir"])
        .arg(&tmp1)
        .args(["-k", "3", "-t", "1", "-o"])
        .arg(&index1)
        .status()
        .unwrap();
    assert!(status.success(), "--file-colors build failed");

    let status = themisto2()
        .args(["build", "--seq-colors", "tests/data/seq-colors.fna", "--temp-dir"])
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

fn build_sbwt_file(fasta_path: &PathBuf, sbwt_path: &PathBuf, k: usize) {
    use sbwt::{BitPackedKmerSortingMem, SbwtIndexBuilder, SbwtIndexVariant, write_sbwt_index_variant};
    let (sbwt, _lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
        .k(k)
        .n_threads(1)
        .add_rev_comp(true)
        .build_lcs(false)
        .algorithm(BitPackedKmerSortingMem::default())
        .run_from_fasta(std::fs::File::open(fasta_path).unwrap());
    let mut out = std::io::BufWriter::new(std::fs::File::create(sbwt_path).unwrap());
    write_sbwt_index_variant(&SbwtIndexVariant::SubsetMatrix(sbwt), &mut out).unwrap();
}

#[test]
fn import_sbwt_builds_single_colored_index() {
    let dir = tmp_dir();
    let sbwt_path = dir.join("input.sbwt");
    let index_path = dir.join("single.thm2");
    let query_file = "example/C1.fna";
    let color_name = "my_color";

    build_sbwt_file(&PathBuf::from(PROJECT_DIR).join(query_file), &sbwt_path, 3);

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

#[test]
fn merged_single_colored_indexes_match_normal_build() {
    let dir = tmp_dir();
    let k = 3;

    // Step 1: build a single-colored Themisto 2 index for each input file
    // by first dumping an SBWT and then running `import-sbwt`.
    let mut single_colored_indexes: Vec<PathBuf> = Vec::new();
    for fname in ["C1.fna", "C2.fna"] {
        let fasta = PathBuf::from(PROJECT_DIR).join("example").join(fname);
        let sbwt = dir.join(format!("{}.sbwt", fname));
        let index = dir.join(format!("{}.thm2", fname));
        build_sbwt_file(&fasta, &sbwt, k);

        // The normal build uses the input path string as the color name, so
        // use the same convention here to keep the color names aligned with
        // the normal build for byte-by-byte output comparison below.
        let color_name = format!("example/{}", fname);

        let status = themisto2()
            .args(["import-sbwt", "-s"]).arg(&sbwt)
            .args(["-o"]).arg(&index)
            .args(["--color-name", &color_name, "-t", "1"])
            .status()
            .unwrap();
        assert!(status.success(), "import-sbwt failed for {}", fname);
        single_colored_indexes.push(index);
    }

    // Step 2: merge the two single-colored indexes.
    let index_list = dir.join("indexes.txt");
    {
        let mut s = String::new();
        for p in &single_colored_indexes {
            s.push_str(p.to_str().unwrap());
            s.push('\n');
        }
        write_file(&index_list, &s);
    }
    let merge_temp = dir.join("merge_temp");
    std::fs::create_dir_all(&merge_temp).unwrap();
    let merged_index = dir.join("merged.thm2");
    let status = themisto2()
        .args(["merge", "--index-file-list"]).arg(&index_list)
        .args(["--temp-dir"]).arg(&merge_temp)
        .args(["-o"]).arg(&merged_index)
        .args(["-t", "1"])
        .status()
        .unwrap();
    assert!(status.success(), "merge failed");

    // Step 3: normal two-color build for comparison.
    let fof = dir.join("fof_two.txt");
    write_file(&fof, "example/C1.fna\nexample/C2.fna\n");
    let normal_temp = dir.join("normal_temp");
    std::fs::create_dir_all(&normal_temp).unwrap();
    let normal_index = dir.join("normal.thm2");
    let status = themisto2()
        .args(["build", "--file-colors"]).arg(&fof)
        .args(["--temp-dir"]).arg(&normal_temp)
        .args(["-k", &k.to_string(), "-t", "1", "-o"]).arg(&normal_index)
        .status()
        .unwrap();
    assert!(status.success(), "normal build failed");

    // Step 4: color names should match between the two indexes.
    let merged_names = themisto2()
        .args(["dump-color-names", "-i"]).arg(&merged_index)
        .output().unwrap();
    assert!(merged_names.status.success());
    let normal_names = themisto2()
        .args(["dump-color-names", "-i"]).arg(&normal_index)
        .output().unwrap();
    assert!(normal_names.status.success());
    assert_eq!(merged_names.stdout, normal_names.stdout, "color names differ between merged and normal indexes");

    // Step 5: pseudoalignment results must match for every input file
    // (covers both single-color and shared-color k-mers).
    for query_file in ["example/C1.fna", "example/C2.fna"] {
        let merged_out = dir.join(format!("merged-{}.jsonl", query_file.replace('/', "_")));
        let normal_out = dir.join(format!("normal-{}.jsonl", query_file.replace('/', "_")));

        for (idx, out) in [(&merged_index, &merged_out), (&normal_index, &normal_out)] {
            let status = themisto2()
                .args(["intersection-pseudoalign", "-i"]).arg(idx)
                .args(["-q", query_file, "-t", "1", "--sort-output", "-o"]).arg(out)
                .status().unwrap();
            assert!(status.success(), "intersection-pseudoalign failed for {} against {}", query_file, idx.display());
        }

        assert_eq!(
            std::fs::read(&merged_out).unwrap(),
            std::fs::read(&normal_out).unwrap(),
            "pseudoalignment of {} differs between merged and normal indexes", query_file,
        );

        // Step 6: per-k-mer color sets must match too.
        let merged_sets = dir.join(format!("merged-{}.colorsets.txt", query_file.replace('/', "_")));
        let normal_sets = dir.join(format!("normal-{}.colorsets.txt", query_file.replace('/', "_")));
        for (idx, out) in [(&merged_index, &merged_sets), (&normal_index, &normal_sets)] {
            let status = themisto2()
                .args(["print-color-sets", "-i"]).arg(idx)
                .args(["-q", query_file, "-p", "-o"]).arg(out)
                .status().unwrap();
            assert!(status.success(), "print-color-sets failed for {} against {}", query_file, idx.display());
        }

        assert_eq!(
            std::fs::read(&merged_sets).unwrap(),
            std::fs::read(&normal_sets).unwrap(),
            "per-k-mer color sets of {} differ between merged and normal indexes", query_file,
        );
    }
}

#[test]
fn report_smoke_test() {
    let dir = tmp_dir();
    let tmp = dir.join("tmp");
    std::fs::create_dir_all(&tmp).unwrap();

    // Build a small index with 3 colors.
    let index = dir.join("index.thm2");
    let status = themisto2()
        .args(["build", "--file-colors", "example/fof.txt", "--temp-dir"])
        .arg(&tmp)
        .args(["-k", "3", "-t", "1", "-o"])
        .arg(&index)
        .status()
        .unwrap();
    assert!(status.success(), "build failed");

    // Pseudoalign all three input files together to get a JSONL file.
    let pa_out = dir.join("pa.jsonl");
    let status = themisto2()
        .args(["intersection-pseudoalign", "-i"])
        .arg(&index)
        .args(["-q", "example/C1.fna", "-t", "1", "-o"])
        .arg(&pa_out)
        .status()
        .unwrap();
    assert!(status.success(), "pseudoalign failed");

    // Run report with an explicit -p and -o.
    let report_out = dir.join("report.json");
    let status = themisto2()
        .args(["report", "-i"])
        .arg(&index)
        .arg("-p").arg(&pa_out)
        .arg("-o").arg(&report_out)
        .status()
        .unwrap();
    assert!(status.success(), "report failed");

    let text = std::fs::read_to_string(&report_out).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("report output is not valid JSON: {e}\n---\n{text}"));

    // Structural checks.
    assert!(v["n_reads"].is_number(), "missing n_reads");
    assert!(v["n_positive_reads"].is_number(), "missing n_positive_reads");
    assert!(v["positive_by_color"].is_array(), "missing positive_by_color");
    assert!(v["unique_positive_by_color"].is_array(), "missing unique_positive_by_color");

    let n_reads = v["n_reads"].as_u64().unwrap();
    let n_positive = v["n_positive_reads"].as_u64().unwrap();
    assert!(n_reads > 0, "expected at least one read in report");
    assert!(n_positive <= n_reads, "n_positive_reads ({}) > n_reads ({})", n_positive, n_reads);

    // Without -p and -o the same output should arrive on stdout.
    let pa_bytes = std::fs::read(&pa_out).unwrap();
    let mut child = themisto2()
        .args(["report", "-i"]).arg(&index)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.as_mut().unwrap().write_all(&pa_bytes).unwrap();
    let stdout_run = child.wait_with_output().unwrap();
    assert!(stdout_run.status.success(), "report (stdin→stdout) failed");
    assert_eq!(
        stdout_run.stdout,
        std::fs::read(&report_out).unwrap(),
        "report output via stdin/stdout differs from -p/-o output",
    );
}

/// Run `cmd`, piping the contents of `query_file` to its stdin, and return the
/// captured stdout bytes. Asserts the process exited successfully.
fn run_with_stdin(cmd: &mut Command, query_file: &str) -> Vec<u8> {
    let query_bytes = std::fs::read(PathBuf::from(PROJECT_DIR).join(query_file)).unwrap();
    let mut child = cmd
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped()) // capture stderr so it doesn't leak into test output
        .spawn()
        .unwrap();
    child.stdin.as_mut().unwrap().write_all(&query_bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "stdin-fed run failed for {}", query_file);
    output.stdout
}
