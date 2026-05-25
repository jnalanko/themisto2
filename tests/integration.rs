use std::path::PathBuf;
use std::process::Command;
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
    }
}
