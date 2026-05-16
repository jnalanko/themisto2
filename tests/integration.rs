use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

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
        .args(["-k", "3", "--index-type", "sparse-dense", "-t", "1", "-o"])
        .arg(&index1)
        .status()
        .unwrap();
    assert!(status.success(), "--file-colors build failed");

    let status = themisto2()
        .args(["build", "--seq-colors", "example/seq-colors.fna", "--temp-dir"])
        .arg(&tmp2)
        .args(["-k", "3", "--index-type", "sparse-dense", "-t", "1", "-o"])
        .arg(&index2)
        .status()
        .unwrap();
    assert!(status.success(), "--seq-colors build failed");

    let bytes1 = std::fs::read(&index1).unwrap();
    let bytes2 = std::fs::read(&index2).unwrap();
    assert_eq!(bytes1, bytes2, "Indexes built with --file-colors and --seq-colors differ");
}
