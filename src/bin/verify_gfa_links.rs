use std::fs::File;
use std::io::{BufRead, BufReader};
use rustc_hash::{FxHashMap, FxHashSet};
use jseqio::reverse_complement;

fn parse_seq_id(name: &str) -> usize {
    name.strip_prefix('u')
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| panic!("unexpected segment name {:?} (expected u{{id}})", name))
}

fn parse_sign(s: &str) -> bool {
    match s {
        "+" => true,
        "-" => false,
        _ => panic!("Error parsing sign: {s}")
    }
}

#[derive(Eq, PartialEq, Hash, Debug)]
struct Edge {
    from: usize,
    to: usize,
    from_sign: bool,
    to_sign: bool,
}

fn main() {
    let cli = clap::Command::new("verify_gfa_links")
        .about("Verify that every link in a Themisto 2 GFA export is correct and no link is missing")
        .arg(clap::Arg::new("gfa")
            .required(true)
            .help("GFA file to verify"))
        .arg(clap::Arg::new("k")
            .short('k')
            .long("kmer-length")
            .required(true)
            .value_parser(clap::value_parser!(usize))
            .help("k-mer length used during index construction"));

    let args = cli.get_matches();
    let gfa_path = args.get_one::<String>("gfa").unwrap();
    let k = *args.get_one::<usize>("k").unwrap();
    assert!(k >= 2, "k must be at least 2");

    eprintln!("Parsing GFA (k={})...", k);

    // Id -> sequence
    let mut seqs: FxHashMap<usize, Vec<u8>> = FxHashMap::default();

    let mut parsed_edges: FxHashSet<Edge> = FxHashSet::default();

    // Maps: (k-1)-mer-string -> unitig ids
    let mut start_kmers = FxHashMap::<Vec<u8>, Vec<usize>>::default();
    let mut start_rc_kmers = FxHashMap::<Vec<u8>, Vec<usize>>::default();

    let reader = BufReader::new(
        File::open(gfa_path).unwrap_or_else(|e| panic!("cannot open {}: {}", gfa_path, e)),
    );
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        if line.is_empty() { continue; }
        let f: Vec<&str> = line.splitn(7, '\t').collect();
        match f[0] {
            "H" => {}
            "S" => {
                assert!(f.len() >= 3, "line {}: S-line needs at least 3 fields", lineno + 1);
                let id = parse_seq_id(f[1]);
                let seq = f[2].as_bytes().to_vec();
                assert!(seq.len() >= k,
                    "line {}: segment u{} has length {} < k={}", lineno + 1, id, seq.len(), k);
                seqs.insert(id, seq.clone());

                let first = &seq[0..k-1];
                if !start_kmers.contains_key(first) {
                    start_kmers.insert(first.to_owned(), vec![]);
                }
                start_kmers.get_mut(first).unwrap().push(id);

                let first_rc = reverse_complement(&seq[seq.len()-(k-1)..]);
                if !start_rc_kmers.contains_key(&first_rc) {
                    start_rc_kmers.insert(first_rc.to_owned(), vec![]);
                }

                start_rc_kmers.get_mut(&first_rc).unwrap().push(id);
            }
            "L" => {
                assert!(f.len() >= 6, "line {}: L-line needs at least 6 fields", lineno + 1);
                let a  = parse_seq_id(f[1]);
                let sa = parse_sign(f[2]);
                let b  = parse_seq_id(f[3]);
                let sb = parse_sign(f[4]);
                let ov: usize = f[5].trim_end_matches('M').parse()
                    .unwrap_or_else(|_| panic!("line {}: bad overlap field {:?}", lineno + 1, f[5]));
                assert_eq!(ov, k - 1,
                    "line {}: overlap {} != k-1={}", lineno + 1, ov, k - 1);
                let e = Edge {
                    from: a,
                    to: b,
                    from_sign: sa,
                    to_sign: sb,
                };
                parsed_edges.insert(e);
            }
            _ => {}
        }
    }

    eprintln!("Parsed {} unitigs and {} edges", seqs.len(), parsed_edges.len());

    // Build the true edges
    let mut true_edges: FxHashSet<Edge> = FxHashSet::default();
    for (&from_id, from_string) in seqs.iter() {

        // ++ edges (and the flipped twin -- edges)
        let last = &from_string[from_string.len()-(k-1)..];
        if let Some(list) = start_kmers.get(last) { for &to_id in list {
            let e = Edge { from: from_id, to: to_id, from_sign: true, to_sign: true };
            true_edges.insert(e);
            let e2 = Edge { from: to_id, to: from_id, from_sign: false, to_sign: false };
            true_edges.insert(e2);
        }};

        // +- edges (and the flipped +- twin edges)
        if let Some(list) = start_rc_kmers.get(last) { for &to_id in list {
            let e = Edge { from: from_id, to: to_id, from_sign: true, to_sign: false };
            true_edges.insert(e);
            let e2 = Edge { from: to_id, to: from_id, from_sign: true, to_sign: false };
            true_edges.insert(e2);
        }};

        // -+ edges (and the flipped -+ twin edges)
        let from_string_rc = reverse_complement(from_string);
        let rc_last = &from_string_rc[from_string_rc.len()-(k-1)..];
        if let Some(list) = start_kmers.get(rc_last) { for &to_id in list {
            let e = Edge { from: from_id, to: to_id, from_sign: false, to_sign: true };
            true_edges.insert(e);
            let e2 = Edge { from: to_id, to: from_id, from_sign: false, to_sign: true };
            true_edges.insert(e2);
        }};
    }

    eprintln!("Computed {} true edges", true_edges.len());

    for e in true_edges.iter() {
        assert!(parsed_edges.contains(e), "Edge {e:?} missing from GFA");
    }

    for e in parsed_edges.iter() {
        assert!(true_edges.contains(e), "Edge {e:?} should not be in GFA");
    }

    eprintln!("Everything ok");
}
