// Counts the number of distinct color sets of k'-mers for k' = k_min..=k, where
// the color set of a k'-mer is the union of the color sets of the k-mers that
// have that k'-mer as a suffix.
//
// The k-mers sharing a k'-suffix form a contiguous colex interval: a maximal run
// where lcs[i] >= k'. These intervals nest like LCP intervals of a suffix array,
// so we compute the unions for all k' with a single bottom-up stack traversal
// of the LCS array. An interval at depth d whose parent interval is at depth p
// is the interval of a k'-mer for every k' in (p, d], with the same union.
//
// k'-mers whose union is empty (only dummy k-mers have it as a suffix, or the
// k'-mer contains $) are not counted.

use std::hash::BuildHasher;
use std::collections::hash_map::RandomState;
use std::ops::Range;

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rustc_hash::FxHashSet;

use crate::colex_colored_kmers::CompactColexKmers;
use crate::coloring_interface::{ColorSetStorage, ColorSetView};

// 128-bit fingerprint of a color set from two independently keyed SipHashes.
struct Fingerprinter {
    h1: RandomState,
    h2: RandomState,
}

impl Fingerprinter {
    fn new() -> Self {
        Self { h1: RandomState::new(), h2: RandomState::new() }
    }

    fn fingerprint(&self, set: &[u32]) -> u128 {
        ((self.h1.hash_one(set) as u128) << 64) | self.h2.hash_one(set) as u128
    }
}

// Per-k' sets of fingerprints. fingerprints[i] is for k' = k_min + i.
type FingerprintSets = Vec<FxHashSet<u128>>;

// a := a ∪ b. Both must be sorted. Uses buf as scratch space.
fn union_into(a: &mut Vec<u32>, b: &[u32], buf: &mut Vec<u32>) {
    if b.is_empty() || a.as_slice() == b {
        return;
    }
    if a.is_empty() {
        a.extend_from_slice(b);
        return;
    }
    buf.clear();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => { buf.push(a[i]); i += 1; },
            std::cmp::Ordering::Greater => { buf.push(b[j]); j += 1; },
            std::cmp::Ordering::Equal => { buf.push(a[i]); i += 1; j += 1; },
        }
    }
    buf.extend_from_slice(&a[i..]);
    buf.extend_from_slice(&b[j..]);
    std::mem::swap(a, buf);
}

struct ChunkProcessor<'a, CSS: ColorSetStorage> {
    index: &'a CompactColexKmers<CSS>,
    dummy_marks: &'a bitvec::vec::BitVec,
    fingerprinter: &'a Fingerprinter,
    k_min: usize,
}

impl<CSS: ColorSetStorage> ChunkProcessor<'_, CSS> {

    // Records the union of an interval at depth d with parent depth p, i.e.
    // the union of the k'-mers for k' in (p, d].
    fn report(&self, union: &[u32], p: usize, d: usize, out: &mut FingerprintSets) {
        let first = std::cmp::max(p + 1, self.k_min);
        if union.is_empty() || first > d {
            return;
        }
        let fp = self.fingerprinter.fingerprint(union);
        for k_prime in first..=d {
            out[k_prime - self.k_min].insert(fp);
        }
    }

    // Processes the colex range. The caller must guarantee that no interval of a
    // k'-mer with k' >= k_min crosses the boundaries of the range, i.e. lcs[range.start] < k_min
    // and lcs[range.end] < k_min (where defined).
    fn process(&self, range: Range<usize>, out: &mut FingerprintSets) {
        let k = self.index.get_k();
        let lcs = self.index.lcs();

        // Stack of (depth, union). The bottom is the root at depth 0, which is never reported.
        let mut stack: Vec<(usize, Vec<u32>)> = vec![(0, vec![])];
        let mut buf = Vec::<u32>::new();
        let mut prev_leaf: Option<(usize, Vec<u32>)> = None; // (set id, set) of the previous non-dummy k-mer

        for i in range.start..=range.end {
            // Treat the boundaries of the range as LCS 0
            let l = if i == range.start || i == range.end { 0 } else { lcs.access(i) };

            // Close the intervals deeper than l
            while stack.last().unwrap().0 > l {
                let (d, u) = stack.pop().unwrap();
                let top = stack.last_mut().unwrap();
                let p = std::cmp::max(l, top.0);
                self.report(&u, p, d, out);
                if top.0 >= l {
                    if top.1.is_empty() {
                        top.1 = u;
                    } else {
                        union_into(&mut top.1, &u, &mut buf);
                    }
                } else {
                    // New internal node at depth l that starts with the popped child
                    stack.push((l, u));
                }
            }

            if i == range.end {
                break;
            }

            // Push the k-mer at colex i as a leaf at depth k
            let leaf = if self.dummy_marks[i] {
                vec![]
            } else {
                let id = self.index.colex_to_set_id(i);
                match &prev_leaf {
                    Some((prev_id, prev_set)) if *prev_id == id => prev_set.clone(),
                    _ => {
                        let set: Vec<u32> = self.index.set_id_to_set(id).iter().map(|c| c as u32).collect();
                        prev_leaf = Some((id, set.clone()));
                        set
                    }
                }
            };
            stack.push((k, leaf));
        }
    }
}

// Splits 0..n into roughly n_pieces ranges such that every split point i has lcs[i] < k_min.
fn split_range<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, k_min: usize, n_pieces: usize) -> Vec<Range<usize>> {
    let n = index.sbwt().n_sets();
    let lcs = index.lcs();
    let piece_len = n.div_ceil(n_pieces.max(1)).max(1);
    let mut ranges = vec![];
    let mut start = 0;
    while start < n {
        let mut end = std::cmp::min(start + piece_len, n);
        while end < n && lcs.access(end) >= k_min {
            end += 1;
        }
        ranges.push(start..end);
        start = end;
    }
    ranges
}

/// Returns pairs (k', number of distinct non-empty color sets of k'-mers) for k' = k_min..=k.
pub fn count_distinct_suffix_color_sets<CSS: ColorSetStorage + Sync>(index: &CompactColexKmers<CSS>, k_min: usize, n_threads: usize) -> Vec<(usize, usize)> {
    let k = index.get_k();
    assert!(k_min >= 1 && k_min <= k, "k_min must be in the range [1, k]");
    let n_k_values = k - k_min + 1;

    log::info!("Marking dummy nodes");
    let dummy_marks = index.sbwt().compute_dummy_node_marks();

    let fingerprinter = Fingerprinter::new();
    let processor = ChunkProcessor { index, dummy_marks: &dummy_marks, fingerprinter: &fingerprinter, k_min };

    let ranges = split_range(index, k_min, n_threads * 64);
    log::info!("Processing {} colex ranges", ranges.len());

    let bar = indicatif::ProgressBar::new(index.sbwt().n_sets() as u64);
    let empty = || vec![FxHashSet::<u128>::default(); n_k_values];
    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let fingerprints = pool.install(|| {
        ranges.into_par_iter().fold(empty, |mut acc, range| {
            let len = range.len() as u64;
            processor.process(range, &mut acc);
            bar.inc(len);
            acc
        }).reduce(empty, |mut a, mut b| {
            for (x, y) in a.iter_mut().zip(b.iter_mut()) {
                if x.len() < y.len() {
                    std::mem::swap(x, y);
                }
                x.extend(y.drain());
            }
            a
        })
    });
    bar.finish();

    fingerprints.iter().enumerate().map(|(i, s)| (k_min + i, s.len())).collect()
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeSet, HashMap, HashSet};

    use rand::{Rng, SeedableRng};

    use super::*;
    use crate::sparse_dense_storage::SparseDenseStorage;

    fn brute_force(index: &CompactColexKmers<SparseDenseStorage>, k_min: usize) -> Vec<(usize, usize)> {
        let k = index.get_k();
        let mut kmers = vec![];
        for colex in 0..index.sbwt().n_sets() {
            let kmer = index.sbwt().access_kmer(colex);
            if kmer.contains(&b'$') {
                continue;
            }
            let set: BTreeSet<usize> = index.colex_to_set(colex).iter().collect();
            kmers.push((kmer, set));
        }

        (k_min..=k).map(|k_prime| {
            let mut unions = HashMap::<Vec<u8>, BTreeSet<usize>>::new();
            for (kmer, set) in kmers.iter() {
                unions.entry(kmer[k - k_prime..].to_vec()).or_default().extend(set.iter().copied());
            }
            let distinct: HashSet<&BTreeSet<usize>> = unions.values().collect();
            (k_prime, distinct.len())
        }).collect()
    }

    fn random_seq(rng: &mut impl Rng, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
    }

    // Mutates random positions to get related sequences with partially shared k-mers
    fn mutate(rng: &mut impl Rng, seq: &[u8], n_mutations: usize) -> Vec<u8> {
        let mut seq = seq.to_vec();
        for _ in 0..n_mutations {
            let i = rng.gen_range(0..seq.len());
            seq[i] = b"ACGT"[rng.gen_range(0..4)];
        }
        seq
    }

    #[test]
    fn matches_brute_force() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1234);
        for &(k, n_colors, seq_len) in &[(5, 3, 50), (7, 6, 200), (11, 10, 300), (31, 8, 500)] {
            let base = random_seq(&mut rng, seq_len);
            let seqs: Vec<Vec<u8>> = (0..n_colors).map(|_| mutate(&mut rng, &base, seq_len / 20)).collect();
            let colored_seqs: Vec<(&[u8], usize)> = seqs.iter().enumerate().map(|(c, s)| (s.as_slice(), c)).collect();
            let index = CompactColexKmers::<SparseDenseStorage>::new_from_small_input(&colored_seqs, k, 3, 2);

            for k_min in [1, 2, k / 2, k] {
                let expected = brute_force(&index, k_min);
                for n_threads in [1, 3] {
                    let got = count_distinct_suffix_color_sets(&index, k_min, n_threads);
                    assert_eq!(got, expected, "k = {}, k_min = {}, n_threads = {}", k, k_min, n_threads);
                }
            }
        }
    }
}
