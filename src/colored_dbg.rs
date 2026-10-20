//! A colored de Bruijn graph view on top of a Themisto 2 index.
//!
//! The graph is the **doubled** node-centric de Bruijn graph: every k-mer in the index is a node,
//! and a k-mer and its reverse complement are separate nodes. There is an edge from x to y iff
//! x[1..k) = y[0..k-1), labeled with the last character of y. Reverse complements are not modeled
//! by the edges. Since the index contains both strands, the graph contains both orientations of
//! every unitig.
//!
//! Nodes are identified by their colex rank in the SBWT ([Node::id]).

use std::ops::Range;

use sbwt::dbg::Dbg;
use sbwt::{SubsetMatrix, SubsetSeq};
use simple_sds_sbwt::ops::BitVec;

use crate::colex_colored_kmers::CompactColexKmers;
use crate::sparse_dense_storage::{SparseDenseColorSetView, SparseDenseStorage};

pub use sbwt::dbg::Node;

/// A color set of a node. Iterate the colors with [crate::ColorSetView::iter].
pub type ColorSet<'a> = SparseDenseColorSetView<'a>;

/// A maximal run of consecutive k-mers in a unitig that share the same color set.
#[derive(Clone, Debug)]
pub struct ColorRun<'a> {
    /// Range of k-mer positions in the unitig (indices into [Unitig::nodes]).
    pub kmers: Range<usize>,
    pub color_set_id: usize,
    pub color_set: ColorSet<'a>,
}

/// A maximal unitig of the de Bruijn graph, with its color runs.
#[derive(Clone, Debug)]
pub struct Unitig<'u, 'a> {
    /// Nodes of the unitig in path order.
    pub nodes: &'u [Node],
    /// The spelled sequence: length nodes.len() + k - 1.
    pub seq: &'u [u8],
    /// Color runs covering [Unitig::nodes] from left to right.
    pub color_runs: &'u [ColorRun<'a>],
}

/// A maximal path of a unitig where all k-mers have the same color set.
#[derive(Clone, Debug)]
pub struct ColoredUnitig<'u, 'a> {
    /// Nodes of the colored unitig in path order.
    pub nodes: &'u [Node],
    /// The spelled sequence: length nodes.len() + k - 1.
    pub seq: &'u [u8],
    pub color_set_id: usize,
    pub color_set: ColorSet<'a>,
}

/// Colored de Bruijn graph operations on a borrowed Themisto 2 index.
pub struct ColoredDbg<'a> {
    index: &'a CompactColexKmers<SparseDenseStorage>,
    dbg: Dbg<'a, SubsetMatrix>,
}

impl<'a> ColoredDbg<'a> {

    /// Initializes the graph structure. Takes 2 bits of space per SBWT set on top of the index.
    /// Panics if the index was loaded without select support
    /// (see [CompactColexKmers::load]).
    pub fn new(index: &'a CompactColexKmers<SparseDenseStorage>, n_threads: usize) -> Self {
        assert!(index.sbwt().sbwt().has_select_support(), "ColoredDbg requires an index loaded with select support");
        let dbg = Dbg::new(index.sbwt(), Some(index.lcs()), n_threads);
        Self{index, dbg}
    }

    pub fn index(&self) -> &'a CompactColexKmers<SparseDenseStorage> {
        self.index
    }

    pub fn k(&self) -> usize {
        self.index.get_k()
    }

    pub fn n_nodes(&self) -> usize {
        self.index.sbwt().n_kmers()
    }

    pub fn n_colors(&self) -> usize {
        self.index.get_color_names().len()
    }

    pub fn color_names(&self) -> &'a [String] {
        self.index.get_color_names()
    }

    /// Iterates all nodes in colex order.
    pub fn nodes(&self) -> impl Iterator<Item = Node> + '_ {
        (0..self.index.sbwt().n_sets())
            .filter(|&colex| !self.dbg.is_dummy_colex_position(colex))
            .map(|id| Node{id})
    }

    /// Returns the node of the given k-mer, if it exists. Panics if `kmer.len() != k`.
    pub fn get_node(&self, kmer: &[u8]) -> Option<Node> {
        self.dbg.get_node(kmer)
    }

    pub fn get_kmer(&self, node: Node) -> Vec<u8> {
        self.dbg.get_kmer(node)
    }

    /// Appends the k-mer of `node` to `buf`.
    pub fn push_kmer(&self, node: Node, buf: &mut Vec<u8>) {
        self.dbg.push_node_kmer(node, buf)
    }

    pub fn outdegree(&self, node: Node) -> usize {
        self.dbg.outdegree(node)
    }

    pub fn indegree(&self, node: Node) -> usize {
        self.dbg.indegree(node)
    }

    /// Appends (target, edge label) pairs of the outgoing edges of `node` to `out`.
    pub fn push_out_neighbors(&self, node: Node, out: &mut Vec<(Node, u8)>) {
        self.dbg.push_out_neighbors(node, out)
    }

    /// Appends (source, edge label) pairs of the incoming edges of `node` to `out`.
    /// The label is the same for all of them: the last character of `node`.
    pub fn push_in_neighbors(&self, node: Node, out: &mut Vec<(Node, u8)>) {
        self.dbg.push_in_neighbors(node, out)
    }

    /// Follows the outgoing edge with the given label, if it exists.
    pub fn follow_outedge(&self, node: Node, label: u8) -> Option<Node> {
        self.dbg.follow_outedge(node, label)
    }

    /// Takes time proportional to the distance to the next sampled k-mer.
    /// For many nodes along a path, prefer the unitig iterators.
    pub fn color_set_id(&self, node: Node) -> usize {
        self.index.colex_to_set_id(node.id)
    }

    pub fn color_set(&self, node: Node) -> ColorSet<'a> {
        self.index.colex_to_set(node.id)
    }

    pub fn color_set_by_id(&self, color_set_id: usize) -> ColorSet<'a> {
        self.index.set_id_to_set(color_set_id)
    }

    /// Calls `callback` on every maximal unitig in parallel, in no particular order.
    /// The starting node of cyclic unitigs is arbitrary.
    pub fn par_for_each_unitig<F: Fn(&Unitig<'_, 'a>) + Send + Sync>(&self, callback: F, n_threads: usize) {
        self.dbg.iter_unitigs_with_callback(|nodes| {
            let mut seq = Vec::<u8>::with_capacity(nodes.len() + self.k() - 1);
            self.dbg.push_unitig_string(nodes, &mut seq);
            let color_runs = self.color_runs(nodes);
            callback(&Unitig{nodes, seq: &seq, color_runs: &color_runs});
        }, n_threads);
    }

    /// Calls `callback` on every maximal unitig piece with a single color set, in parallel,
    /// in no particular order. These are the unitigs from [ColoredDbg::par_for_each_unitig]
    /// split at color set changes.
    pub fn par_for_each_colored_unitig<F: Fn(&ColoredUnitig<'_, 'a>) + Send + Sync>(&self, callback: F, n_threads: usize) {
        let k = self.k();
        self.par_for_each_unitig(|unitig| {
            for run in unitig.color_runs {
                callback(&ColoredUnitig{
                    nodes: &unitig.nodes[run.kmers.clone()],
                    seq: &unitig.seq[run.kmers.start .. run.kmers.end + k - 1],
                    color_set_id: run.color_set_id,
                    color_set: run.color_set,
                });
            }
        }, n_threads);
    }

    // The color set of a node is the color set of the next sampled node forward in the graph,
    // so we scan the path backward. If the last node is not sampled (can happen on a cycle),
    // we look its set up by walking forward around the cycle.
    fn color_runs(&self, nodes: &[Node]) -> Vec<ColorRun<'a>> {
        let sampling = &self.index.get_map().sampling;
        let mut runs = Vec::<ColorRun>::new();
        if nodes.is_empty() { return runs }

        let mut cur_id = self.color_set_id(*nodes.last().unwrap());
        let mut run_end = nodes.len();
        for (pos, v) in nodes.iter().enumerate().rev() {
            if sampling.get(v.id) {
                let id = self.color_set_id(*v);
                if id != cur_id {
                    runs.push(ColorRun{kmers: pos+1..run_end, color_set_id: cur_id, color_set: self.color_set_by_id(cur_id)});
                    run_end = pos + 1;
                    cur_id = id;
                }
            }
        }
        runs.push(ColorRun{kmers: 0..run_end, color_set_id: cur_id, color_set: self.color_set_by_id(cur_id)});
        runs.reverse();
        runs
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Mutex;

    use super::*;
    use crate::coloring_interface::ColorSetView;

    fn colors(set: &ColorSet) -> Vec<usize> {
        set.iter().collect()
    }

    fn rc(seq: &[u8]) -> Vec<u8> {
        jseqio::reverse_complement(seq)
    }

    // Brute-force k-mer -> colors from the input, both strands
    fn true_colors(seqs: &[(&[u8], usize)], k: usize) -> HashMap<Vec<u8>, Vec<usize>> {
        let mut map = HashMap::<Vec<u8>, Vec<usize>>::new();
        for (seq, color) in seqs {
            for s in [seq.to_vec(), rc(seq)] {
                for kmer in s.windows(k) {
                    let e = map.entry(kmer.to_vec()).or_default();
                    if !e.contains(color) { e.push(*color) }
                }
            }
        }
        for v in map.values_mut() { v.sort() }
        map
    }

    fn test_input() -> Vec<(&'static [u8], usize)> {
        vec![
            (b"ACGTTGCATGCAAGTCCGATAGCTAGGCTA".as_slice(), 0),
            (b"TGCATGCAAGTCCGATTTTTTACG".as_slice(), 1),
            (b"GGGGGCCCCCAAAAATTTTT".as_slice(), 2),
            (b"CAGTCAGTCAGTCAGT".as_slice(), 0), // cyclic in the dBG
        ]
    }

    fn build(seqs: &[(&[u8], usize)], k: usize) -> CompactColexKmers<SparseDenseStorage> {
        let both_strands: Vec<(Vec<u8>, usize)> = seqs.iter().flat_map(|(s, c)| [(s.to_vec(), *c), (rc(s), *c)]).collect();
        let refs: Vec<(&[u8], usize)> = both_strands.iter().map(|(s, c)| (s.as_slice(), *c)).collect();
        CompactColexKmers::new_from_small_input(&refs, k, 3, 2)
    }

    #[test]
    fn node_ops_and_colors() {
        let k = 5;
        let seqs = test_input();
        let index = build(&seqs, k);
        let g = ColoredDbg::new(&index, 2);
        let truth = true_colors(&seqs, k);

        assert_eq!(g.n_nodes(), truth.len());
        assert_eq!(g.nodes().count(), truth.len());
        for v in g.nodes() {
            let kmer = g.get_kmer(v);
            assert_eq!(g.get_node(&kmer), Some(v));
            assert_eq!(colors(&g.color_set(v)), truth[&kmer]);

            let mut out = vec![];
            g.push_out_neighbors(v, &mut out);
            assert_eq!(out.len(), g.outdegree(v));
            for (u, c) in out {
                assert_eq!(g.get_kmer(u)[..k-1], kmer[1..]);
                assert_eq!(*g.get_kmer(u).last().unwrap(), c);
                assert_eq!(g.follow_outedge(v, c), Some(u));
            }

            let mut inn = vec![];
            g.push_in_neighbors(v, &mut inn);
            assert_eq!(inn.len(), g.indegree(v));
            for (u, _) in inn {
                assert_eq!(g.get_kmer(u)[1..], kmer[..k-1]);
            }
        }
        assert_eq!(g.get_node(b"AAAAC"), None);
    }

    #[test]
    fn unitigs_cover_graph_with_correct_colors() {
        let k = 5;
        let seqs = test_input();
        let index = build(&seqs, k);
        let g = ColoredDbg::new(&index, 2);
        let truth = true_colors(&seqs, k);

        let seen = Mutex::new(Vec::<Vec<u8>>::new());
        g.par_for_each_unitig(|u| {
            assert_eq!(u.seq.len(), u.nodes.len() + k - 1);
            assert_eq!(u.color_runs.first().unwrap().kmers.start, 0);
            assert_eq!(u.color_runs.last().unwrap().kmers.end, u.nodes.len());
            for w in u.color_runs.windows(2) {
                assert_eq!(w[0].kmers.end, w[1].kmers.start);
                assert_ne!(w[0].color_set_id, w[1].color_set_id);
            }
            for run in u.color_runs {
                for i in run.kmers.clone() {
                    let kmer = &u.seq[i..i+k];
                    assert_eq!(g.get_kmer(u.nodes[i]), kmer);
                    assert_eq!(colors(&run.color_set), truth[kmer]);
                }
            }
            seen.lock().unwrap().extend(u.seq.windows(k).map(|x| x.to_vec()));
        }, 2);
        let mut seen = seen.into_inner().unwrap();
        seen.sort();
        let mut expected: Vec<Vec<u8>> = truth.keys().cloned().collect();
        expected.sort();
        assert_eq!(seen, expected); // Every k-mer exactly once

        let seen = Mutex::new(Vec::<Vec<u8>>::new());
        g.par_for_each_colored_unitig(|u| {
            assert_eq!(u.seq.len(), u.nodes.len() + k - 1);
            for kmer in u.seq.windows(k) {
                assert_eq!(colors(&u.color_set), truth[kmer]);
            }
            seen.lock().unwrap().extend(u.seq.windows(k).map(|x| x.to_vec()));
        }, 2);
        let mut seen = seen.into_inner().unwrap();
        seen.sort();
        assert_eq!(seen, expected);
    }

    #[test]
    fn cyclic_unitig_with_color_change() {
        // AACCG repeated is a cycle of 5 4-mers. Two of them also have color 1.
        let k = 4;
        let seqs: Vec<(&[u8], usize)> = vec![(b"AACCGAACCGAACCG", 0), (b"AACCG", 1)];
        let index = build(&seqs, k);
        let g = ColoredDbg::new(&index, 1);
        let truth = true_colors(&seqs, k);

        let n_unitigs = Mutex::new(0);
        g.par_for_each_unitig(|u| {
            *n_unitigs.lock().unwrap() += 1;
            assert_eq!(u.nodes.len(), 5);
            assert!(u.color_runs.len() >= 2);
            // Every rotation of the cycle, so that the last node is sometimes not sampled
            for r in 0..u.nodes.len() {
                let mut rotated = u.nodes.to_vec();
                rotated.rotate_left(r);
                for run in g.color_runs(&rotated) {
                    for i in run.kmers.clone() {
                        assert_eq!(colors(&run.color_set), truth[&g.get_kmer(rotated[i])]);
                    }
                }
            }
        }, 1);
        assert_eq!(n_unitigs.into_inner().unwrap(), 2); // Both strands
    }
}
