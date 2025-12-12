#!/usr/bin/env python3
"""
Node-centric de Bruijn graph visualization from a FASTA file.

- Nodes are k-mers (length k).
- Directed edges connect consecutive k-mers in the input sequence(s),
  i.e., overlaps of length (k-1).
- Node weight = number of times the k-mer occurs across all sequences.
- Edge label = (k-1)-mer overlap (suffix of source == prefix of target).

Usage:
  python debruijn_nodecentric.py input.fasta -k 21 --max-nodes 500 --layout spring --out graph.png
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from typing import Dict, Iterable, Iterator, List, Tuple

import networkx as nx
import matplotlib.pyplot as plt


def read_fasta_sequences(path: str) -> List[str]:
    """Minimal FASTA parser. Returns sequences as uppercase strings (no whitespace)."""
    seqs: List[str] = []
    buf: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf).upper())
                    buf = []
            else:
                buf.append(line)
        if buf:
            seqs.append("".join(buf).upper())
    return seqs


def iter_kmers(seq: str, k: int) -> Iterator[str]:
    """Yield k-mers from seq, skipping any k-mer containing non-ACGT characters."""
    valid = set("ACGT")
    n = len(seq)
    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        if all(c in valid for c in kmer):
            yield kmer


def build_nodecentric_debruijn(
    sequences: Iterable[str], k: int
) -> Tuple[nx.DiGraph, Counter]:
    """
    Build node-centric de Bruijn graph:
      nodes = k-mers
      edges = consecutive k-mer transitions (overlap k-1)
      node weight = k-mer count
      edge attributes:
        - overlap: (k-1)-mer label
        - weight: number of times this transition occurs
    """
    if k < 2:
        raise ValueError("k must be >= 2 (so k-1 overlaps exist).")

    node_counts: Counter = Counter()
    edge_counts: Dict[Tuple[str, str], int] = defaultdict(int)

    for seq in sequences:
        kmers = list(iter_kmers(seq, k))
        node_counts.update(kmers)

        for a, b in zip(kmers, kmers[1:]):
            # by construction of consecutive kmers, they overlap by k-1 if no invalid chars were skipped;
            # still, guard for safety
            if a[1:] == b[:-1]:
                edge_counts[(a, b)] += 1

    G = nx.DiGraph()

    # add nodes with weight
    for kmer, c in node_counts.items():
        G.add_node(kmer, weight=int(c))

    # add edges with overlap label and transition weight
    for (a, b), w in edge_counts.items():
        overlap = a[1:]  # == b[:-1], length k-1
        G.add_edge(a, b, overlap=overlap, weight=int(w))

    return G, node_counts


def prune_to_top_nodes(G: nx.DiGraph, max_nodes: int) -> nx.DiGraph:
    """Keep only the top max_nodes nodes by node weight (descending)."""
    if max_nodes <= 0 or G.number_of_nodes() <= max_nodes:
        return G

    nodes_sorted = sorted(G.nodes, key=lambda n: G.nodes[n].get("weight", 1), reverse=True)
    keep = set(nodes_sorted[:max_nodes])

    # include neighbors to avoid isolated arrows? keep strictly top nodes for determinism
    H = G.subgraph(keep).copy()
    return H


def layout_positions(G: nx.DiGraph, layout: str, seed: int = 42):
    if layout == "spring":
        return nx.spring_layout(G, seed=seed, k=None)
    if layout == "kamada_kawai":
        return nx.kamada_kawai_layout(G)
    if layout == "circular":
        return nx.circular_layout(G)
    raise ValueError(f"Unknown layout: {layout}")


def draw_graph(
    G: nx.DiGraph,
    out_path: str | None,
    show: bool,
    with_labels: bool,
    label_max_len: int,
    figsize: Tuple[int, int],
    layout: str,
):
    if G.number_of_nodes() == 0:
        raise SystemExit("Graph is empty (no valid k-mers found).")

    pos = layout_positions(G, layout=layout)

    weights = [G.nodes[n].get("weight", 1) for n in G.nodes]
    w_min, w_max = min(weights), max(weights)

    # scale node sizes (area) for visibility
    # avoid huge blow-ups for very large counts
    def scale(w: int) -> float:
        if w_max == w_min:
            return 800.0
        # sqrt scaling tends to look good
        import math
        return 200.0 + 1800.0 * (math.sqrt(w) - math.sqrt(w_min)) / (math.sqrt(w_max) - math.sqrt(w_min))

    node_sizes = [scale(w) for w in weights]

    edge_weights = [G.edges[e].get("weight", 1) for e in G.edges]
    ew_min, ew_max = (min(edge_weights), max(edge_weights)) if edge_weights else (1, 1)

    def scale_edge(w: int) -> float:
        if ew_max == ew_min:
            return 1.5
        return 0.5 + 3.0 * (w - ew_min) / (ew_max - ew_min)

    widths = [scale_edge(w) for w in edge_weights] if edge_weights else 1.5

    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.set_axis_off()

    nx.draw_networkx_edges(
        G, pos,
        arrowstyle="-|>",
        arrowsize=12,
        width=widths,
        alpha=0.7,
        connectionstyle="arc3,rad=0.08",
    )
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        alpha=0.9,
    )

    if with_labels:
        def short(s: str) -> str:
            return s if len(s) <= label_max_len else s[: label_max_len - 1] + "…"

        labels = {n: f"{short(n)}\n({G.nodes[n]['weight']})" for n in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

    # Edge labels: (k-1)-mer overlaps
    #edge_labels = {(u, v): G.edges[u, v].get("overlap", "") for u, v in G.edges}
    #nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=7, rotate=False)

    plt.tight_layout()
    if out_path:
        plt.savefig(out_path, dpi=200, bbox_inches="tight")
        print(f"Wrote: {out_path}")

    if show:
        plt.show()
    else:
        plt.close()


def main():
    p = argparse.ArgumentParser(description="Visualize a node-centric de Bruijn graph from FASTA.")
    p.add_argument("fasta", help="Input FASTA file")
    p.add_argument("-k", type=int, required=True, help="k-mer size (>=2)")
    p.add_argument("--max-nodes", type=int, default=300,
                   help="Prune to top N nodes by k-mer count (0 = no pruning)")
    p.add_argument("--layout", choices=["spring", "kamada_kawai", "circular"],
                   default="spring", help="Layout algorithm")
    p.add_argument("--labels", action="store_true", help="Show node labels (k-mer + count)")
    p.add_argument("--label-max-len", type=int, default=18, help="Max chars for node labels")
    p.add_argument("--figsize", type=int, nargs=2, default=[14, 10], metavar=("W", "H"))
    p.add_argument("--out", default="debruijn.png", help="Output image path (set empty to skip saving)")
    p.add_argument("--no-show", action="store_true", help="Do not open a window; just save")
    args = p.parse_args()

    seqs = read_fasta_sequences(args.fasta)
    if not seqs:
        raise SystemExit("No sequences found in FASTA.")

    G, _ = build_nodecentric_debruijn(seqs, args.k)
    G = prune_to_top_nodes(G, args.max_nodes)

    out_path = args.out if args.out else None
    draw_graph(
        G,
        out_path=out_path,
        show=not args.no_show,
        with_labels=args.labels,
        label_max_len=args.label_max_len,
        figsize=(args.figsize[0], args.figsize[1]),
        layout=args.layout,
    )


if __name__ == "__main__":
    main()

