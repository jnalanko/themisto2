# Themisto 2

Themisto 2 is a **colored k-mer pseudoalignment** tool for microbial genomics and the successor to [Themisto 1](https://github.com/algbio/Themisto). It extends Themisto 1 with:

* **A low-memory and disk-free** [construction algorithm](https://www.biorxiv.org/content/10.64898/2026.02.16.706153v1) from the SBWT.
* **Faster queries** enabled by streaming k-mer search via the LCS array.
* **Index updates** via index merging.
* **Export and import** of the index in a text-based colored unitig format.
* **Easy installation** from source via the Rust package manager Cargo. The build toolchain can be installed without the need for root privileges.

Themisto 2 builds a colored k-mer index over a collection of reference sequences. Each input file is assigned a color (an integer identifier). For every distinct k-mer in the collection, the index records the set of colors, i.e. the set of input files that contain it.

Given a query sequence, **pseudoalignment** looks up each of the query's k-mers in the color matrix and aggregates the results across all k-mers:

- **Intersection pseudoalignment** (`intersection-pseudoalign`) reports the colors present in *every* k-mer of the query that matches the index.
- **Threshold pseudoalignment** (`threshold-pseudoalign`) reports colors that appear in at least a given fraction of the query's k-mers, controlled by `--threshold`. This is more robust to sequencing errors, small variants and incomplete references, but is less sensitive to small strain-level variation.

Both modes output one JSONL record per query sequence, with the sequence name and the list of matching color indices.

# Installation

First, install Rust and the Cargo package manager using rustup (no root privileges required): [https://rust-lang.org/tools/install/](https://rust-lang.org/tools/install/). Then, clone this repository and compile with: 

```
git clone --recursive https://github.com/jnalanko/themisto2`
cd themisto2
cargo build --release
```

The executable will be stored at `./target/release/themisto2`.
# Quick start: Constructing small indices

For small datasets less than a few gigabytes in size, you can easily construct the index with a single command. You can try it with the included example data (run from the repository root):

```
themisto2 build --file-colors example/fof.txt -o index.thm2 --temp-dir temp -k 5 -t 4
```

Here `example/fof.txt` is a file with one fasta/fastq filename per line (each file represents one color), `index.thm2` is the output index, `temp` is a directory for temporary files, and `-k 5` sets the k-mer length. The first file in the input list is assigned to color id 0, the next file to color id 1, and so on.

Once the index is built, you can pseudoalign queries against it:

```
themisto2 threshold-pseudoalign --threshold 0.7 --denominator all --min-hits 1 -i index.thm2 -q example/C1.fna -t 4
```

This outputs to the standard output, for each query sequence, the set of colors whose who have at least 70% of the k-mers of the query:

```
{"name": "C1.1", "colors": [0,1]}
{"name": "C1.2", "colors": [0,1,2]}
```
# Constructing large indices

For large indexing runs (more than a few gigabytes in size) with repetitive data, it is important to deduplicate the k-mers first. The final index will be the same as if the index was built directly as described above, but the running time, peak RAM and disk can be orders of magnitude better.

Here we describe a construction pipeline that should suit most large microbial genomics indexing tasks. The pipeline requires two external tools. Pleases install these two before proceeding:

* [GGCAT](https://github.com/algbio/ggcat)
* [SBWT-rs-cli](https://github.com/jnalanko/sbwt-rs-cli)
## Step 1: Unitigs

Deduplicate the k-mers by reducing them to de Bruijn graph unitigs. We recommend using [GGCAT](https://github.com/algbio/ggcat) for this, but tools like BCALM2 and Cuttlefish 3 will also work. You can run GGCAT as follows:

```
ggcat build -l input_file_list.txt -o unitigs-k31.fna -s 1 -k 31 -t temp -j 32 -m 64 -p
```

This uses 64 GiB of RAM and 32 threads, with k = 31. The option `-s 1` is very important: otherwise GGCAT  discards k-mers occurring only once. Giving more RAM generally speeds up the construction.

## Step 2: SBWT

Built the SBWT using the SBWT Rust tool.

```
sbwt build -i unitigs-k31.fna -t 32 -m 64 -v -o unitigs-k31 -l -k 31 -r --temp-dir temp
```

This uses 32 threads and 64 GiB memory (the memory usage may be larger if the final index is larger than 64 GiB). The option `-o unitigs-k31` is the prefix for the output files. The directory `temp` will be used for temporary storage. This should point to a fast file system with a lot of disk space (up to 8n bytes per distinct k-mer). Two files will be built: `unitigs-k31.sbwt` and `unitigs-k31.lcs`.
## Step 3: Build the Themisto 2 index

Build the Themisto 2 index using the SBWT and LCS array.

```
themisto2 build --file-colors input_file_list.txt -s unitigs-k31.sbwt -l unitigs-k31.lcs -o index.thm2 -k 31 -t 32
```

This will use 32 threads.

# Citation

A Themisto 2 paper is forthcoming, but in the meantime, you can cite the preprint describing the construction algorithm:

```
@article{alanko2026construction,
  title={Construction of distinct k-mer color sets via set fingerprinting},
  author={Alanko, Jarno N and Puglisi, Simon J},
  journal={bioRxiv},
  pages={2026--02},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}
```

# Related tools

Other tools implementing k-mer pseudoalignment include [Bifrost](https://github.com/pmelsted/bifrost), [Metagraph](https://github.com/ratschlab/metagraph) and [Fulgor](https://github.com/jermp/fulgor).  Themisto 2 is at the cutting edge of BWT-based k-mer indexing, benefitting from latest research on SBWT. This enables in particular, fast and very parallel index merging.

See also: [HKS](https://github.com/jnalanko/HKS), an SBWT-based sequence annotation tool with hierarchical labels.