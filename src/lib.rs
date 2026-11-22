//! Library interface to Themisto 2 indexes.

#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material
#![allow(clippy::len_zero, clippy::len_without_is_empty)] // !is_empty reads as "not is empty" which is not English 
#![allow(clippy::manual_is_multiple_of)] // Oh please

// We assume that usize is 64 bits in many places
#[cfg(not(target_pointer_width = "64"))]
compile_error!("This crate requires a 64-bit usize (target_pointer_width = 64).");

// Public API
pub use colex_colored_kmers::{Coloring, DbgIndex, IndexStats};
pub use coloring_interface::{ColorSetStorage, ColorSetView};
pub use sparse_dense_storage::SparseDenseStorage;
pub use index_dir::{Component, IndexDir, IndexWriter, Metadata};
pub use kmer_lookup::KmerLookupIndex;
pub use pseudoalignment::PseudoalignmentIndex;

pub mod index_dir;
pub mod kmer_lookup;
pub mod set_operations;

// Internals. Public so that the themisto2 binary can use them, but hidden from the docs.
#[doc(hidden)] pub mod EM;
#[doc(hidden)] pub mod bitmap_storage;
#[doc(hidden)] pub mod index_import;
#[doc(hidden)] pub mod compatibility_criteria;
#[doc(hidden)] pub mod colex_colored_kmers;
#[doc(hidden)] pub mod coloring_interface;
#[doc(hidden)] pub mod sparse_dense_storage;
#[doc(hidden)] pub mod io;
#[doc(hidden)] pub mod set_of_sets_construction;
#[doc(hidden)] pub mod iterators;
#[doc(hidden)] pub mod parallel_ms_iteration;
#[doc(hidden)] pub mod atomic_bitmap;
#[doc(hidden)] pub mod int_vec;
#[doc(hidden)] pub mod util;
#[doc(hidden)] pub mod pseudoalignment;
#[doc(hidden)] pub mod pseudoalignment_metrics;
#[doc(hidden)] pub mod sparse_dense_storage_to_disk;
#[doc(hidden)] pub mod report;
#[doc(hidden)] pub mod filter_reads;
#[doc(hidden)] pub mod work_dispatcher;
#[doc(hidden)] pub mod unitig_export;
