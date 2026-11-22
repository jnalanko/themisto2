//! Set operations on colored indexes: merging (union), intersecting and subtracting two indexes,
//! and taking the subset of an index for some of its colors.

pub(crate) mod colors;
pub(crate) mod key_kmers;
pub(crate) mod result_positions;
pub(crate) mod color_sets;
pub(crate) mod pipeline;
pub mod merge;
pub mod intersect;
pub mod difference;
pub mod subset;
#[cfg(test)]
mod test_utils;

