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

use sbwt::{dbg::Dbg, SubsetMatrix};

use crate::{colex_colored_kmers::DbgIndex, coloring_interface::ColorSetStorage};

/// The DBG of the SBWT of `coloring`, an input of a set operation, from its stored marks.
pub(crate) fn input_dbg<CSS: ColorSetStorage>(coloring: &DbgIndex<CSS>) -> Dbg<'_, SubsetMatrix> {
    coloring.dbg()
}
