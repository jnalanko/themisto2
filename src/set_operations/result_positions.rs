//! Where the k-mers of the input indexes of a set operation are in its result.

use simple_sds_sbwt::ops::{BitVec, Rank, Select};

use crate::util::u64_bitvec_to_simple_sds;

/// Maps the colex positions of the nodes of one input index to their colex positions in the result
/// of a set operation, through the interleaving that the result consists of: input position i is
/// at interleaving position p = select(s, i), where s is the input's bit vector of the interleaving
/// (s1 or s2). The result consists of the interleaving positions that were not removed, in order,
/// so p is at result position p minus the number of removed positions before it, or not in the
/// result if p itself was removed.
pub(crate) struct InputToResult<'a> {
    in_input: simple_sds_sbwt::bit_vector::BitVector, // s1 or s2 of the interleaving, with select support
    removed: Option<&'a simple_sds_sbwt::bit_vector::BitVector>, // With rank support. None if nothing was removed
}

impl<'a> InputToResult<'a> {
    pub(crate) fn new(in_input: &bitvec::vec::BitVec<u64, bitvec::order::Lsb0>, removed: Option<&'a simple_sds_sbwt::bit_vector::BitVector>) -> Self {
        // TODO: This copies s1 or s2 (n bits for n interleaving positions) only to get select
        // support from simple-sds. A select structure directly on the bitvec of the interleaving
        // would save that copy, at the cost of maintaining our own select code.
        let mut in_input = u64_bitvec_to_simple_sds(in_input);
        in_input.enable_select();
        Self { in_input, removed }
    }

    /// The result colex position of input colex position `input_colex`, or None if it was removed.
    pub(crate) fn result_colex(&self, input_colex: usize) -> Option<usize> {
        let pos = self.in_input.select(input_colex).unwrap();
        match self.removed {
            Some(removed) if removed.get(pos) => None,
            Some(removed) => Some(pos - removed.rank(pos)),
            None => Some(pos),
        }
    }
}
