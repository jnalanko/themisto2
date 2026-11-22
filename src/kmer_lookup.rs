//! Looking up the colex ranks of the k-mers of a query in an SBWT, without the LCS array.

use sbwt::{SbwtIndex, SubsetMatrix};

/// Finds the colex ranks of the k-mers of query sequences using the SBWT and its (k-1)-suffix
/// group marks.
pub struct KmerLookupIndex<'a> {
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    #[allow(dead_code)] // Used by the streaming implementation (TODO)
    sufgroups: &'a bitvec::vec::BitVec,
}

impl<'a> KmerLookupIndex<'a> {

    /// `sufgroups` are the (k-1)-suffix group marks of `sbwt` (see [sbwt::dbg::Dbg::k_minus_1_marks]).
    pub fn new(sbwt: &'a SbwtIndex<SubsetMatrix>, sufgroups: &'a bitvec::vec::BitVec) -> Self {
        assert_eq!(sufgroups.len(), sbwt.n_sets());
        Self { sbwt, sufgroups }
    }

    pub fn k(&self) -> usize {
        self.sbwt.k()
    }

    pub fn sbwt(&self) -> &'a SbwtIndex<SubsetMatrix> {
        self.sbwt
    }

    /// Pushes the colex rank of every k-mer of `query` to `out`, in order, or None if the k-mer is
    /// not in the index (including k-mers with characters other than ACGT). Pushes nothing if the
    /// query is shorter than k. Does not clear `out`.
    pub fn push_colex_ranks(&self, query: &[u8], out: &mut Vec<Option<usize>>) {
        // TODO: streaming lookup using the suffix group marks. This stub searches every k-mer
        // separately, which is correct but takes O(k) time per k-mer.
        let k = self.k();
        if query.len() < k {
            return;
        }
        for kmer in query.windows(k) {
            let colex = if kmer.iter().all(|&c| crate::util::IS_DNA[c as usize]) {
                self.sbwt.search(kmer).map(|range| {
                    assert_eq!(range.len(), 1);
                    range.start
                })
            } else {
                None
            };
            out.push(colex);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colex_ranks_match_search() {
        let seqs: Vec<&[u8]> = vec![b"ACGTTGCATGCAAGTCCGATAG", b"TGCATGCAAGTCCGTTTTTT"];
        let (mut sbwt, lcs) = sbwt::BitPackedKmerSortingMem::new_from_slices(&seqs, 5).build_lcs(true).run();
        sbwt.build_select();
        let dbg = sbwt::dbg::Dbg::new(&sbwt, lcs.as_ref(), 1);
        let sufgroups = dbg.k_minus_1_marks().clone();
        let index = KmerLookupIndex::new(&sbwt, &sufgroups);

        let query = b"GGACGTTGCANGCAAGTCCGTTTTTTAA";
        let mut ranks = Vec::new();
        index.push_colex_ranks(query, &mut ranks);
        assert_eq!(ranks.len(), query.len() - 5 + 1);
        for (kmer, rank) in query.windows(5).zip(ranks.iter()) {
            let expected = if kmer.contains(&b'N') { None } else { sbwt.search(kmer).map(|r| r.start) };
            assert_eq!(*rank, expected, "{}", String::from_utf8_lossy(kmer));
        }
        assert!(ranks.iter().any(|r| r.is_some()));
        assert!(ranks.iter().any(|r| r.is_none()));
    }
}
