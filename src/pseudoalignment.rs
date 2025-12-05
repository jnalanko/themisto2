use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage};

trait PseudoalignmentAlgorithm<CSS: ColorSetStorage> {

    // The &mut self is to access and modify thread-local buffers
    // owned by the algorithm.
    fn process<'a>(&'a mut self, seq: &[u8]) -> CSS::SetView<'a>;
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Denominator { // Options for the CLI
    All,
    Relevant,
    MaxHits,
}

struct ThresholdPseudoaligner<'a, CSS: ColorSetStorage> {
    counts: Vec<usize>,
    nonzero_count_indices: Vec<usize>,
    threshold: f64,
    denominator: Denominator,
    answer: CSS::OwnedSet, // Answer to the current query
    storage: &'a CSS,
}

impl<'a, CSS: ColorSetStorage> PseudoalignmentAlgorithm<CSS> for ThresholdPseudoaligner<'a, CSS> {
    fn process<'b>(&'b mut self, seq: &[u8]) -> <CSS as ColorSetStorage>::SetView<'b> {
        self.storage.owned_to_view(&self.answer)
    }
}

fn test <CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>) {
    let mut aligner = ThresholdPseudoaligner {
        counts: vec![0; index.get_set_storage().n_colors()],
        nonzero_count_indices: vec![],
        threshold: 0.7,
        denominator: Denominator::Relevant,
        storage: index.get_set_storage(),
    };
    let ans1 = aligner.process(b"ACGTAGCTGAC");
    let ans2 = aligner.process(b"ACGTAGCTGAC");
}