use std::ops::Range;

use jseqio::seq_db::SeqDB;
use sbwt::SeqStream;

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}};

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

#[derive(Clone)]
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

struct PseudoalignmentBatch {
    seqs: SeqDB,
    seq_ranks: Range<usize>
}

impl PseudoalignmentBatch {
    fn new(first_seq_rank: usize) -> Self {
        Self { seqs: SeqDB::new(), seq_ranks: first_seq_rank..first_seq_rank }
    }

    fn push(&mut self, seq: &[u8]) {
        self.seqs.push_seq(seq);
        self.seq_ranks.end += 1;
    }
}

fn run_all_queries<CSS: ColorSetStorage>(index: &CompactColexKmers<CSS>, mut queries: impl SeqStream + Send + 'static, algorithm: Box<dyn PseudoalignmentAlgorithm<CSS>>, n_workers: usize) {
    let mut aligner = ThresholdPseudoaligner {
        counts: vec![0; index.get_set_storage().n_colors()],
        nonzero_count_indices: vec![],
        threshold: 0.7,
        denominator: Denominator::Relevant,
        answer: index.get_set_storage().get_empty_set(),
        storage: index.get_set_storage(),
    };

    let batch_size = 10_000_usize;
    let (work_send, work_recv) = crossbeam::channel::bounded::<PseudoalignmentBatch>(10);

    std::thread::scope(|scope| {
        let parser_handle = scope.spawn(move || {
            let mut n_seqs_read = 0_usize;
            let mut cur_batch = PseudoalignmentBatch::new(0);
            while let Some(q) = queries.stream_next() {
                n_seqs_read += 1;
                cur_batch.push(q);
                if cur_batch.seqs.total_seq_len() >= batch_size {
                    work_send.send(cur_batch).unwrap();
                    cur_batch = PseudoalignmentBatch::new(n_seqs_read);
                }
            }
            if cur_batch.seqs.total_seq_len() > 0 { // Last batch
                work_send.send(cur_batch).unwrap();
            }
            drop(work_send);
        });

        let mut worker_handles = vec![];
        for worker_id in 0..n_workers {
            let handle = scope.spawn(|| {
                while let Ok(batch) = work_recv.recv() {

                }
            });
            worker_handles.push(handle);
        }
    }); 

    let ans1 = aligner.process(b"ACGTAGCTGAC");
    for c in ans1.iter() {
        println!("{}", c);
    }
    drop(ans1);
    let ans2 = aligner.process(b"ACAATGCTGATCA");
    for c in ans2.iter() {
        println!("{}", c);
    }
}