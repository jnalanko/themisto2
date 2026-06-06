// Benchmark the I/O and producer-consumer overhead of MsElementGenerator,
// with the matching-statistics work replaced by a no-op.
//
// Usage: cargo run --release --bin bench_producer_consumer_io -- <file_of_files> [n_threads]
//
// <file_of_files> is a text file with one fasta/fastq path per line, matching
// the --file-colors input format in main.rs.  Each file is one color.
//
// Reads through the same batching pipeline used by MsElementGenerator
// (8 MiB batches, crossbeam channel, n_threads consumers) but the consumers do
// nothing — batches are immediately dropped.  The throughput reported is the
// ceiling imposed by I/O + batching alone.

mod io;

use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering::Relaxed},
    time::Instant,
};

use io::RewindableSeqStreamGenerator;
use sbwt::SeqStream as _;

struct Batch {
    seq_concat: Vec<u8>,
    seq_ends: Vec<usize>,
    colors: Vec<usize>,
}

impl Batch {
    fn new() -> Self {
        Self { seq_concat: vec![], seq_ends: vec![], colors: vec![] }
    }

    fn push(&mut self, seq: &[u8], color: usize) {
        self.seq_concat.extend_from_slice(seq);
        self.seq_ends.push(self.seq_concat.len());
        self.colors.push(color);
    }

    fn size_in_bytes(&self) -> usize {
        self.seq_concat.len() + (self.seq_ends.len() + self.colors.len()) * size_of::<usize>()
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <file_of_files> [n_threads]", args[0]);
        std::process::exit(1);
    }
    let fof_path = PathBuf::from(&args[1]);
    let n_threads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);

    // One fasta/fastq file per line, one file per color — matching --file-colors in main.rs.
    let input_paths: Vec<PathBuf> = BufReader::new(File::open(&fof_path).unwrap())
        .lines()
        .map(|l| PathBuf::from(l.unwrap()))
        .collect();
    let n_colors = input_paths.len();

    let mut gen: Box<dyn RewindableSeqStreamGenerator + Sync + Send> =
        Box::new(io::SeqStreamGeneratorFromFiles::new(input_paths));

    let total_bases = AtomicU64::new(0);
    let total_seqs = AtomicU64::new(0);
    let total_bases_ref = &total_bases;
    let total_seqs_ref = &total_seqs;

    let t0 = Instant::now();

    let (sender, receiver) = crossbeam::channel::bounded::<Batch>(2 * n_threads);
    let receiver_ref = &receiver;

    std::thread::scope(|scope| {
        let producer = scope.spawn(|| {
            const BATCH_SIZE: usize = 1 << 23; // 8 MiB, same as MsElementGenerator
            let mut color = 0_usize;
            let mut cur_batch = Batch::new();
            while let Some(mut stream) = gen.next() {
                while let Some(seq) = stream.stream_next() {
                    cur_batch.push(seq, color);
                    if cur_batch.size_in_bytes() >= BATCH_SIZE {
                        sender.send(cur_batch).unwrap();
                        cur_batch = Batch::new();
                    }
                }
                color += 1;
            }
            sender.send(cur_batch).unwrap(); // flush last (possibly empty) batch
            drop(sender);
        });

        let consumers: Vec<_> = (0..n_threads)
            .map(|_| {
                scope.spawn(|| {
                    while let Ok(batch) = receiver_ref.recv() {
                        total_bases_ref.fetch_add(batch.seq_concat.len() as u64, Relaxed);
                        total_seqs_ref.fetch_add(batch.seq_ends.len() as u64, Relaxed);
                        // No work: batch is dropped here
                    }
                })
            })
            .collect();

        producer.join().unwrap();
        for h in consumers {
            h.join().unwrap();
        }
    });

    let elapsed = t0.elapsed();
    let total_bases = total_bases.load(Relaxed);
    let total_seqs = total_seqs.load(Relaxed);
    let secs = elapsed.as_secs_f64();
    let gib = total_bases as f64 / (1u64 << 30) as f64;

    println!("Threads:      {n_threads}");
    println!("Colors:       {n_colors}");
    println!("Sequences:    {total_seqs}");
    println!("Total bases:  {total_bases} ({gib:.3} GiB)");
    println!("Time:         {secs:.3} s");
    println!("Throughput:   {:.2} GiB/s", gib / secs);
}
