use std::sync::{Arc, Mutex};

use crate::io::RewindableSeqStreamGenerator;

pub trait Worker {
    fn process(&mut self, seq: &[u8], color: usize); // Takes pairs (sequence, sequence color)
}

pub fn dispatch_work(
    seq_stream_gen: &mut Box<dyn RewindableSeqStreamGenerator + Sync + Send>,
    workers: Vec<impl Worker + Sync + Send>,
    n_parsers: usize,
    batch_flush_threshold: usize
) {

    let (sender, receiver) = crossbeam::channel::bounded::<WorkBatch>(2 * workers.len());
    let receiver_ref = &receiver;

    let gen_mutex = Arc::new(Mutex::new(seq_stream_gen));

    std::thread::scope(|scope| {

        let mut producers = vec![];
        for _ in 0..n_parsers { // Create parsers 
            let sender_clone = sender.clone();
            let gen_mutex_clone = gen_mutex.clone();
            let handle = scope.spawn(move || {
                let mut cur_batch = WorkBatch::new();
                loop {
                    let next_stream_maybe = { 
                        let mut gen = gen_mutex_clone.lock().unwrap();
                        gen.next()
                    }; // gen_mutex_clone lock is freed here

                    if let Some((mut stream, color)) = next_stream_maybe {
                        log::info!("Processing color {color}");
                        while let Some(seq) = stream.stream_next() {
                            cur_batch.push_seq(seq, color);
                            if cur_batch.size_in_bytes() >= batch_flush_threshold {
                                sender_clone.send(cur_batch).unwrap();
                                cur_batch = WorkBatch::new();
                            }
                        }
                    } else {
                        sender_clone.send(cur_batch).unwrap(); // flush last (possibly empty) batch
                        break;
                    }
                }
                drop(sender_clone); // Drop our copy of the generator mutex. When all copies are dropped, the channel is closed
            });
            producers.push(handle);
        }
        drop(sender); // Drop this clone of the sender. Now the only copies are within the producers.

        let consumers: Vec<_> = workers.into_iter().map(|mut worker| {
            scope.spawn(move || { // Takes ownership of the worker
                while let Ok(batch) = receiver_ref.recv() {
                    batch.process(&mut worker);
                }
            })
        }).collect();


        for c in producers {
            c.join().unwrap();
        }
        for h in consumers {
            h.join().unwrap();
        }
    });
}

struct WorkBatch {
    seq_concat: Vec<u8>,
    seq_ends: Vec<usize>,
    seq_colors: Vec<usize>,
}

impl WorkBatch {

    fn push_seq(&mut self, seq: &[u8], color: usize) {
        self.seq_concat.extend_from_slice(seq);
        self.seq_ends.push(self.seq_concat.len());
        self.seq_colors.push(color);
    }

    fn size_in_bytes(&self) -> usize {
        self.seq_concat.len() + self.seq_ends.len()*size_of::<usize>() + self.seq_colors.len()*size_of::<usize>()
    }

    fn new() -> WorkBatch {
        WorkBatch { seq_concat: vec![], seq_ends: vec![], seq_colors: vec![] }
    }

    // Process the batch with the given worker
    fn process(self, worker: &mut impl Worker) {
        let n_seqs = self.seq_ends.len();
        assert_eq!(n_seqs, self.seq_colors.len());

        let mut seq_start = 0_usize;
        for seq_idx in 0..n_seqs {
            let seq_end = self.seq_ends[seq_idx];
            let seq = &self.seq_concat[seq_start..seq_end];
            let color = self.seq_colors[seq_idx];
            worker.process(seq, color);
            seq_start = seq_end;
        }
    }
}
