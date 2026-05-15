use std::path::PathBuf;

use crossbeam::channel::Receiver;
use jseqio::reader::DynamicFastXReader;
use sbwt::{reverse_complement_in_place, SeqStream};

pub struct ChainedInputStream{
    paths: Vec<PathBuf>,
    cur_file: Option<DynamicFastXReader>,
    seq_buf: Vec<u8>, // Local buffer from which we can borrow (can not use the buffer of cur_file for lifetime reasons)
    cur_file_idx : usize, // Index of the db currently being iterated over
}

impl ChainedInputStream {
    pub fn new(filenames: Vec<PathBuf>) -> Self {
        let first_file = filenames.first().map(|f| DynamicFastXReader::from_file(f).unwrap());
        Self {paths: filenames, cur_file: first_file, seq_buf: vec![], cur_file_idx: 0}
    }

    #[allow(dead_code)]
    pub fn cur_file_idx(&self) -> usize {
        self.cur_file_idx
    }

    #[allow(dead_code)]
    pub fn get_seq_buf(&self) -> &[u8] {
        &self.seq_buf
    }

    #[allow(dead_code)]
    pub fn get_seq_buf_mut(&mut self) -> &mut [u8] {
        &mut self.seq_buf
    }

    #[allow(dead_code)]
    pub fn done(&self) -> bool {
        self.cur_file_idx == self.paths.len()
    }
}

impl SeqStream for ChainedInputStream {
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.seq_buf.clear();
        if let Some(f) = self.cur_file.as_mut() {
            if let Some(rec) = f.read_next().unwrap() {
                self.seq_buf.extend_from_slice(rec.seq);
                Some(&self.seq_buf)
            } else {
                // File is finished -> open the next file
                self.cur_file_idx += 1;
                self.cur_file = if self.cur_file_idx == self.paths.len() {
                    None // All files procesed
                } else {
                    let new_file = DynamicFastXReader::from_file(&self.paths[self.cur_file_idx]).unwrap();
                    Some(new_file)
                };

                self.stream_next()
            }
        } else {
            None    
        }
    }
}

pub struct ChainedInputStreamWithRevComp{
    inner: ChainedInputStream,
    rev_comp_next: bool,
}

impl ChainedInputStreamWithRevComp {
    pub fn new(filenames: Vec<PathBuf>) -> Self {
        let inner = ChainedInputStream::new(filenames);
        Self{inner, rev_comp_next: false}
    }
}

impl SeqStream for ChainedInputStreamWithRevComp {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.rev_comp_next {
            reverse_complement_in_place(&mut self.inner.seq_buf); 
            self.rev_comp_next = false;
            Some(&self.inner.seq_buf)
        } else {
            self.rev_comp_next = true;
            self.inner.stream_next()
        }
    }
}

impl ChainedInputStreamWithRevComp {
    #[allow(dead_code)]
    pub fn cur_file_idx(&self) -> usize {
        self.inner.cur_file_idx
    }

    #[allow(dead_code)]
    pub fn get_seq_buf(&self) -> &[u8] {
        &self.inner.seq_buf
    }

    #[allow(dead_code)]
    pub fn get_seq_buf_mut(&mut self) -> &mut [u8] {
        &mut self.inner.seq_buf
    }

    #[allow(dead_code)]
    pub fn done(&self) -> bool {
        self.inner.cur_file_idx == self.inner.paths.len()
    }
}

// ── Per-color sequence source ────────────────────────────────────────────────

/// One color's sequence input, either a whole file or a single in-memory sequence.
#[derive(Clone)]
pub enum ColorSource {
    File(PathBuf),
    SingleSeq(Vec<u8>),
}

impl ColorSource {
    pub fn open(&self) -> ColorSourceIter {
        match self {
            Self::File(p) => ColorSourceIter::File(DynamicFastXReader::from_file(p).unwrap()),
            Self::SingleSeq(seq) => ColorSourceIter::Single { seq: seq.clone(), done: false },
        }
    }
}

/// Iterator yielding owned sequences for one color.
pub enum ColorSourceIter {
    File(DynamicFastXReader),
    Single { seq: Vec<u8>, done: bool },
}

impl Iterator for ColorSourceIter {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Vec<u8>> {
        match self {
            Self::File(reader) => reader.read_next().unwrap().map(|rec| rec.seq.to_vec()),
            Self::Single { seq, done } => {
                if *done { None } else { *done = true; Some(std::mem::take(seq)) }
            }
        }
    }
}

// ── SeqStream over all color sources (for SBWT / key-kmer phases) ────────────

pub struct AllColorSeqs {
    rx: Receiver<ColorSource>,
    cur_iter: Option<ColorSourceIter>,
    seq_buf: Vec<u8>,
    include_rev_comp: bool,
    rev_comp_pending: bool,
}

impl AllColorSeqs {
    pub fn new(rx: Receiver<ColorSource>, include_rev_comp: bool) -> Self {
        Self { rx, cur_iter: None, seq_buf: vec![], include_rev_comp, rev_comp_pending: false }
    }
}

impl SeqStream for AllColorSeqs {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.rev_comp_pending {
            reverse_complement_in_place(&mut self.seq_buf);
            self.rev_comp_pending = false;
            return Some(&self.seq_buf);
        }

        loop {
            let got_seq = if let Some(it) = self.cur_iter.as_mut() {
                if let Some(seq) = it.next() {
                    self.seq_buf = seq;
                    true
                } else {
                    false
                }
            } else {
                false
            };

            if got_seq {
                self.rev_comp_pending = self.include_rev_comp;
                return Some(&self.seq_buf);
            }

            match self.rx.recv() {
                Ok(source) => self.cur_iter = Some(source.open()),
                Err(_) => return None, // sender dropped, channel exhausted
            }
        }
    }
}