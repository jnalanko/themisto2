use std::path::PathBuf;

use jseqio::{reader::DynamicFastXReader, reverse_complement};
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

    pub fn cur_file_idx(&self) -> usize {
        self.cur_file_idx
    }

    pub fn get_seq_buf(&self) -> &[u8] {
        &self.seq_buf
    }

    pub fn get_seq_buf_mut(&mut self) -> &mut [u8] {
        &mut self.seq_buf
    }

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
    pub fn cur_file_idx(&self) -> usize {
        self.inner.cur_file_idx
    }

    pub fn get_seq_buf(&self) -> &[u8] {
        &self.inner.seq_buf
    }

    pub fn get_seq_buf_mut(&mut self) -> &mut [u8] {
        &mut self.inner.seq_buf
    }

    pub fn done(&self) -> bool {
        self.inner.cur_file_idx == self.inner.paths.len()
    }
}