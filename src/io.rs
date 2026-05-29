use std::path::PathBuf;

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

pub trait RewindableSeqStreamGenerator {
	fn next(&mut self) -> Option<Box<dyn SeqStream + Send + Sync>>;
	fn rewind(&mut self);
}

pub struct SeqStreamGeneratorFromFiles {
    files: Vec<PathBuf>,
    cur_file_idx: usize,
}

impl SeqStreamGeneratorFromFiles {
    pub fn new(files: Vec<PathBuf>) -> Self {
        Self {files, cur_file_idx: 0}
    }
}

pub struct JSeqIOWrapper { // So that we can implement sbwt::SeqStream for jseqio::reader
    inner: jseqio::reader::DynamicFastXReader,
    cur_buf: Vec<u8>,
}

impl sbwt::SeqStream for JSeqIOWrapper {
    fn stream_next(&mut self) -> Option<&[u8]> {
        let maybe_rec = self.inner.read_next().unwrap(); // Unwrap the IO Result
        let rec = maybe_rec?; // If None -> end of stream
        self.cur_buf.clear();
        self.cur_buf.extend_from_slice(rec.seq);
        Some(&self.cur_buf)
    }
}

impl RewindableSeqStreamGenerator for SeqStreamGeneratorFromFiles {
    fn next(&mut self) -> Option<Box<dyn SeqStream + Send + Sync>> {
        if self.cur_file_idx == self.files.len() { return None; }

        let reader = jseqio::reader::DynamicFastXReader::from_file(&self.files[self.cur_file_idx]).unwrap();
        let reader = JSeqIOWrapper {inner: reader, cur_buf: vec![]};
        let reader: Box<dyn SeqStream + Send+ Sync> = Box::new(reader);

        self.cur_file_idx += 1;
        Some(reader)
    }

    fn rewind(&mut self) {
        self.cur_file_idx = 0;
    }
}

pub struct SeqStreamGeneratorFromSingleFile {
    file: PathBuf,
    cur_stream: jseqio::reader::DynamicFastXReader,
}

impl SeqStreamGeneratorFromSingleFile {
    pub fn new(file: PathBuf) -> Self {
        let cur_stream = jseqio::reader::DynamicFastXReader::from_file(&file).unwrap();
        Self {file, cur_stream}
    }
}

pub struct SingleSeqStream {
    seq: Vec<u8>,
    done: bool,
}

impl SeqStream for SingleSeqStream {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.done { return None };

        self.done = true;
        Some(&self.seq)
    }
}

impl RewindableSeqStreamGenerator for SeqStreamGeneratorFromSingleFile {
    fn next(&mut self) -> Option<Box<dyn SeqStream + Sync + Send>> {
        let rec = self.cur_stream.read_next().unwrap()?;
        let seq = SingleSeqStream { seq: rec.seq.to_vec(), done: false };
        let seq: Box<dyn SeqStream + Sync + Send> = Box::new(seq);
        Some(seq)
    }

    fn rewind(&mut self) {
        let mut new_reader = jseqio::reader::DynamicFastXReader::from_file(&self.file).unwrap();
        std::mem::swap(&mut self.cur_stream, &mut new_reader);

        // Thd old reader is dropped here
    }
}

// Chains multiple fasta/fastq files (gzipped or not), uppercases sequences,
// and emits each sequence followed by its reverse complement.
pub struct NeedletailSeqStreamWithRevComp {
    paths: Vec<PathBuf>,
    cur_idx: usize,
    cur_reader: Option<Box<dyn needletail::FastxReader>>,
    seq_buf: Vec<u8>,
    rev_comp_next: bool,
}

impl NeedletailSeqStreamWithRevComp {
    pub fn new(paths: Vec<PathBuf>) -> Self {
        let cur_reader = paths.first().map(|p| needletail::parse_fastx_file(p).unwrap());
        Self { paths, cur_idx: 0, cur_reader, seq_buf: vec![], rev_comp_next: false }
    }
}

impl SeqStream for NeedletailSeqStreamWithRevComp {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.rev_comp_next {
            reverse_complement_in_place(&mut self.seq_buf);
            self.rev_comp_next = false;
            return Some(&self.seq_buf);
        }
        loop {
            let has_record = match self.cur_reader.as_mut() {
                None => return None,
                Some(reader) => match reader.next() {
                    None => false,
                    Some(Err(e)) => panic!("Error reading sequence: {e}"),
                    Some(Ok(rec)) => {
                        let seq = rec.seq();
                        self.seq_buf.clear();
                        self.seq_buf.extend(seq.iter().map(|&b| b.to_ascii_uppercase()));
                        true
                    }
                }
            };
            if has_record {
                self.rev_comp_next = true;
                return Some(&self.seq_buf);
            }
            self.cur_idx += 1;
            self.cur_reader = self.paths.get(self.cur_idx)
                .map(|p| needletail::parse_fastx_file(p).unwrap());
        }
    }
}

pub struct EmptyRewindableSeqStreamGenerator { // Generates nothing

}

impl RewindableSeqStreamGenerator for EmptyRewindableSeqStreamGenerator {
    fn next(&mut self) -> Option<Box<dyn SeqStream + Send + Sync>> {
        None
    }

    fn rewind(&mut self) {}
}