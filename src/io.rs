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
    #[allow(dead_code)]
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
    // Gives a stream and the index of the stream
	fn next(&mut self) -> Option<(Box<dyn SeqStream + Send + Sync>, usize)>;
    #[allow(dead_code)]
	fn rewind(&mut self);
}

// Generates one stream per group of files. The files of a group are read one after
// the other as a single stream.
pub struct SeqStreamGeneratorFromFiles {
    file_groups: Vec<Vec<PathBuf>>,
    cur_group_idx: usize,
}

impl SeqStreamGeneratorFromFiles {
    pub fn new(file_groups: Vec<Vec<PathBuf>>) -> Self {
        Self {file_groups, cur_group_idx: 0}
    }
}

impl RewindableSeqStreamGenerator for SeqStreamGeneratorFromFiles {
    fn next(&mut self) -> Option<(Box<dyn SeqStream + Send + Sync>, usize)> {
        if self.cur_group_idx == self.file_groups.len() { return None; }

        let reader = ChainedInputStream::new(self.file_groups[self.cur_group_idx].clone());
        let reader: Box<dyn SeqStream + Send+ Sync> = Box::new(reader);

        let stream_idx = self.cur_group_idx;
        self.cur_group_idx += 1;
        Some((reader, stream_idx))
    }

    fn rewind(&mut self) {
        self.cur_group_idx = 0;
    }
}

/// Checks that `record_colors` (the color of each record of a file, in file order) has the
/// records of each color consecutive, and the colors numbered 0, 1, 2, ... in file order.
/// Otherwise returns the index of the first record that breaks this.
pub fn check_consecutive_colors(record_colors: &[usize]) -> Result<(), usize> {
    let mut expected_next = 0; // The color of a record must be that of the previous record, or this
    for (record_idx, &color) in record_colors.iter().enumerate() {
        if color + 1 == expected_next {
            continue; // Same color as the previous record
        }
        if color != expected_next {
            return Err(record_idx);
        }
        expected_next += 1;
    }
    Ok(())
}

// Generates one stream per color from a single file, where record_colors[i] is the color
// of the i-th record. The records of each color must be consecutive in the file, and colors
// numbered in file order (see check_consecutive_colors); otherwise this panics.
// The stream of a color holds all of its records in memory, so that consumers can read the
// streams of different colors in parallel.
pub struct SeqStreamGeneratorFromSingleFile {
    file: PathBuf,
    cur_stream: jseqio::reader::DynamicFastXReader,
    record_colors: Vec<usize>,
    n_seqs_read: usize,
}

impl SeqStreamGeneratorFromSingleFile {
    pub fn new(file: PathBuf, record_colors: Vec<usize>) -> Self {
        if let Err(record_idx) = check_consecutive_colors(&record_colors) {
            panic!("Sequence {} of {} has a color whose sequences are not consecutive in the file. All sequences of a color must be consecutive.", record_idx + 1, file.display());
        }
        let cur_stream = jseqio::reader::DynamicFastXReader::from_file(&file).unwrap();
        Self {file, cur_stream, record_colors, n_seqs_read: 0}
    }
}

impl RewindableSeqStreamGenerator for SeqStreamGeneratorFromSingleFile {
    fn next(&mut self) -> Option<(Box<dyn SeqStream + Sync + Send>, usize)> {
        let color = *self.record_colors.get(self.n_seqs_read)?;

        let mut seqs = Vec::<Vec<u8>>::new();
        while self.record_colors.get(self.n_seqs_read) == Some(&color) {
            let rec = self.cur_stream.read_next().unwrap()
                .unwrap_or_else(|| panic!("{} has fewer sequences than expected", self.file.display()));
            seqs.push(rec.seq.to_vec());
            self.n_seqs_read += 1;
        }

        Some((Box::new(crate::util::VecVecSeqStream::new(seqs)), color))
    }

    fn rewind(&mut self) {
        let mut new_reader = jseqio::reader::DynamicFastXReader::from_file(&self.file).unwrap();
        self.n_seqs_read = 0;
        std::mem::swap(&mut self.cur_stream, &mut new_reader);

        // Thd old reader is dropped here
    }
}

// Chains multiple fasta/fastq files (gzipped or not), uppercases sequences,
// and emits each sequence followed by its reverse complement.
#[allow(dead_code)]
pub struct NeedletailSeqStreamWithRevComp {
    paths: Vec<PathBuf>,
    cur_idx: usize,
    cur_reader: Option<Box<dyn needletail::FastxReader>>,
    seq_buf: Vec<u8>,
    rev_comp_next: bool,
}

impl NeedletailSeqStreamWithRevComp {
    #[allow(dead_code)]
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

// Chains multiple fasta/fastq files (gzipped or not), uppercases sequences.
#[allow(dead_code)]
pub struct NeedletailSeqStream {
    paths: Vec<PathBuf>,
    cur_idx: usize,
    cur_reader: Option<Box<dyn needletail::FastxReader>>,
    seq_buf: Vec<u8>,
}

impl NeedletailSeqStream {
    #[allow(dead_code)]
    pub fn new(paths: Vec<PathBuf>) -> Self {
        let cur_reader = paths.first().map(|p| needletail::parse_fastx_file(p).unwrap());
        Self { paths, cur_idx: 0, cur_reader, seq_buf: vec![] }
    }
}

impl SeqStream for NeedletailSeqStream {
    fn stream_next(&mut self) -> Option<&[u8]> {
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
                return Some(&self.seq_buf);
            }
            self.cur_idx += 1;
            self.cur_reader = self.paths.get(self.cur_idx)
                .map(|p| needletail::parse_fastx_file(p).unwrap());
        }
    }
}

#[allow(dead_code)]
pub struct EmptyRewindableSeqStreamGenerator { // Generates nothing

}

impl RewindableSeqStreamGenerator for EmptyRewindableSeqStreamGenerator {
    fn next(&mut self) -> Option<(Box<dyn SeqStream + Send + Sync>, usize)> {
        None
    }

    fn rewind(&mut self) {}
}
#[cfg(test)]
mod tests {
    use super::*;

    // Writes the sequences as a fasta file in the temp directory of the project and returns its path
    fn write_fasta(name: &str, seqs: &[&str]) -> PathBuf {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("temp").join("io-unit-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        let content: String = seqs.iter().enumerate().map(|(i, s)| format!(">seq{}\n{}\n", i, s)).collect();
        std::fs::write(&path, content).unwrap();
        path
    }

    // Returns the (color, sequences) pairs generated by the generator
    fn collect(gen: &mut SeqStreamGeneratorFromSingleFile) -> Vec<(usize, Vec<String>)> {
        let mut out = vec![];
        while let Some((mut stream, color)) = gen.next() {
            let mut seqs = vec![];
            while let Some(seq) = stream.stream_next() {
                seqs.push(String::from_utf8(seq.to_vec()).unwrap());
            }
            out.push((color, seqs));
        }
        out
    }

    fn owned(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn single_file_generator_one_color_per_record() {
        let path = write_fasta("one-per-record.fna", &["AAA", "CCC", "GGG"]);
        let mut gen = SeqStreamGeneratorFromSingleFile::new(path, vec![0, 1, 2]);
        let expected = vec![(0, owned(&["AAA"])), (1, owned(&["CCC"])), (2, owned(&["GGG"]))];
        assert_eq!(collect(&mut gen), expected);
    }

    #[test]
    fn single_file_generator_groups_consecutive_records() {
        let path = write_fasta("grouped.fna", &["AAA", "CCC", "GGG", "TTT", "ACG"]);
        let mut gen = SeqStreamGeneratorFromSingleFile::new(path, vec![0, 0, 1, 2, 2]);
        let expected = vec![(0, owned(&["AAA", "CCC"])), (1, owned(&["GGG"])), (2, owned(&["TTT", "ACG"]))];
        assert_eq!(collect(&mut gen), expected);

        // Rewinding must generate the same streams again
        gen.rewind();
        assert_eq!(collect(&mut gen), expected);
    }

    #[test]
    #[should_panic(expected = "not consecutive")]
    fn single_file_generator_rejects_non_consecutive_records() {
        let path = write_fasta("non-consecutive.fna", &["AAA", "CCC", "GGG"]);
        SeqStreamGeneratorFromSingleFile::new(path, vec![0, 1, 0]);
    }

    #[test]
    fn check_consecutive_colors_cases() {
        assert_eq!(check_consecutive_colors(&[]), Ok(()));
        assert_eq!(check_consecutive_colors(&[0, 0, 1, 2, 2, 2]), Ok(()));
        assert_eq!(check_consecutive_colors(&[0, 1, 0]), Err(2)); // Color 0 comes back
        assert_eq!(check_consecutive_colors(&[0, 2, 1]), Err(1)); // Not numbered in file order
        assert_eq!(check_consecutive_colors(&[1]), Err(0));
    }

    #[test]
    fn single_file_generator_all_records_one_color() {
        let path = write_fasta("all-one-color.fna", &["AAA", "CCC", "GGG"]);
        let mut gen = SeqStreamGeneratorFromSingleFile::new(path, vec![0, 0, 0]);
        assert_eq!(collect(&mut gen), vec![(0, owned(&["AAA", "CCC", "GGG"]))]);
    }

    #[test]
    #[should_panic(expected = "not consecutive")]
    fn single_file_generator_rejects_colors_out_of_order() {
        let path = write_fasta("out-of-order.fna", &["AAA", "CCC"]);
        SeqStreamGeneratorFromSingleFile::new(path, vec![1, 0]);
    }
}
