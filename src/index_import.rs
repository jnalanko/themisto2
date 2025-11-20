use std::io::BufRead;

use bitvec::prelude::*;
use bitvec::vec::BitVec;
use clap::builder::styling::Color;
use sbwt::SubsetSeq;
use simple_sds_sbwt::raw_vector::AccessRaw;

use crate::{coloring_interface::{ColorSetIterStream, ColorSetStream}, iterators::{USizeIterator, USizeIteratorGenerator}};

pub fn ascii_to_int(ascii: &[u8]) -> usize {
    std::str::from_utf8(ascii)
    .expect("Unitig id is not valid utf-8").parse()
    .expect("Could not convert unitig id string to unsigned integer")
}

pub fn get_color_set_id_from_fasta_header(fasta_header: &[u8]) -> usize {
    // The fasta header should look like this " unitig_id=0 color_set_id=0".
    // Note the space at the start.

    let part = fasta_header[1..].split(|c| *c == b' ').nth(1).expect("Color set id missing");
    let mut tokens = part.split(|c| *c == b'=');
    assert_eq!(tokens.next().expect("Color set id missing"), b"color_set_id");
    ascii_to_int(tokens.next().expect("Color set id missing"))
}

pub struct ColorSetDumpIterGenerator<B: BufRead> {
    input: B,
    line: String,
    n_sets_read: usize,
    set_buf: Vec<usize>,
}

pub struct ColorSetDumpSetStream<'a> {
    // Todo: stream over? It gets really complicated with the types involed in std::slice::Split
    set: &'a [usize], 
    pos: usize,
}

impl<'a> USizeIterator<'a> for ColorSetDumpSetStream<'a> {
    fn next(&mut self) -> Option<usize> {
        if self.pos == self.set.len() {
            None
        } else {
            let x = self.set[self.pos];
            self.pos += 1;
            Some(x)
        }
    }
}

fn is_space(c: &u8) -> bool {
    *c == b' '
}

impl<B: BufRead> USizeIteratorGenerator for ColorSetDumpIterGenerator<B> {

    // The type of this iterator gets quite complex because it involves
    // a slice split that is generic over the split predicate.
    type Iter<'a> = ColorSetDumpSetStream<'a> where B: 'a;
    
    fn next<'b>(&'b mut self) -> Option<Self::Iter<'b>> {
        // Lines should look like this:
        // color_set_id=9 size=7 3 4 9 12 14 15 16

        self.set_buf.clear();
        self.line.clear();
        if self.input.read_line(&mut self.line).unwrap() > 0 {

            let line_bytes = self.line.trim_end().as_bytes();
            let mut tokens = line_bytes.split(is_space);

            let first_token = tokens.next().unwrap();
            assert_eq!(&first_token[0..13], b"color_set_id=");
            let color_set_id: usize = ascii_to_int(&first_token[13..]);
            assert!(color_set_id == self.n_sets_read, "Error reading color dump: color set ids are not in order");

            let second_token = tokens.next().unwrap();
            assert_eq!(&second_token[0..5], b"size=");
            let _ = ascii_to_int(&second_token[5..]); // Length of the color set

            self.set_buf.extend(tokens.map(ascii_to_int));
            
            self.n_sets_read += 1;

            Some(Self::Iter{set: &self.set_buf, pos: 0})
        } else {
            None
        }
    }
}

impl<B: BufRead> ColorSetDumpIterGenerator<B> {
    pub fn new(input: B) -> Self {
        Self { input, line: String::new(), n_sets_read: 0, set_buf: vec![]}
    }
}

// Parses the set and pushed it to `push_set_to`. Returns the color set id.
pub fn parse_color_set_line(line: &String, push_set_to: &mut Vec<usize>) -> usize {
    // Lines should look like this:
    // color_set_id=9 size=7 3 4 9 12 14 15 16

    let line_bytes = line.trim_end().as_bytes();
    let mut tokens = line_bytes.split(|c| *c == b' ');

    let first_token = tokens.next().unwrap();
    assert_eq!(&first_token[0..13], b"color_set_id=");
    let color_set_id: usize = ascii_to_int(&first_token[13..]);

    let second_token = tokens.next().unwrap();
    assert_eq!(&second_token[0..5], b"size=");
    let _ = ascii_to_int(&second_token[5..]); // Length of the color set

    push_set_to.extend(tokens.map(ascii_to_int));

    color_set_id
}

pub fn sbwt_ascii_dump_to_sbwt_index(mut sbwt_ascii_dump: impl std::io::BufRead, precalc_prefix_length: usize) -> sbwt::SbwtIndex<sbwt::SubsetMatrix> {

    let mut buf = String::new();

    let parse_key_value_from_buf = |buf: &mut String| {
        let tokens = buf.split(' ').collect::<Vec<&str>>();
        assert_eq!(tokens.len(), 2);
        (tokens[0].to_owned(), tokens[1].strip_suffix('\n').unwrap().to_owned())
    };

    if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
        let (key, value) = parse_key_value_from_buf(&mut buf);
        assert_eq!(key, "version:");
        assert_eq!(value, "v0.1");
    } else {
        panic!("Error reading SBWT ascii dump");
    }
    buf.clear();

    let k: usize = if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
        let (key, value) = parse_key_value_from_buf(&mut buf);
        assert_eq!(key, "k:");
        value.parse().unwrap()
    } else {
        panic!("Error reading SBWT ascii dump");
    };
    buf.clear();

    let n_sets: usize = if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
        let (key, value) = parse_key_value_from_buf(&mut buf);
        assert_eq!(key, "number_of_sets:");
        value.parse().unwrap()
    } else {
        panic!("Error reading SBWT ascii dump");
    };
    buf.clear();

    let n_kmers: usize = if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
        let (key, value) = parse_key_value_from_buf(&mut buf);
        assert_eq!(key, "number_of_kmers:");
        value.parse().unwrap()
    } else {
        panic!("Error reading SBWT ascii dump");
    };
    buf.clear();

    let mut rows: Vec<BitVec<u64, Lsb0>> = vec![bitvec![u64, Lsb0; n_sets; 0]; 4];

    // Read from sbwt_ascii_dump byte by byte
    let mut one_byte = [0_u8; 1];
    let mut colex = 0_usize;
    while sbwt_ascii_dump.read(&mut one_byte).unwrap() > 0 {
        let mut c = one_byte[0];
        if c == b'\n' {
            break
        } else if c == b'$' {
            colex += 1;
            // Empty set
        } else {
            let end_of_set = c.is_ascii_lowercase();
            c.make_ascii_uppercase();
            match c {
                b'A' => {
                    rows[0].set(colex, true);
                },
                b'C' => {
                    rows[1].set(colex, true);
                },
                b'G' => {
                    rows[2].set(colex, true);
                },
                b'T' => {
                    rows[3].set(colex, true);
                },
                _ => panic!("Invalid character in SBWT ascii dump"),
            }
            if end_of_set {
                colex += 1;
            }
        }
    }

    assert_eq!(colex, n_sets); // There should be one set for each colex position

    let mut subsetseq = sbwt::SubsetMatrix::new_from_bit_vectors(rows);
    subsetseq.build_rank();

    sbwt::SbwtIndex::from_subset_seq(subsetseq, n_kmers, k, precalc_prefix_length)
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct IndexDumpMetadata {
    pub num_unitigs: usize,
    pub num_colors: usize,
    pub num_color_sets: usize,
    pub k: usize
}

pub fn read_index_dump_metadata(mut reader: impl BufRead) -> IndexDumpMetadata {

    // File should look like this:
    // num_colors=3682
    // num_unitigs=9314735
    // num_color_sets=5591009
    // k=31

    let mut line = String::new();

    let mut num_unitigs = None;
    let mut num_colors = None;
    let mut num_color_sets = None;
    let mut k = None;

    while reader.read_line(&mut line).unwrap() > 0 {
        let line_bytes = line.trim_end().as_bytes();
        let mut tokens = line_bytes.split(|c| *c == b'=');

        let first_token = tokens.next().unwrap();
        let second_token = tokens.next().unwrap();

        match first_token {
            b"num_colors" => num_colors = Some(ascii_to_int(second_token)),
            b"num_unitigs" => num_unitigs = Some(ascii_to_int(second_token)),
            b"num_color_sets" => num_color_sets = Some(ascii_to_int(second_token)),
            b"k" => k = Some(ascii_to_int(second_token)),
            _ => panic!("Unknown metadata field: {}", line)
        }

        line.clear();
    }

    IndexDumpMetadata {
        num_unitigs: num_unitigs.expect("num_unitigs missing from metadata"),
        num_colors: num_colors.expect("num_colors missing from metadata"),
        num_color_sets: num_color_sets.expect("num_color_sets missing from metadata"),
        k: k.expect("k missing from metadata")
    }

}

pub fn build_colex_to_color_set_mapping(unitig_input: impl BufRead + Send + 'static, sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray) -> Vec<usize> {
    let mut reader = jseqio::reader::DynamicFastXReader::new(unitig_input).unwrap();
    let index = sbwt::StreamingIndex::new(sbwt, lcs);
    let mut colex_to_color_set_id = vec![0_usize; sbwt.n_sets()];
    while let Some(rec) = reader.read_next().unwrap() {
        let color_set_id = get_color_set_id_from_fasta_header(rec.head);
        for (match_len, colex_range) in index.matching_statistics(rec.seq) {
            if match_len == sbwt.k() {
                assert_eq!(colex_range.len(), 1);
                colex_to_color_set_id[colex_range.start] = color_set_id;
            }
        }
    }
    colex_to_color_set_id
}
