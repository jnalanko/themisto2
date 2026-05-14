use std::io::BufRead;

use crate::iterators::{USizeIterator, USizeIteratorGenerator};

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
