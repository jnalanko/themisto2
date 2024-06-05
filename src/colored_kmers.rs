use std::io::Read;

use rand_distr::num_traits::ToBytes;
use sbwt::{self, SbwtIndex, SubsetMatrix};
use bitvec::prelude::*;

#[derive(Debug)]
pub struct ColoredKmers {
    kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    color_matrix: BitVec, // Concatenated rows
    empty_set: BitVec, // So that we can return a bitslice to an empty set
    n_colors: usize,
}

impl ColoredKmers {
    pub fn new_from_themisto_color_dump(kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, mut themisto_color_dump_stream: impl std::io::BufRead, n_colors: usize) -> Self {
        let mut color_matrix: BitVec = bitvec![0; kmers.n_sets() * n_colors]; // Concatenated rows
        let mut line = String::new();

        // Todo: this operates with unicode strings even though the input is just ASCII and there is no reason to ever
        // have unicode characters in it.
        let bar = indicatif::ProgressBar::new(kmers.n_sets() as u64);
        bar.set_style(indicatif::ProgressStyle::with_template("[ETA: {eta_precise}] {wide_bar} {pos:>7}/{len:7} {msg}")
        .unwrap());

        while themisto_color_dump_stream.read_line(&mut line).unwrap() > 0 {
            bar.inc(1);
            assert!(line.ends_with('\n'));
            line.pop(); // Trim newline at the end

            let line_bytes = line.as_bytes();
            let mut tokens = line_bytes.split(|c| *c == b' ');
            
            let kmer = tokens.next().unwrap();
            let bits = tokens.next().unwrap();
            assert!(tokens.next() == None);
            assert_eq!(kmer.len(), kmers.k());

            let mut bitvec = bitvec![0; bits.len()];
            for (i, &b) in bits.iter().enumerate() {
                assert!(b == b'0' || b == b'1');
                bitvec.set(i, b == b'1');
            }

            let row = match kmers.search(kmer) {
                Some(range) => range.start,
                None => panic!("Kmer {} not found in sbwt", String::from_utf8(kmer.to_vec()).unwrap()),
            };

            color_matrix[row*n_colors..(row+1)*n_colors].copy_from_bitslice(&bitvec);
            
            line.clear();
        }
        bar.finish();

        ColoredKmers{kmers, lcs, color_matrix, n_colors, empty_set: bitvec![0; n_colors]}
    }

    pub fn get_color_set(&self, kmer: &[u8]) -> &BitSlice {
        match self.kmers.search(kmer){
            Some(range) => {
                let row = range.start;
                &self.color_matrix[row*self.n_colors..(row+1)*self.n_colors]
            }
            None => &self.empty_set,
        }
    }

    fn get_color_set_by_row(&self, row: usize) -> &BitSlice {
        &self.color_matrix[row*self.n_colors..(row+1)*self.n_colors]
    }

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) {
        self.kmers.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();

        out.write_all(&(self.n_colors as u64).to_le_bytes()).unwrap();
        bincode::serialize_into(out, &self.color_matrix).unwrap();
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> Self {
        let kmers = SbwtIndex::<SubsetMatrix>::load(input).unwrap();
        let lcs = sbwt::LcsArray::load(input).unwrap();

        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf).unwrap();
        let n_colors = u64::from_le_bytes(buf);

        let color_matrix: BitVec = bincode::deserialize_from(input).unwrap();

        ColoredKmers{kmers, lcs, n_colors: n_colors as usize, color_matrix, empty_set: bitvec![0; n_colors as usize]}
    }

    pub fn intersection_pseudoalignment(self, query: &[u8]) -> BitVec {
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);
        let mut intersection = bitvec![1; self.n_colors]; // Set with all elements (identity element of intersection).
        let mut at_least_one = false;
        for (match_len, colex_range) in index.matching_statistics(query) {
            if match_len == self.kmers.k() {
                at_least_one = true;
                assert_eq!(colex_range.len(), 1);
                intersection &= self.get_color_set_by_row(colex_range.start);
            }

        }
        
        // Return the intersection if there was at least one match of length k
        if at_least_one {
            intersection
        } else {
            self.empty_set.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use sbwt::BitPackedKmerSorting;

    use super::*;

    #[test]
    fn from_themisto_color_dump(){
        let dump = 
"\
AGATTAGAGTGTCTTTTTCTTTTGCGAGTAG 0000000001001101010000000000000000000000000000000000000000000000000
AGATTAGGGTGTCTTTTTCTTTTGCGAGTAG 0000000011111011101010000000000000000000000000000000000000000000000
GGATTAGGGTGTCTTTTTCTTTTGCGAGTAG 0000000000000001000000000000000000000000000000000000000000000000000
GTACATATCCAGCGCCGCGTTTTGCGAGTAG 0000000000000000000000100000000000000000000000000000000000000000000
GTACATGTCCAGCGCCGCGTTTTGCGAGTAG 0000000000000000000000000000000000000011000000001000000000000000100
ATACATATCCAGCGGCGCGTTTTGCGAGTAG 0000000000000000000000000001111111111111111111111111111111111111111
GAGTAAACAACCTCTGACTTTTTGCGAGTAG 0000000000000000000000000000000000000000000000001000010000000000000
TATATCTTTTTCATACGCTTTTTGCGAGTAG 0000000100000000000000000000000000000000000000000000000000000000000
TCAGTTTTTTACCATGGCTTTTTGCGAGTAG 1000000000000000000000000000000000000000000000000000000000000000000
";

        eprintln!("{}", dump);

        let bitvec_strings = ["0000000001001101010000000000000000000000000000000000000000000000000", "0000000011111011101010000000000000000000000000000000000000000000000", "0000000000000001000000000000000000000000000000000000000000000000000", "0000000000000000000000100000000000000000000000000000000000000000000", "0000000000000000000000000000000000000011000000001000000000000000100", "0000000000000000000000000001111111111111111111111111111111111111111", "0000000000000000000000000000000000000000000000001000010000000000000", "0000000100000000000000000000000000000000000000000000000000000000000", "1000000000000000000000000000000000000000000000000000000000000000000"];

        let kmers_data = [b"AGATTAGAGTGTCTTTTTCTTTTGCGAGTAG", b"AGATTAGGGTGTCTTTTTCTTTTGCGAGTAG", b"GGATTAGGGTGTCTTTTTCTTTTGCGAGTAG", b"GTACATATCCAGCGCCGCGTTTTGCGAGTAG", b"GTACATGTCCAGCGCCGCGTTTTGCGAGTAG", b"ATACATATCCAGCGGCGCGTTTTGCGAGTAG", b"GAGTAAACAACCTCTGACTTTTTGCGAGTAG", b"TATATCTTTTTCATACGCTTTTTGCGAGTAG", b"TCAGTTTTTTACCATGGCTTTTTGCGAGTAG"];
        let kmers_slices = kmers_data.map(|x| x.as_slice());

        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(kmers_slices.first().unwrap().len()).build_lcs(true).run_from_slices(&kmers_slices);

        let serialized_bytes = { // Build index and serialize to also test serialization
            let colored_kmers = ColoredKmers::new_from_themisto_color_dump(sbwt, lcs.unwrap(), dump.as_bytes(), bitvec_strings.first().unwrap().len());
            let mut serialized_bytes = Vec::<u8>::new();
            colored_kmers.serialize(&mut std::io::Cursor::new(&mut serialized_bytes));
            serialized_bytes
        };

        let colored_kmers = ColoredKmers::load(&mut std::io::Cursor::new(serialized_bytes)); // Load back

        for (i, kmer) in kmers_slices.iter().enumerate() {
            let color_set = colored_kmers.get_color_set(kmer);
            eprintln!("{}, {:?}", String::from_utf8(kmer.to_vec()).unwrap(), color_set);
            let mut color_set_string = String::new();
            for b in color_set {
                color_set_string.push(match *b {true => '1', false => '0'});
            }
            assert_eq!(color_set_string, bitvec_strings[i]);
        }
    }

    #[test]
    fn bitvec_serialization() {
        // Just checking how much overhead there is
        let bv = bitvec![0; 0];
        let bytes = bincode::serialize(&bv).unwrap(); 
        eprintln!("{}", bytes.len());
    }
}
