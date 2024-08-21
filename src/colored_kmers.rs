use std::io::Read;

use rand_distr::num_traits::ToBytes;
use sbwt::{self, SbwtIndex, SubsetMatrix, SubsetSeq};
use bitvec::prelude::*;
use simple_sds_sbwt::{self, raw_vector::AccessRaw};

#[derive(Debug)]
pub struct ColoredKmers {
    kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    color_matrix: BitVec, // Concatenated rows
    empty_set: BitVec, // So that we can return a bitslice to an empty set
    n_colors: usize,
}

fn build_sbwt_from_ascii_dump(sbwt_ascii_dump: impl std::io::BufRead){

}

impl ColoredKmers {

    pub fn n_colors(&self) -> usize {
        self.n_colors
    }

    pub fn new_from_new_themisto_index_dump(mut sbwt_ascii_dump: impl std::io::BufRead, themisto_unitig_dump: impl std::io::BufRead, themisto_color_dump: impl std::io::BufRead) {
        let mut buf = String::new();

        let parse_key_value_from_buf = |buf: &mut String| {
            let tokens = buf.split(' ').collect::<Vec<&str>>();
            assert_eq!(tokens.len(), 2);
            (tokens[0].to_owned(), tokens[1].to_owned())
        };

        if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
            let (key, value) = parse_key_value_from_buf(&mut buf);
            assert_eq!(key, "version: ");
            assert_eq!(value, "v0.1");
        } else {
            panic!("Error reading SBWT ascii dump");
        }

        let k: usize = if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
            let (key, value) = parse_key_value_from_buf(&mut buf);
            assert_eq!(key, "k: ");
            value.parse().unwrap()
        } else {
            panic!("Error reading SBWT ascii dump");
        };

        let n_sets: usize = if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
            let (key, value) = parse_key_value_from_buf(&mut buf);
            assert_eq!(key, "number_of_sets: ");
            value.parse().unwrap()
        } else {
            panic!("Error reading SBWT ascii dump");
        };

        let n_kmers: usize = if sbwt_ascii_dump.read_line(&mut buf).unwrap() > 0 {
            let (key, value) = parse_key_value_from_buf(&mut buf);
            assert_eq!(key, "number_of_sets: ");
            value.parse().unwrap()
        } else {
            panic!("Error reading SBWT ascii dump");
        };

        let mut rows: Vec<simple_sds_sbwt::raw_vector::RawVector> = vec![simple_sds_sbwt::raw_vector::RawVector::with_len(n_sets, false); 4];

        // Read from sbwt_ascii_dump byte by byte
        let mut one_byte = [0_u8; 1];
        let mut colex = 0_usize;
        while sbwt_ascii_dump.read(&mut one_byte).unwrap() > 0 {
            let mut c = one_byte[0];
            if c == b'$' {
                colex += 1;
                // Empty set
            } else {
                let end_of_set = c.is_ascii_lowercase();
                c.make_ascii_uppercase();
                match c {
                    b'A' => {
                        rows[0].set_bit(colex, true);
                    },
                    b'C' => {
                        rows[1].set_bit(colex, true);
                    },
                    b'G' => {
                        rows[2].set_bit(colex, true);
                    },
                    b'T' => {
                        rows[3].set_bit(colex, true);
                    },
                    _ => panic!("Invalid character in SBWT ascii dump"),
                }
                if end_of_set {
                    colex += 1;
                }
            }
        }

        assert_eq!(colex, n_sets); // There should be one set for each colex position

        let mut subsetseq = sbwt::SubsetMatrix::new_from_bit_vectors(rows.into_iter().map(simple_sds_sbwt::bit_vector::BitVector::from).collect());
        subsetseq.build_rank();

        todo!();

        
    }

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

    pub fn intersection_pseudoalignment(&self, query: &[u8], minimum_hits: usize) -> BitVec {
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);
        let mut intersection = bitvec![1; self.n_colors]; // Set with all elements (identity element of intersection).
        let mut hit_count = 0_usize;
        for (match_len, colex_range) in index.matching_statistics(query) {
            if match_len == self.kmers.k() {
                hit_count += 1;
                assert_eq!(colex_range.len(), 1);
                intersection &= self.get_color_set_by_row(colex_range.start);
            }

        }
        
        // Return the intersection if there was at least one match of length k
        if hit_count >= minimum_hits {
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
