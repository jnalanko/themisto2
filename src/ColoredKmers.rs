use sbwt;
use bitvec::prelude::*;

#[derive(Debug)]
struct ColoredKmers {
    kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    color_matrix: Vec<BitVec>,
}

impl ColoredKmers {
    fn new_from_themisto_color_dump(kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>, mut themisto_color_dump_stream: impl std::io::BufRead) -> Self {
        let mut color_matrix: Vec<BitVec> = vec![bitvec![]; kmers.n_sets()];
        let mut line = String::new();

        // Todo: this operates with unicode strings even though the input is just ASCII and there is no reason to ever
        // have unicode characters in it.
        while themisto_color_dump_stream.read_line(&mut line).unwrap() > 0 {
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

            color_matrix[row] = bitvec;
            
            line.clear();
        }

        ColoredKmers{kmers, color_matrix}
    }

    fn get_color_set(&self, kmer: &[u8]) -> Option<&BitVec> {
        let row = self.kmers.search(kmer)?.start;        
        Some(&self.color_matrix[row])
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

        let (sbwt, _) = sbwt::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(kmers_slices.first().unwrap().len()).run_from_slices(&kmers_slices);
        let colored_kmers = ColoredKmers::new_from_themisto_color_dump(sbwt, dump.as_bytes());

        for (i, kmer) in kmers_slices.iter().enumerate() {
            let color_set = colored_kmers.get_color_set(kmer).unwrap();
            eprintln!("{}, {:?}", String::from_utf8(kmer.to_vec()).unwrap(), color_set);
            let mut color_set_string = String::new();
            for b in color_set {
                color_set_string.push(match *b {true => '1', false => '0'});
            }
            assert_eq!(color_set_string, bitvec_strings[i]);
        }


    }
}
