use sbwt;
use bitvec::prelude::*;

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
}