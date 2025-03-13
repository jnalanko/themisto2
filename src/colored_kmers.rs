use std::{io::BufRead, path::{Path, PathBuf}};

use clap::builder::styling::Color;
use sbwt::{self, BitPackedKmerSorting, SbwtIndex, SeqStream, SubsetMatrix, SubsetSeq};
use bitvec::prelude::*;
use simple_sds_sbwt::{self, raw_vector::AccessRaw};

#[derive(Debug)]
pub struct ColoredKmers {
    kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    distinct_color_sets: BitVec, // Concatenation of distinct color sets
    colex_to_color_set_id: Vec<usize>,
    empty_set: BitVec, // So that we can return a bitslice to an empty set
    n_colors: usize,
}

fn ascii_to_int(ascii: &[u8]) -> usize {
    std::str::from_utf8(ascii)
    .expect("Unitig id is not valid utf-8").parse()
    .expect("Could not convert unitig id string to unsigned integer")
}

fn get_color_set_id_from_fasta_header(fasta_header: &[u8]) -> usize {
    // The fasta header should look like this " unitig_id=0 color_set_id=0".
    // Note the space at the start.

    let part = fasta_header[1..].split(|c| *c == b' ').nth(1).expect("Color set id missing");
    let mut tokens = part.split(|c| *c == b'=');
    assert_eq!(tokens.next().expect("Color set id missing"), b"color_set_id");
    ascii_to_int(tokens.next().expect("Color set id missing"))
}

// Returns the concatenation of distinct color sets
fn read_color_sets(mut reader: impl BufRead, num_color_sets: usize, num_colors: usize) -> BitVec {

    // Lines should look like this:
    // color_set_id=9 size=7 3 4 9 12 14 15 16

    let mut color_sets = bitvec![0; num_color_sets * num_colors]; // Concatenation of distinct color sets

    let mut line = String::new();

    let bar = indicatif::ProgressBar::new(num_color_sets as u64);
    while reader.read_line(&mut line).unwrap() > 0 {
        let line_bytes = line.trim_end().as_bytes();
        let mut tokens = line_bytes.split(|c| *c == b' ');

        let first_token = tokens.next().unwrap();
        assert_eq!(&first_token[0..13], b"color_set_id=");
        let color_set_id: usize = ascii_to_int(&first_token[13..]);
        assert!(color_set_id < num_color_sets);

        let second_token = tokens.next().unwrap();
        assert_eq!(&second_token[0..5], b"size=");
        let _ = ascii_to_int(&second_token[5..]); // Length of the color set

        for color in tokens.map(ascii_to_int) {
            color_sets.set(color_set_id*num_colors + color, true);
        }

        line.clear();
        bar.inc(1);
    }
    bar.finish();

    color_sets
}

fn sbwt_ascii_dump_to_sbwt_index(mut sbwt_ascii_dump: impl std::io::BufRead, precalc_prefix_length: usize) -> sbwt::SbwtIndex<sbwt::SubsetMatrix> {

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

    let mut rows: Vec<simple_sds_sbwt::raw_vector::RawVector> = vec![simple_sds_sbwt::raw_vector::RawVector::with_len(n_sets, false); 4];

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

    sbwt::SbwtIndex::from_subset_seq(subsetseq, n_kmers, k, precalc_prefix_length)
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
struct ThemistoDumpMetadata {
    num_unitigs: usize,
    num_colors: usize,
    num_color_sets: usize,
    k: usize
}

fn read_themisto_dump_metadata(mut reader: impl BufRead) -> ThemistoDumpMetadata {

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

    ThemistoDumpMetadata {
        num_unitigs: num_unitigs.expect("num_unitigs missing from metadata"),
        num_colors: num_colors.expect("num_colors missing from metadata"),
        num_color_sets: num_color_sets.expect("num_color_sets missing from metadata"),
        k: k.expect("k missing from metadata")
    }

}

fn build_colex_to_color_set_mapping(unitig_input: impl BufRead + Send + 'static, sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray) -> Vec<usize> {
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

#[derive(serde::Serialize)]
pub struct PseudoalignmentData {
    pub hit_counts: Vec<usize>,
    pub distinguishing_hit_counts: Vec<usize>,
    pub n_relevant_kmers: usize,
    pub n_all_kmers: usize,
}

impl PseudoalignmentData {
    pub fn new_empty(n_colors: usize) -> Self {
        Self {
            hit_counts: vec![0; n_colors],
            distinguishing_hit_counts: vec![0; n_colors],
            n_relevant_kmers: 0,
            n_all_kmers: 0,
        }
    }
}

impl ColoredKmers {

    pub fn n_colors(&self) -> usize {
        self.n_colors
    }

    pub fn n_kmers(&self) -> usize {
        self.kmers.n_kmers()
    }

    pub fn get_k(&self) -> usize {
        self.kmers.k()
    }

    pub fn new_from_new_themisto_index_dump(sbwt_ascii_dump: impl std::io::BufRead, themisto_metadata_dump: impl std::io::BufRead, themisto_unitig_dump: impl std::io::BufRead + Send + 'static, themisto_color_dump: impl std::io::BufRead, precalc_prefix_length: usize) -> Self {
        log::info!("Reading metadata");
        let metadata = read_themisto_dump_metadata(themisto_metadata_dump);
        log::info!("Reading SBWT dump");
        let sbwt_index = sbwt_ascii_dump_to_sbwt_index(sbwt_ascii_dump, precalc_prefix_length);
        log::info!("Reading Distinct color sets");
        let distinct_color_sets = read_color_sets(themisto_color_dump, metadata.num_color_sets, metadata.num_colors); 
        log::info!("Reading Building LCS array");
        let lcs = sbwt::LcsArray::from_sbwt(&sbwt_index);
        log::info!("Building colex to color set id mapping");
        let colex_to_color_set_id = build_colex_to_color_set_mapping(themisto_unitig_dump, &sbwt_index, &lcs);
        Self { kmers: sbwt_index, lcs, distinct_color_sets, colex_to_color_set_id, empty_set: bitvec![0; metadata.num_colors], n_colors: metadata.num_colors}
    }


    pub fn get_color_set(&self, kmer: &[u8]) -> &BitSlice {
        match self.kmers.search(kmer){
            Some(range) => {
                let row = range.start;
                &self.distinct_color_sets[row*self.n_colors..(row+1)*self.n_colors]
            }
            None => &self.empty_set,
        }
    }

    fn get_color_set_by_id(&self, id: usize) -> &BitSlice {
        &self.distinct_color_sets[id*self.n_colors..(id+1)*self.n_colors]
    }

    pub fn serialize<W: std::io::Write>(&self, mut out: &mut W) {
        self.kmers.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();

        out.write_all(&(self.n_colors as u64).to_le_bytes()).unwrap();
        out.write_all(&(self.distinct_color_sets.len() as u64).to_le_bytes()).unwrap();

        bincode::serialize_into(&mut out, &self.distinct_color_sets).unwrap();
        bincode::serialize_into(&mut out, &self.colex_to_color_set_id).unwrap();
    }

    pub fn load<R: std::io::Read>(mut input: &mut R) -> Self {
        let kmers = SbwtIndex::<SubsetMatrix>::load(input).unwrap();
        let lcs = sbwt::LcsArray::load(input).unwrap();

        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf).unwrap();
        let n_colors = u64::from_le_bytes(buf);

        input.read_exact(&mut buf).unwrap();
        let _ = u64::from_le_bytes(buf); // Total length of distinct color sets

        let distinct_color_sets: BitVec = bincode::deserialize_from(&mut input).unwrap();
        let colex_to_color_set_id: Vec<usize> = bincode::deserialize_from(&mut input).unwrap(); // Todo: u64 instead of usize

        ColoredKmers{kmers, lcs, colex_to_color_set_id, n_colors: n_colors as usize, distinct_color_sets, empty_set: bitvec![0; n_colors as usize]}
    }

    pub fn intersection_pseudoalignment(&self, query: &[u8], minimum_hits: usize) -> BitVec {
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);
        let mut intersection = bitvec![1; self.n_colors]; // Set with all elements (identity element of intersection).
        let mut hit_count = 0_usize;
        for (match_len, colex_range) in index.matching_statistics(query) {
            if match_len == self.kmers.k() {
                hit_count += 1;
                assert_eq!(colex_range.len(), 1);
                let id = self.colex_to_color_set_id[colex_range.start];
                intersection &= self.get_color_set_by_id(id);
            }
        }
        
        // Return the intersection if there was at least one match of length k
        if hit_count >= minimum_hits {
            intersection
        } else {
            self.empty_set.clone()
        }
    }

    // All counts smaller than the compatibility threshold are set to zero.
    // Distinguishing k-mers are determined among colors whose hits exceed the compatibility threshold. 
    pub fn compute_pseudoalignment_data(&self, query: &[u8], compatibility_threshold: usize) -> PseudoalignmentData {
        if query.len() < self.kmers.k() {
            return PseudoalignmentData::new_empty(self.n_colors);
        }

        let mut hit_counts: Vec<usize> = vec![0; self.n_colors];
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);
        let mut color_sets: Vec<Option<&BitSlice>> = vec![None; query.len() - self.kmers.k() + 1];
        let mut n_relevant_kmers = 0_usize;

        // Retrieve color sets
        for (i, (match_len, colex_range)) in index.matching_statistics(query).iter().enumerate() {
            if *match_len == self.kmers.k() {
                assert_eq!(colex_range.len(), 1);
                let id = self.colex_to_color_set_id[colex_range.start];
                color_sets[i + 1 - self.kmers.k()] = Some(self.get_color_set_by_id(id));
                n_relevant_kmers += 1;
            }
        }

        if n_relevant_kmers == 0 {
            return PseudoalignmentData::new_empty(self.n_colors); 
        }

        // Count hits (flatten removes nones).
        color_sets.iter().flatten().for_each(|bitmap| {
            for (color, bit) in bitmap.iter().enumerate() {
                if *bit {
                    hit_counts[color] += 1;
                }
            }
        });

        // Set all counts smaller than the compatibility threshold to zero.
        hit_counts.iter_mut().for_each(|n| *n *= (*n >= compatibility_threshold) as usize);

        let mut candidate_set = bitvec!(0; self.n_colors); // All colors with at least one hit
        for (color, count) in hit_counts.iter().enumerate() {
            if *count >= compatibility_threshold {
                candidate_set.set(color, true);
            }
        }

        let is_distinguishing_color_set = |x: &BitSlice| {
            let y = x.to_bitvec() & (&candidate_set); // Set intersection
            y.count_ones() > 0 && (y & (&candidate_set)) != candidate_set // Proper non-empty subset of the candidate set
        };

        let mut distinguishing_hit_counts: Vec<usize> = vec![0; self.n_colors];
        for color_set in color_sets.iter().flatten().filter(|x| is_distinguishing_color_set(x)) {
            for (i, x) in color_set.iter().enumerate() {
                if hit_counts[i] >= compatibility_threshold && *x {
                    distinguishing_hit_counts[i] += 1;
                }
            }
        }

        let n_all_kmers = query.len() as i64 - self.get_k() as i64 + 1;
        let n_all_kmers = std::cmp::max(0, n_all_kmers) as usize; 
        PseudoalignmentData{hit_counts, distinguishing_hit_counts, n_relevant_kmers, n_all_kmers}
    }

    // Returns pairs (color id, score). Not all colors are necessarily present in the output.
    pub fn compute_distinguishing_scores(&self, query: &[u8]) -> Vec<(usize, f64)> {
        if query.len() < self.kmers.k() {
            return vec![];
        }

        //let mut hit_counts: Vec<usize> = vec![0; self.n_colors];
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);
        let mut color_sets: Vec<Option<&BitSlice>> = vec![None; query.len() - self.kmers.k() + 1];
        let mut n_relevant_kmers = 0_usize;

        // Retrieve color sets
        for (i, (match_len, colex_range)) in index.matching_statistics(query).iter().enumerate() {
            if *match_len == self.kmers.k() {
                assert_eq!(colex_range.len(), 1);
                let id = self.colex_to_color_set_id[colex_range.start];
                color_sets[i + 1 - self.kmers.k()] = Some(self.get_color_set_by_id(id));
                n_relevant_kmers += 1;
            }
        }

        if n_relevant_kmers == 0 {
            return vec![];
        }

        let mut union = bitvec!(0; self.n_colors); // All colors with at least one hit

        // Take the union of the found color sets. flatten() filters out None values.
        color_sets.iter().flatten().for_each(|x| { union |= *x; });

        let mut distinguishing_hit_counts: Vec<usize> = vec![0; self.n_colors];
        for color_set in color_sets.iter().flatten().filter(|&&x| x != union) {
            for (i, x) in color_set.iter().enumerate() {
                if *x {
                    distinguishing_hit_counts[i] += 1;
                }
            }
        }

        let max_distinguishing_hits = distinguishing_hit_counts.iter().max().unwrap();
        distinguishing_hit_counts.iter().enumerate().filter(|(_, hits)| **hits > 0).map(|(color, hits)| (color, *hits as f64 / *max_distinguishing_hits as f64)).collect()

    }
}

struct InputStream {
    dbs: Vec<jseqio::seq_db::SeqDB>,
    cur_db_idx: usize, // Index of the db currently being iterated over
    seq_idx_in_cur_db: usize,
}

impl InputStream {
    fn new(filenames: &[&Path]) -> InputStream {
        let mut dbs: Vec<jseqio::seq_db::SeqDB> = vec![];
        for path in filenames {
            let reader = jseqio::reader::DynamicFastXReader::from_file(path).unwrap();
            let (mut fw, rc) = reader.into_db_with_revcomp().unwrap();

            if fw.sequence_count() == 0 {
                panic!("No sequences found in file {}", path.display());
            }

            // Append reverse complement records to the forward database
            for rec in rc.iter() {
                fw.push_record(rec);
            }
            dbs.push(fw);
        }
        Self {dbs, cur_db_idx: 0, seq_idx_in_cur_db: 0}
    }

    fn reset(&mut self) {
        self.cur_db_idx = 0;
        self.seq_idx_in_cur_db = 0;
    }
}

impl SeqStream for InputStream {

    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.cur_db_idx == self.dbs.len() {
            return None; // Done
        }

        // Fetch the next sequence
        let db = &self.dbs[self.cur_db_idx];
        assert!(db.sequence_count() > 0);
        let seq = db.get(self.seq_idx_in_cur_db).seq;

        // Update the "cursor"
        self.seq_idx_in_cur_db += 1;
        if self.seq_idx_in_cur_db == db.sequence_count() {
            self.cur_db_idx += 1;
            self.seq_idx_in_cur_db = 0;
        }
        Some(seq)
    } 
}

impl ColoredKmers {
    pub fn new(filenames: &[&Path], k: usize, n_threads: usize, temp_dir: &Path) -> Self {
        let input_stream = InputStream::new(filenames);
        let sbwt_builder = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false) // Already added in the input stream
            .k(k)
            .build_lcs(true)
            .n_threads(n_threads)
            .precalc_length(8)
            .algorithm(BitPackedKmerSorting::new()
                .dedup_batches(true)
                .temp_dir(temp_dir)
        );
        let (sbwt, lcs) = sbwt_builder.run(input_stream);
        todo!();
    } 
}

#[cfg(test)]
mod tests {
    use sbwt::BitPackedKmerSorting;

    use super::*;

    /*
        kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
        lcs: sbwt::LcsArray,
        distinct_color_sets: BitVec, // Concatenation of distinct color sets
        colex_to_color_set_id: Vec<usize>,
        empty_set: BitVec, // So that we can return a bitslice to an empty set
        n_colors: usize,
    */

    #[test]
    fn test_distinguishing_scores(){
        let red_seq: &[u8] = b"AAACATCGATCGTACGTACGTCAGCTACTGCA";
        let blue_seq: &[u8] = b"CACTCTATCGCGTTATCTTACGATCATGCTAGC";
        let green_seq: &[u8] = b"ACATCGGCGTATCTATCTACGATCGTACGTCA";
        let uncolored_seq: &[u8] = b"GGATTCGGATCTATCGTAGCTGTACGTGCTGAC";
        let red_and_blue_seq: &[u8] = b"TTAGCTATCGTATCCGATCACGTACGTAGTCAA";
        let red_and_blue_and_green_seq: &[u8] = b"CCGTTATCGGCCTATACTATCGACTACGTAGC";
        let k = 12; // Let's hope I didn't accidentally repeat any kmers

        let all_seqs: &[&[u8]] = &[&red_seq, &blue_seq, &green_seq, &uncolored_seq, &red_and_blue_seq, &red_and_blue_and_green_seq];
        let distinct_color_sets = bitvec![0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,1,1];
        let seq_color_set_ids: Vec<usize> = vec![1, 2, 3, 0, 4, 5];

        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(k).build_lcs(true).run_from_slices(all_seqs);
        let lcs = lcs.unwrap();

        // Build colex to color id mapping
        let mut colex_to_color_set_id = vec![0_usize; sbwt.n_sets()];
        for (seq_id, seq) in all_seqs.iter().enumerate() {
            for (match_len, colex_range) in sbwt::StreamingIndex::new(&sbwt, &lcs).matching_statistics(seq).iter() {
                if *match_len == k {
                    assert_eq!(colex_range.len(), 1);
                    colex_to_color_set_id[colex_range.start] = seq_color_set_ids[seq_id];
                }
            }
        }

        let index = ColoredKmers{kmers: sbwt, lcs, distinct_color_sets, colex_to_color_set_id, empty_set: bitvec![0,0,0], n_colors: 3};

        // Concatenate the sequences
        let all_seqs_concatenated: Vec<u8> = all_seqs.iter().flat_map(|x| x.iter().copied()).collect();

        //let total_kmers = all_seqs_concatenated.len() - k + 1;
        //let n_distinguishing = red_seq.len() - k + 1 + blue_seq.len() - k + 1 + green_seq.len() - k + 1 + red_and_blue_seq.len() - k + 1;
        let n_red_distinguishing_hits = red_seq.len() - k + 1 + red_and_blue_seq.len() - k + 1;
        let n_blue_distinguishing_hits = blue_seq.len() - k + 1 + red_and_blue_seq.len() - k + 1;
        let n_green_distinguishing_hits = green_seq.len() - k + 1;
        let max_distinguishing_hits = n_red_distinguishing_hits.max(n_blue_distinguishing_hits).max(n_green_distinguishing_hits);
        let true_scores = [n_red_distinguishing_hits as f64 / max_distinguishing_hits as f64, n_blue_distinguishing_hits as f64 / max_distinguishing_hits as f64, n_green_distinguishing_hits as f64 / max_distinguishing_hits as f64];

        let scores = index.compute_distinguishing_scores(&all_seqs_concatenated);

        eprintln!("{:?}", scores);
        eprintln!("{:?}", true_scores);
        
        let epsilon = 1e-6;
        assert_eq!(true_scores.len(), scores.len());
        for (color, score) in scores {
            assert!((score - true_scores[color]).abs() < epsilon);
        }

    }

    /*
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
*/

    #[test]
    fn bitvec_serialization() {
        // Just checking how much overhead there is
        let bv = bitvec![0; 0];
        let bytes = bincode::serialize(&bv).unwrap(); 
        eprintln!("{}", bytes.len());
    }
}
