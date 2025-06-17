use std::{ops::DerefMut, path::Path, sync::{Arc, Mutex}};

use crossbeam::channel::{Receiver, RecvError, Sender};
use sbwt::{self, BitPackedKmerSorting, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use bitvec::prelude::*;

use crate::{coloring::CompactColexColoring, themisto1_compatibility::{build_colex_to_color_set_mapping, read_color_sets, read_themisto_dump_metadata, sbwt_ascii_dump_to_sbwt_index}};

#[derive(Debug, Clone)]
pub struct ColoredKmers {
    kmers: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    distinct_color_sets: BitVec, // Concatenation of distinct color sets
    colex_to_color_set_id: Vec<usize>,
    empty_set: BitVec, // So that we can return a bitslice to an empty set
    n_colors: usize,
}

#[derive(serde::Serialize)]
pub struct PseudoalignmentData {
    pub hit_counts: Vec<usize>,
    pub distinguishing_hit_counts: Vec<usize>,
    pub unique_hit_counts: Vec<usize>,
    pub n_unique_kmers: usize,
    pub n_relevant_kmers: usize,
    pub n_all_kmers: usize,
}

impl PseudoalignmentData {
    pub fn new_empty(n_colors: usize) -> Self {
        Self {
            hit_counts: vec![0; n_colors],
            distinguishing_hit_counts: vec![0; n_colors],
            unique_hit_counts: vec![0; n_colors],
            n_unique_kmers: 0,
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
        let lcs = sbwt::LcsArray::from_sbwt(&sbwt_index, 1);
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

    pub fn get_all_color_sets(&self, seq: &[u8]) -> Vec<&BitSlice> {
        let mut sets = Vec::new();
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);

        // Start iterating from index k-1 because those before do not correspond to a full k-mer
        for (match_len, colex_range) in index.matching_statistics(seq)[self.kmers.k()-1 ..].iter() {
            if *match_len == self.kmers.k() {
                assert_eq!(colex_range.len(), 1);
                let id = self.colex_to_color_set_id[colex_range.start];
                sets.push(self.get_color_set_by_id(id));
            } else {
                sets.push(&self.empty_set);
            }
        }
        sets
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
        let mut unique_hit_counts: Vec<usize> = vec![0; self.n_colors];
        let index = sbwt::StreamingIndex::new(&self.kmers, &self.lcs);
        let mut color_sets: Vec<Option<&BitSlice>> = vec![None; query.len() - self.kmers.k() + 1];
        let mut n_relevant_kmers = 0_usize;
        let mut n_unique_kmers = 0_usize;

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
            let is_singleton = bitmap.count_ones() == 1;
            n_unique_kmers += is_singleton as usize;
            for (color, bit) in bitmap.iter().enumerate() {
                if *bit {
                    hit_counts[color] += 1;
                    if is_singleton {
                        unique_hit_counts[color] += 1;
                    }
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
        PseudoalignmentData{hit_counts, unique_hit_counts, distinguishing_hit_counts, n_unique_kmers, n_relevant_kmers, n_all_kmers}
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

    pub fn sbwt(&self) -> &SbwtIndex<SubsetMatrix> {
        &self.kmers
    }

    pub fn build_sbwt_select_support(&mut self) {
        self.kmers.build_select();
    }

    pub fn compress_colors(mut self, sample_distance: usize, n_threads: usize) -> CompactColexColoring {
        self.build_sbwt_select_support(); // Required in the compactification
        let sbwt = Arc::new(self.kmers); // Move to heap
        CompactColexColoring::new(sbwt, &self.distinct_color_sets, self.n_colors, sample_distance, n_threads)
    }
}

struct InputStream {
    dbs: Vec<jseqio::seq_db::SeqDB>,
    cur_db_idx: usize, // Index of the db currently being iterated over
    seq_idx_in_cur_db: usize,
}

impl InputStream {
    fn new<P: AsRef<Path>>(filenames: &[P]) -> InputStream {
        let mut dbs: Vec<jseqio::seq_db::SeqDB> = vec![];
        for path in filenames {
            let reader = jseqio::reader::DynamicFastXReader::from_file(path).unwrap();
            let (mut fw, rc) = reader.into_db_with_revcomp().unwrap();

            if fw.sequence_count() == 0 {
                panic!("No sequences found in file {}", path.as_ref().display());
            }

            // Append reverse complement records to the forward database
            for rec in rc.iter() {
                fw.push_record(rec);
            }
            dbs.push(fw);
        }
        Self {dbs, cur_db_idx: 0, seq_idx_in_cur_db: 0}
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

fn mark_bits(bv: &mut BitVec, color: usize, num_colors: usize, to_mark: Vec<usize>) {
    for i in to_mark {
        bv.set(i*num_colors + color, true);
    }

}

fn mark_all_kmers_of_seq(bv: Arc<Mutex<BitVec>>, num_colors: usize, color: usize, seq: &[u8], k: usize, mark_buffer_size: usize, index: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>){
    // Search all k-mers
    let mut marking_buffer: Vec<usize> = Vec::new(); // These bits should be marked
    for (len, colex) in index.matching_statistics(seq) {
        if len == k {
            // Full kmer -> set the bit in the color set of the k-mer
            assert!(colex.len() == 1);
            marking_buffer.push(colex.start);
            if marking_buffer.len() == mark_buffer_size {
                mark_bits(bv.lock().unwrap().deref_mut(), color, num_colors, marking_buffer);
                marking_buffer = Vec::new();
            }
        }
    }

    if !marking_buffer.is_empty() { 
        // Mark the rest
        mark_bits(bv.lock().unwrap().deref_mut(), color, num_colors, marking_buffer);
    }
} 

impl ColoredKmers {

    #[allow(clippy::type_complexity)]
    pub fn new<P: AsRef<Path> + Send + Sync>(filenames: &[P], k: usize, n_threads: usize, temp_dir: &Path) -> Self {

        log::info!("Loading {} sequence files (colors) into memory", filenames.len());
        let input_stream = InputStream::new(filenames);
        let num_colors = input_stream.dbs.len();
        log::info!("Building SBWT");
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
        let lcs = lcs.unwrap(); // Ok since used build_lcs(true) above

        let sbwt_len = sbwt.n_sets();
        let streaming_index_owned = StreamingIndex::new(&sbwt, &lcs);
        let streaming_index = &streaming_index_owned; // Pass by reference into the scope

        log::info!("Reading input sequences back into memory again");
        let input_stream = InputStream::new(filenames); // Read the data into memory again
        let dbs = &input_stream.dbs;

        let color_sets = std::thread::scope(|scope| {

            log::info!("Building colors");

            let work_input_queue: (Sender<(usize, &[u8])>, Receiver<(usize, &[u8])>) = crossbeam::channel::unbounded();

            // Push work to the input queue
            for color in 0..filenames.len() {
                for rec in dbs[color].iter() {
                    work_input_queue.0.send((color, rec.seq)).unwrap(); 
                }
            }
            drop(work_input_queue.0); // Close the channel

            // Spawn worker threads
            let mut worker_handles = Vec::new();
            let color_sets_lock = Arc::new(Mutex::new(bitvec![0; num_colors*sbwt_len])); // Concatenation of color sets
            for thread_id in 0..n_threads {
                let recv_clone = work_input_queue.1.clone();
                let color_sets_lock_clone = color_sets_lock.clone();
                let consumer_handle = scope.spawn(move || {
                    loop {
                        match recv_clone.recv() {
                            Ok((color, seq)) => {
                                mark_all_kmers_of_seq(color_sets_lock_clone.clone(), num_colors, color, seq, k, 100000, streaming_index);
                            },
                            Err(RecvError) => {
                                log::info!("Thread {} finished", thread_id);
                                break;
                            }
                        }
                    }
                });
                worker_handles.push(consumer_handle);
            }

            // Wait for all workers to finish
            for handle in worker_handles {
                handle.join().unwrap();
            }

            // Since we have joined the workers, there should be only one clone of the
            // Arc<Mutex> (the one owned by this thread), so we can consume the lock and return the data.
            Arc::try_unwrap(color_sets_lock).unwrap().into_inner().unwrap()

        }); // End of thread scope 

        // Todo: deduplicate color sets

        let colex_to_color_set_id: Vec<usize> = (0..sbwt_len).collect(); // Identity mapping

        ColoredKmers{kmers: sbwt, lcs, distinct_color_sets: color_sets, empty_set: BitVec::new(), colex_to_color_set_id, n_colors: filenames.len()}
    }
}

#[cfg(test)]
mod tests {
    use sbwt::BitPackedKmerSorting;

    use super::*;

    #[test]
    fn test_distinguishing_scores(){
        let red_seq: &[u8] = b"AAACATCGATCGTACGTACGTCAGCTACTGCA";
        let blue_seq: &[u8] = b"CACTCTATCGCGTTATCTTACGATCATGCTAGC";
        let green_seq: &[u8] = b"ACATCGGCGTATCTATCTACGATCGTACGTCA";
        let uncolored_seq: &[u8] = b"GGATTCGGATCTATCGTAGCTGTACGTGCTGAC";
        let red_and_blue_seq: &[u8] = b"TTAGCTATCGTATCCGATCACGTACGTAGTCAA";
        let red_and_blue_and_green_seq: &[u8] = b"CCGTTATCGGCCTATACTATCGACTACGTAGC";
        let k = 12; // Let's hope I didn't accidentally repeat any kmers

        let all_seqs: &[&[u8]] = &[red_seq, blue_seq, green_seq, uncolored_seq, red_and_blue_seq, red_and_blue_and_green_seq];
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
