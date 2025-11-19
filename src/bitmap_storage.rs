use std::{ops::DerefMut, path::Path, sync::{Arc, Mutex}};

use crossbeam::channel::{Receiver, RecvError, Sender};
use sbwt::{self, LcsArray, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use bitvec::prelude::*;

use crate::coloring_interface::{self, ColorSetStream};

/*
 *
 * Structs
 * 
 */

#[derive(Debug)]
pub struct BitmapStorage {
    pub bitmap: BitVec, // Concatenation of distinct color sets
    pub n_colors: usize,
}

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct BitSetView<'storage> {
    bs: &'storage BitSlice<usize, Lsb0>,
}
pub struct BitSetViewIter<'storage> {
    it: bitvec::slice::IterOnes<'storage, usize, Lsb0>,
}

#[derive(Debug, Clone)]
pub struct BitSetOwned {
    bv: BitVec<usize, Lsb0>,
}
pub struct BitSetOwnedIter<'a> {
    it: bitvec::slice::IterOnes<'a, usize, Lsb0>,
}


/*
 *
 * Trait implementations for above structs 
 * 
 */


impl<'storage> crate::coloring_interface::ColorSetView<'storage> for BitSetView<'storage> {
    type Iter = BitSetViewIter<'storage>;

    fn iter(&self) -> Self::Iter {
        BitSetViewIter{it: self.bs.iter_ones()}
    }

    fn len(&self) -> usize {
        self.bs.count_ones()
    }
}

impl<'a> Iterator for BitSetViewIter<'a> {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        self.it.next()
    }

}

impl coloring_interface::ColorSetOwned for BitSetOwned {
    type Iter<'a> = BitSetOwnedIter<'a> where Self: 'a;

    fn iter(&self) -> Self::Iter<'_> {
        BitSetOwnedIter{it: self.bv.iter_ones()}
    }
}

impl<'a> Iterator for BitSetOwnedIter<'a> {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        self.it.next()
    }
}

impl crate::coloring_interface::ColorSetStorage for BitmapStorage {
    type SetView<'a> = BitSetView<'a> where Self: 'a;
    type OwnedSet = BitSetOwned;

    fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow> {
        BitSetView{bs: &self.bitmap[id*self.n_colors..(id+1)*self.n_colors]}
    }

    fn new(mut sets: impl ColorSetStream, n_colors: usize) -> Box<Self> {
        let mut bitmap = bitvec![];
        let empty_set = bitvec![0; n_colors];
        let mut id = 0_usize;

        while let Some(set) = sets.next() {
            bitmap.extend_from_bitslice(&empty_set);
            for color in set {
                bitmap.set(id*n_colors + color, true);
            }
            id += 1;
        }

        Box::new(Self{bitmap, n_colors})
    }
    
    fn new_from_transpose(mut set_ids_per_color: impl ColorSetStream, n_colors: usize, _set_sizes: &Vec<u64>) -> Box<Self> {
        // The set sizes are not needed.

        let mut bitmap = bitvec![];
        let empty_set = bitvec![0; n_colors];
        let mut color_id = 0_usize;
        while let Some(set_ids) = set_ids_per_color.next() {
            bitmap.extend_from_bitslice(&empty_set);
            for set_id in set_ids {
                bitmap.set(set_id*n_colors + color_id, true);
            }
            color_id += 1;
        }

        Box::new(Self{bitmap, n_colors})
    }

    fn get_empty_set(&self) -> Self::OwnedSet {
        BitSetOwned{bv: bitvec![]}
    }

    fn get_full_set(&self) -> Self::OwnedSet {
        BitSetOwned{bv: bitvec![usize, Lsb0; 1; self.n_colors]}
    }

    fn serialize<W: std::io::Write>(&self, out: &mut W) {
        out.write_all(&(self.n_colors as u64).to_le_bytes()).unwrap();
        out.write_all(&(self.bitmap.len() as u64).to_le_bytes()).unwrap();
        bincode::serialize_into(out, &self.bitmap).unwrap();
    }

    fn load<R: std::io::Read>(input: &mut R) -> Self {
        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf).unwrap();
        let n_colors = u64::from_le_bytes(buf) as usize;

        input.read_exact(&mut buf).unwrap();
        let _ = u64::from_le_bytes(buf); // Total length of distinct color sets

        let bitmap: BitVec = bincode::deserialize_from(input).unwrap();

        Self{n_colors, bitmap}
    }

    fn view_to_owned(&self, view: &Self::SetView<'_>) -> Self::OwnedSet {
        BitSetOwned{bv: view.bs.to_bitvec()}
    }

    fn owned_to_view<'a>(&self, owned: &'a Self::OwnedSet) -> Self::SetView<'a> {
        BitSetView{bs: &owned.bv}
    }
    
    fn intersect(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
        a.bv &= b.bs;
    }
    
    fn union(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
        a.bv |= b.bs;
    }
    
    fn n_sets(&self) -> usize {
        self.bitmap.len() / self.n_colors
    }
}

/*
 *
 * Construction algorithm
 *
 */

struct InputStream {
    dbs: Arc<Vec<jseqio::seq_db::SeqDB>>, // Arc because we may want to hold onto the dbs even if the stream is consumed
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
        Self {dbs: Arc::new(dbs), cur_db_idx: 0, seq_idx_in_cur_db: 0}
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

/// Note: reverse complements are not added, so if you want them, include them in the dbs.
pub fn build_from_seq_dbs(dbs: Vec<jseqio::seq_db::SeqDB>, sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, n_threads: usize) -> BitmapStorage {
    let dbs = Arc::new(dbs);
    let input_stream = InputStream {
        dbs: dbs.clone(),
        cur_db_idx: 0,
        seq_idx_in_cur_db: 0,
    };

    let num_colors = input_stream.dbs.len();

    let sbwt_len = sbwt.n_sets();
    let k = sbwt.k();
    let streaming_index_owned = StreamingIndex::new(sbwt, lcs);
    let streaming_index = &streaming_index_owned; // Pass by reference into the scope

    // Stream the output again to mark colors
    let color_sets = std::thread::scope(|scope| {

        log::info!("Building colors");

        #[allow(clippy::type_complexity)]
        let work_input_queue: (Sender<(usize, &[u8])>, Receiver<(usize, &[u8])>) = crossbeam::channel::unbounded();

        // Push work to the input queue
        for color in 0..dbs.len() {
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

    BitmapStorage{
        bitmap: color_sets,
        n_colors: num_colors,
    }
}

#[allow(clippy::type_complexity)]
pub fn build_from_files<P: AsRef<Path> + Send + Sync>(filenames: &[P], sbwt: &SbwtIndex<SubsetMatrix>, lcs: &LcsArray, n_threads: usize) -> BitmapStorage {

    log::info!("Loading {} sequence files (colors) into memory", filenames.len());
    let dbs = Arc::try_unwrap(InputStream::new(filenames).dbs).ok().unwrap(); // Also appends reverse complements to the dbs

    log::info!("Indexing");
    build_from_seq_dbs(dbs, sbwt, lcs, n_threads)
}

#[cfg(test)]
mod tests {

    use super::*;

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
