#![allow(non_snake_case)]

use std::{fs::File, io::{BufRead, BufReader}, path::Path};
use jseqio::{reverse_complement, seq_db::SeqDB};
use sha1::{Sha1, Digest};
use rayon::prelude::*;

fn ascii_to_int(ascii: &[u8]) -> usize {
    std::str::from_utf8(ascii)
    .expect("Unitig id is not valid utf-8").parse()
    .expect("Could not convert unitig id string to unsigned integer")
}

fn get_color_set_id(fasta_header: &[u8]) -> usize {
    // The fasta header should look like this " unitig_id=0 color_set_id=0".
    // Note the space at the start.

    let part = fasta_header[1..].split(|c| *c == b' ').nth(1).expect("Color set id missing");
    let mut tokens = part.split(|c| *c == b'=');
    assert_eq!(tokens.next().expect("Color set id missing"), b"color_set_id");
    ascii_to_int(tokens.next().expect("Color set id missing"))
}

fn read_color_sets(filename: impl AsRef<Path>, num_color_sets: usize) -> Vec<Vec<usize>> {
    // Lines should look like this:
    // color_set_id=9 size=7 3 4 9 12 14 15 16

    let mut color_sets = vec![vec![]; num_color_sets];

    let mut reader = BufReader::new(File::open(filename).unwrap());
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
        let list_len = ascii_to_int(&second_token[5..]);

        color_sets[color_set_id] = tokens.map(ascii_to_int).collect();
        assert_eq!(color_sets[color_set_id].len(), list_len);

        line.clear();
        bar.inc(1);
    }
    bar.finish();

    color_sets
}

fn canonicalize_rotation_of_cyclic_unitig(unitig: &mut Vec<u8>, k: usize) {

    assert!(unitig.len() >= k);

    // Find the smallest k-mer
    let min_kmer_start = (0..(unitig.len()-k+1)).min_by(|&i,&j| unitig[i..i+k].cmp(&unitig[j..j+k])).unwrap();

    // Extend unitig so that all unitig rotations are substrings of it
    let original_unitig_len = unitig.len();
    for i in 0..original_unitig_len {
        unitig.push(unitig[k-1+i]);
    }

    // Shift the smallest rotation to the start
    for i in 0..original_unitig_len {
        unitig[i] = unitig[min_kmer_start + i];
    }
    unitig.truncate(original_unitig_len);

}

fn canonicalize_unitig(unitig: &mut Vec<u8>, k: usize) {
    if unitig[0..k-1] == unitig[unitig.len()-(k-1) ..] {
        // Cyclic unitig
        canonicalize_rotation_of_cyclic_unitig(unitig, k);
    }

    let rc = jseqio::reverse_complement(unitig);
    if rc < *unitig {
        jseqio::reverse_complement_in_place(unitig);
    }
} 

fn read_unitigs(filename: impl AsRef<Path>) -> SeqDB {
    jseqio::reader::DynamicFastXReader::from_file(&filename).unwrap().into_db().unwrap()
}

fn read_and_canonicalize_unitigs(filename: impl AsRef<Path>, k: usize) -> SeqDB {
    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&filename).unwrap();
    let mut db = jseqio::seq_db::SeqDB::new();

    while let Some(rec) = reader.read_next().unwrap() {
        assert!(rec.seq.len() >= k);

        let mut rec = rec.to_owned();
        canonicalize_unitig(&mut rec.seq, k);
        db.push_record(rec);
    }

    db

}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
struct Metadata {
    num_unitigs: usize,
    num_colors: usize,
    num_color_sets: usize,
    k: usize
}

fn read_metadata(filename: impl AsRef<Path>) -> Metadata {

    // File should look like this:
    // num_colors=3682
    // num_unitigs=9314735
    // num_color_sets=5591009
    // k=31

    let mut reader = BufReader::new(File::open(filename).unwrap());
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

    Metadata {
        num_unitigs: num_unitigs.expect("num_unitigs missing from metadata"),
        num_colors: num_colors.expect("num_colors missing from metadata"),
        num_color_sets: num_color_sets.expect("num_color_sets missing from metadata"),
        k: k.expect("k missing from metadata")
    }

}

fn sha1(data: &[u8]) -> [u8; 20] {
    let mut hasher = Sha1::new();
    hasher.update(data);
    hasher.finalize().into()
}

fn xor_into<const N: usize>(target: &mut[u8; N], other: &[u8; N]) {
    for i in 0..N {
        target[i] ^= other[i];
    }
}

// This should only be called for odd k because otherwise the rev. comp. of a k-mer may be equal to itself
// Also returns the number of k-mers that were hashed.
fn unitig_checksum(unitig: &[u8], k: usize) -> ([u8; 20], usize) {
    let S = unitig.to_vec();
    let Srev = reverse_complement(unitig);

    let n = S.len();
    assert!(n >= k);
    assert!(Srev.len() == n);
    assert!(k % 2 == 1); // Only works for odd k

    let mut checksum = [0_u8; 20];
    let mut n_hashes = 0_usize;

    for i in 0..(n-k+1) {
        let fw = &S[i..i+k];
        let rc = &Srev[n-k-i..n-i];
        let is_canonical = fw <= rc;

        let hash = if is_canonical {
            sha1(fw)
        } else {
            sha1(rc)
        };

        xor_into(&mut checksum, &hash);

        n_hashes += 1;

    }

    (checksum, n_hashes)
}

fn checksum_unitig_db(unitigs: &SeqDB, k: usize) -> ([u8; 20], usize) {
    let mut checksum = [0_u8; 20];
    let mut n_hashed = 0_usize;

    let bar = indicatif::ProgressBar::new(unitigs.sequence_count() as u64);
    for unitig in unitigs.iter() {
        bar.inc(1);
        let (unitig_checksum, hash_count) = unitig_checksum(&unitig.seq, k);
        xor_into(&mut checksum, &unitig_checksum);
        n_hashed += hash_count;
    }
    bar.finish();

    (checksum, n_hashed)
}

fn hash_color_set(color_set: &[usize]) -> [u8; 20] {
    let mut checksum = [0_u8; 20];
    for color in color_set {
        let hash = sha1(&color.to_le_bytes());
        xor_into(&mut checksum, &hash);
    }
    checksum
}

fn sha1_kmer_into_color_set_hash(kmer: &[u8], color_set_hash: &[u8; 20]) -> [u8; 20] {
    let mut hasher = Sha1::new();
    hasher.update(kmer);
    hasher.update(color_set_hash);
    hasher.finalize().into()
}

// Assumes unitigs are a disjoint spectrum preserving string set (= no duplicate k-mers)
fn checksum_unitig_kmers_and_colorsets(unitig: &[u8], color_set_hash: &[u8; 20], k: usize) -> [u8; 20] {

    let n = unitig.len();
    let unitig = unitig.to_owned();
    let unitig_rc = reverse_complement(&unitig);

    let mut checksum = [0_u8; 20];

    for i in 0..(n-k+1) {
        let fw = &unitig[i..i+k];
        let rc = &unitig_rc[n-k-i..n-i];

        let is_canonical = fw <= rc;

        let combined_hash = if is_canonical {
            sha1_kmer_into_color_set_hash(fw, color_set_hash)
        } else {
            sha1_kmer_into_color_set_hash(rc, color_set_hash)
        };
        xor_into(&mut checksum, &combined_hash);
    }
    checksum
}

// Non-canonical ignored in A
fn compare_color_sets(A_unitigs: &SeqDB, B_unitigs: &SeqDB, A_color_sets: &[Vec<usize>], B_color_sets: &[Vec<usize>], k: usize) {
    eprintln!("Hashing A color sets...");
    let A_color_set_hashes = A_color_sets.par_iter().map(|color_set| hash_color_set(color_set)).collect::<Vec<_>>();

    eprintln!("Hashing B color sets...");
    let B_color_set_hashes = B_color_sets.par_iter().map(|color_set| hash_color_set(color_set)).collect::<Vec<_>>();

    let mut A_checksum = [0_u8; 20];
    let mut B_checksum = [0_u8; 20];

    eprintln!("Computing A checksum...");
    let A_bar = indicatif::ProgressBar::new(A_unitigs.sequence_count() as u64);
    for rec in A_unitigs.iter() {
        let color_set_hash = A_color_set_hashes[get_color_set_id(rec.head)];
        let checksum = checksum_unitig_kmers_and_colorsets(rec.seq, &color_set_hash, k);
        xor_into(&mut A_checksum, &checksum);
        A_bar.inc(1);
    }
    A_bar.finish();

    eprintln!("Computing B checksum...");
    let B_bar = indicatif::ProgressBar::new(B_unitigs.sequence_count() as u64);
    for rec in B_unitigs.iter() {
        let color_set_hash = B_color_set_hashes[get_color_set_id(rec.head)];
        let checksum = checksum_unitig_kmers_and_colorsets(rec.seq, &color_set_hash, k);
        xor_into(&mut B_checksum, &checksum);
        B_bar.inc(1);
    }
    B_bar.finish();

    assert_eq!(A_checksum, B_checksum);
    println!("Checksums match: {:?}", A_checksum);
}

fn main() {
    let mut args = std::env::args();
    args.next().unwrap(); // Program name
    let dump_A_file_prefix = args.next().unwrap();
    let dump_B_file_prefix = args.next().unwrap();

    eprintln!("Reading metadata...");
    let A_metadata = read_metadata(format!("{}.metadata.txt", dump_A_file_prefix));
    let B_metadata = read_metadata(format!("{}.metadata.txt", dump_B_file_prefix));
    assert_eq!(A_metadata.k, B_metadata.k);
    let k = A_metadata.k; 

    eprintln!("Reading unitigs...");
    // We're not canonicalizing the unitigs because the comparison algorithm works anyway.
    let A_unitigs = read_unitigs(format!("{}.unitigs.fa", dump_A_file_prefix));
    let B_unitigs = read_unitigs(format!("{}.unitigs.fa", dump_B_file_prefix)); 

    eprintln!("Computing k-mer checksums...");

    let (A_checksum, A_kmer_count) = &checksum_unitig_db(&A_unitigs, A_metadata.k);
    let (B_checksum, B_kmer_count) = &checksum_unitig_db(&B_unitigs, B_metadata.k);

    assert_eq!(A_kmer_count, B_kmer_count);
    eprintln!("Canonical k-mer counts match: {}", A_kmer_count);

    assert_eq!(A_checksum, B_checksum);
    eprintln!("Canonical k-mer checksums match: {:?}", A_checksum);

    eprintln!("Reading color sets...");
    let A_color_sets = read_color_sets(format!("{}.color_sets.txt", dump_A_file_prefix), A_metadata.num_color_sets);
    let B_color_sets = read_color_sets(format!("{}.color_sets.txt", dump_B_file_prefix), B_metadata.num_color_sets);

    eprintln!("Comparing k-mer color sets...");
    compare_color_sets(&A_unitigs, &B_unitigs, &A_color_sets, &B_color_sets, k);

}


#[cfg(test)]
mod tests {
    use super::*;

    /*
    #[test]
    fn test_unitig_checksum(){
        let S = b"AACTCATACAGCTCTACTTACGACTGCGTCTACTGCTAGCTACA"; // Canonical unitig
        let Srev = reverse_complement(S); // Non-canonical version

        let k = 5;

        let (X, _) = unitig_checksum(S, k, false); // Hashes the canonical version of all k-mers
        
        let (mut Y, _) = unitig_checksum(S, k, true); // Hashes only k-mers that are already canonical 
        let (YRev, _) = unitig_checksum(&Srev, k, true); // Hashes only k-mers that are already canonical 
        xor_into(&mut Y, &YRev);

        assert_eq!(X, Y);
    }

    #[test]
    fn test_cyclic_unitig_rotation(){
        let k = 5;
        let mut S = b"ATCTATCGAAATCACACACTGATCT".to_vec(); // Cyclic: first (k-1)-mer is equal to last (k-1)-mer
        canonicalize_rotation_of_cyclic_unitig(&mut S, k);
        eprintln!("{}", std::str::from_utf8(&S).unwrap());

        assert_eq!(S, b"AAATCACACACTGATCTATCGAAAT");
    }
    */


}
