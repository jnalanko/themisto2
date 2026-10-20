//! Reading and writing Themisto 2 index files.
//!
//! An index file starts with an 8-byte id of the color set storage type, followed by the
//! serialized [CompactColexKmers]. Only the sparse-dense storage is supported. The id
//! [0,0,0,0,0,0,0,1] was used by an older bitmap storage that no longer exists.

use std::fs::File;
use std::io::{BufReader, Error, ErrorKind, Read, Write};
use std::path::Path;

use crate::colex_colored_kmers::CompactColexKmers;
use crate::sparse_dense_storage::SparseDenseStorage;

/// A Themisto 2 index.
pub type Index = CompactColexKmers<SparseDenseStorage>;

const SPARSE_DENSE_ID: [u8; 8] = [0, 0, 0, 0, 0, 0, 0, 2];

fn read_header(input: &mut impl Read) -> std::io::Result<()> {
    let mut id_buf = [0u8; 8];
    input.read_exact(&mut id_buf)?;
    if id_buf == SPARSE_DENSE_ID {
        Ok(())
    } else {
        Err(Error::new(ErrorKind::InvalidData, format!("Unrecognized index serialization ID: {:?}", id_buf)))
    }
}

/// Loads a Themisto 2 index from a file. Select support is needed for de Bruijn graph
/// operations ([crate::ColoredDbg]), merging and exporting, but not for pseudoalignment.
/// Returns an error if the file can not be opened or is not a Themisto 2 index. Panics if
/// the file is corrupt.
pub fn load_index(path: &Path, build_select: bool) -> std::io::Result<Index> {
    read_index(&mut BufReader::new(File::open(path)?), build_select)
}

/// Like [load_index], but reads from a stream.
pub fn read_index(input: &mut impl Read, build_select: bool) -> std::io::Result<Index> {
    read_header(input)?;
    Ok(Index::load(input, build_select))
}

/// Loads only the color names of a Themisto 2 index without reading the rest of the file.
pub fn load_color_names(path: &Path) -> std::io::Result<Vec<String>> {
    let mut input = BufReader::new(File::open(path)?);
    read_header(&mut input)?;
    Ok(Index::load_color_names_only(&mut input))
}

/// Writes the index in the format read by [load_index].
pub fn write_index(index: &Index, out: &mut impl Write) -> std::io::Result<()> {
    out.write_all(&SPARSE_DENSE_ID)?;
    index.serialize(out);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ColorSetView;

    #[test]
    fn write_and_read_roundtrip() {
        let seqs: Vec<(&[u8], usize)> = vec![(b"ACGTTGCATGCAAGTCCGATAG", 0), (b"TGCATGCAAGTCCGTTTTTT", 1)];
        let index = Index::new_from_small_input(&seqs, 5, 3, 1);

        let mut buf = Vec::<u8>::new();
        write_index(&index, &mut buf).unwrap();
        let loaded = read_index(&mut buf.as_slice(), true).unwrap();

        assert_eq!(loaded.get_k(), index.get_k());
        assert_eq!(loaded.get_color_names(), index.get_color_names());
        let g = crate::ColoredDbg::new(&loaded, 1); // Requires that select support was built
        for v in g.nodes() {
            let a: Vec<usize> = index.colex_to_set(v.id).iter().collect();
            let b: Vec<usize> = loaded.colex_to_set(v.id).iter().collect();
            assert_eq!(a, b);
        }

        let mut bad = buf.clone();
        bad[7] = 1;
        assert_eq!(read_index(&mut bad.as_slice(), false).err().unwrap().kind(), ErrorKind::InvalidData);
    }
}
