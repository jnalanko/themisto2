//! Reading and writing Themisto 2 indexes.
//!
//! A Themisto 2 index is a directory with one file per component, so that each operation can load
//! only the components it needs, and the sizes of the files tell how much memory loading them takes:
//!
//! - `metadata.json`: [Metadata]. Written last, so a directory without it is not a complete index.
//! - `kmers.sbwt`: the SBWT, in the format of sbwt-rs-cli ([sbwt::sbwt_index_variant::SbwtIndexVariant]).
//! - `kmers.lcs`: the LCS array of the SBWT.
//! - `kmers.sufgroups`: the (k-1)-suffix group marks of the SBWT ([sbwt::dbg::Dbg::k_minus_1_marks]).
//! - `kmers.dummymarks`: the dummy k-mer marks of the SBWT ([sbwt::dbg::Dbg::dummy_marks]).
//! - `colors.names`: the color names, one per line.
//! - `colors.sampling`: the color set ids of the sampled k-mers.
//! - `colors.storage`: the distinct color sets.
//!
//! All components except the metadata are optional: loading a missing one returns an error that
//! names it.

use std::fs::File;
use std::io::{BufReader, BufWriter, Error, ErrorKind, Write};
use std::path::{Path, PathBuf};

use sbwt::sbwt_index_variant::SbwtIndexVariant;
use sbwt::{LcsArray, SbwtIndex, SubsetMatrix};
use serde::{Deserialize, Serialize};
use simple_sds_sbwt::ops::BitVec;

use crate::colex_colored_kmers::{ColexToColorSetMap, Coloring, DbgIndex};
use crate::coloring_interface::ColorSetStorage;
use crate::sparse_dense_storage::SparseDenseStorage;

/// The version of the index directory format, stored in [Metadata::format_version].
pub const FORMAT_VERSION: u64 = 1;

const METADATA_FILE: &str = "metadata.json";

/// The contents of `metadata.json`.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct Metadata {
    pub format_version: u64,
    pub k: usize,
    pub n_kmers: usize,
    /// The number of positions in the SBWT, including dummy k-mers
    pub n_sbwt_sets: usize,
    pub n_colors: usize,
    pub n_distinct_color_sets: usize,
    /// Whether the index has both orientations of every k-mer. None if not known (an index
    /// imported from an SBWT built elsewhere).
    pub both_orientations: Option<bool>,
    /// The subset rank structure of the SBWT, as named in sbwt-rs-cli
    pub sbwt_variant: String,
}

impl Metadata {
    pub fn new<CSS: ColorSetStorage>(sbwt: &SbwtIndex<SubsetMatrix>, coloring: &Coloring<CSS>, both_orientations: Option<bool>) -> Self {
        Self {
            format_version: FORMAT_VERSION,
            k: sbwt.k(),
            n_kmers: sbwt.n_kmers(),
            n_sbwt_sets: sbwt.n_sets(),
            n_colors: coloring.n_colors(),
            n_distinct_color_sets: coloring.get_set_storage().n_sets(),
            both_orientations,
            sbwt_variant: SUBSET_MATRIX_ID.to_string(),
        }
    }
}

const SUBSET_MATRIX_ID: &str = "SubsetMatrix";

/// A component of an index: one file in the index directory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Component {
    Sbwt,
    Lcs,
    Sufgroups,
    DummyMarks,
    ColorNames,
    ColorSampling,
    ColorStorage,
}

impl Component {
    pub const ALL: [Component; 7] = [Component::Sbwt, Component::Lcs, Component::Sufgroups, Component::DummyMarks, Component::ColorNames, Component::ColorSampling, Component::ColorStorage];

    pub fn file_name(&self) -> &'static str {
        match self {
            Component::Sbwt => "kmers.sbwt",
            Component::Lcs => "kmers.lcs",
            Component::Sufgroups => "kmers.sufgroups",
            Component::DummyMarks => "kmers.dummymarks",
            Component::ColorNames => "colors.names",
            Component::ColorSampling => "colors.sampling",
            Component::ColorStorage => "colors.storage",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Component::Sbwt => "SBWT",
            Component::Lcs => "LCS array",
            Component::Sufgroups => "suffix group marks",
            Component::DummyMarks => "dummy k-mer marks",
            Component::ColorNames => "color names",
            Component::ColorSampling => "color set sampling",
            Component::ColorStorage => "color set storage",
        }
    }
}

fn invalid_data(msg: String) -> Error {
    Error::new(ErrorKind::InvalidData, msg)
}

/// An index directory opened for reading. Opening reads only the metadata, and each component is
/// loaded separately.
pub struct IndexDir {
    path: PathBuf,
    metadata: Metadata,
}

impl IndexDir {

    /// Opens the index directory at `path` and reads its metadata.
    pub fn open(path: &Path) -> std::io::Result<Self> {
        if path.is_file() {
            return Err(invalid_data(format!("{} is a file, but a Themisto 2 index is a directory. Indexes in the single-file format of earlier versions are not supported: rebuild the index", path.display())));
        }
        let metadata_path = path.join(METADATA_FILE);
        let file = File::open(&metadata_path).map_err(|e| {
            Error::new(e.kind(), format!("{} is not a Themisto 2 index: could not open {}: {}", path.display(), metadata_path.display(), e))
        })?;
        let metadata: Metadata = serde_json::from_reader(BufReader::new(file))
            .map_err(|e| invalid_data(format!("Could not parse {}: {}", metadata_path.display(), e)))?;
        if metadata.format_version != FORMAT_VERSION {
            return Err(invalid_data(format!("The Themisto 2 index at {} has format version {}, but this version of Themisto 2 reads version {}", path.display(), metadata.format_version, FORMAT_VERSION)));
        }
        Ok(Self { path: path.to_path_buf(), metadata })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }

    pub fn component_path(&self, component: Component) -> PathBuf {
        self.path.join(component.file_name())
    }

    pub fn has(&self, component: Component) -> bool {
        self.component_path(component).is_file()
    }

    /// The size in bytes of the file of each component, or None if the index does not have it.
    pub fn component_sizes(&self) -> Vec<(Component, Option<u64>)> {
        Component::ALL.iter().map(|&c| (c, std::fs::metadata(self.component_path(c)).ok().map(|m| m.len()))).collect()
    }

    fn open_component(&self, component: Component) -> std::io::Result<BufReader<File>> {
        let path = self.component_path(component);
        match File::open(&path) {
            Ok(f) => {
                let size = f.metadata().map(|m| m.len()).unwrap_or(0);
                log::info!("Loading the {} ({}, {})", component.description(), path.display(), human_bytes::human_bytes(size as f64));
                Ok(BufReader::new(f))
            },
            Err(e) if e.kind() == ErrorKind::NotFound => Err(Error::new(ErrorKind::NotFound, format!(
                "The Themisto 2 index at {} has no {} ({}), which this operation needs",
                self.path.display(), component.description(), component.file_name()))),
            Err(e) => Err(Error::new(e.kind(), format!("Could not open {}: {}", path.display(), e))),
        }
    }

    fn check_len(&self, component: Component, len: usize) -> std::io::Result<()> {
        if len != self.metadata.n_sbwt_sets {
            return Err(invalid_data(format!("The {} of the Themisto 2 index at {} has length {}, but the SBWT has {} positions",
                component.description(), self.path.display(), len, self.metadata.n_sbwt_sets)));
        }
        Ok(())
    }

    /// Loads the SBWT. Select support is needed for de Bruijn graph operations,
    /// set operations and exporting, but not for pseudoalignment.
    pub fn load_sbwt(&self, build_select: bool) -> std::io::Result<SbwtIndex<SubsetMatrix>> {
        let mut input = self.open_component(Component::Sbwt)?;
        let variant = SbwtIndexVariant::load(&mut input)
            .map_err(|e| invalid_data(format!("Could not load the SBWT of the Themisto 2 index at {}: {}", self.path.display(), e)))?;
        let mut sbwt = match variant {
            SbwtIndexVariant::SubsetMatrix(sbwt) => sbwt,
            _ => return Err(invalid_data(format!("The SBWT of the Themisto 2 index at {} is not of type {}, the only type Themisto 2 supports", self.path.display(), SUBSET_MATRIX_ID))),
        };
        self.check_len(Component::Sbwt, sbwt.n_sets())?;
        if build_select {
            log::info!("Building select support");
            sbwt.build_select();
        }
        Ok(sbwt)
    }

    pub fn load_lcs(&self) -> std::io::Result<LcsArray> {
        let lcs = LcsArray::load(&mut self.open_component(Component::Lcs)?)?;
        self.check_len(Component::Lcs, lcs.len())?;
        Ok(lcs)
    }

    pub fn load_sufgroups(&self) -> std::io::Result<bitvec::vec::BitVec> {
        let bits = crate::util::read_bitvec(&mut self.open_component(Component::Sufgroups)?)?;
        self.check_len(Component::Sufgroups, bits.len())?;
        Ok(bits)
    }

    pub fn load_dummy_marks(&self) -> std::io::Result<bitvec::vec::BitVec> {
        let bits = crate::util::read_bitvec(&mut self.open_component(Component::DummyMarks)?)?;
        self.check_len(Component::DummyMarks, bits.len())?;
        Ok(bits)
    }

    pub fn load_color_names(&self) -> std::io::Result<Vec<String>> {
        let names = Coloring::<SparseDenseStorage>::load_color_names(self.open_component(Component::ColorNames)?)?;
        if names.len() != self.metadata.n_colors {
            return Err(invalid_data(format!("The Themisto 2 index at {} has {} color names, but {} colors", self.path.display(), names.len(), self.metadata.n_colors)));
        }
        Ok(names)
    }

    /// Loads the color names, the color set sampling and the color set storage.
    pub fn load_coloring(&self) -> std::io::Result<Coloring<SparseDenseStorage>> {
        // Open all files first, so that a missing one is reported before loading anything big
        let names_input = self.open_component(Component::ColorNames)?;
        let mut sampling_input = self.open_component(Component::ColorSampling)?;
        let mut storage_input = self.open_component(Component::ColorStorage)?;
        drop(names_input);
        let names = self.load_color_names()?;
        let map = ColexToColorSetMap::load(&mut sampling_input);
        self.check_len(Component::ColorSampling, map.sampling.len())?;
        let sets = SparseDenseStorage::load(&mut storage_input);
        if sets.n_colors() != names.len() {
            return Err(invalid_data(format!("The color set storage of the Themisto 2 index at {} has {} colors, but there are {} color names", self.path.display(), sets.n_colors(), names.len())));
        }
        Ok(Coloring::new(map, sets, Some(&names)))
    }

    /// Loads the components that de Bruijn graph operations and set operations need: the SBWT with
    /// select support, the suffix group marks, the dummy marks and the coloring.
    pub fn load_dbg_index(&self) -> std::io::Result<DbgIndex<SparseDenseStorage>> {
        self.require(&[Component::Sbwt, Component::Sufgroups, Component::DummyMarks, Component::ColorNames, Component::ColorSampling, Component::ColorStorage])?;
        let sbwt = self.load_sbwt(true)?;
        let sufgroups = self.load_sufgroups()?;
        let dummy_marks = self.load_dummy_marks()?;
        let coloring = self.load_coloring()?;
        Ok(DbgIndex::new(sbwt, sufgroups, dummy_marks, coloring))
    }

    /// Returns an error naming the first of `components` that the index does not have. Call
    /// before loading anything, to fail fast.
    pub fn require(&self, components: &[Component]) -> std::io::Result<()> {
        for &c in components {
            if !self.has(c) {
                return Err(Error::new(ErrorKind::NotFound, format!(
                    "The Themisto 2 index at {} has no {} ({}), which this operation needs",
                    self.path.display(), c.description(), c.file_name())));
            }
        }
        Ok(())
    }
}

/// Writes an index directory. The components are written one at a time, and [IndexWriter::finish]
/// writes the metadata last.
pub struct IndexWriter {
    path: PathBuf,
}

impl IndexWriter {

    /// Creates the directory `path`. Returns an error if `path` exists and is not an empty directory.
    pub fn create(path: &Path) -> std::io::Result<Self> {
        if path.exists() {
            if !path.is_dir() {
                return Err(Error::new(ErrorKind::AlreadyExists, format!("Can not write a Themisto 2 index to {}: it exists and is not a directory", path.display())));
            }
            if std::fs::read_dir(path)?.next().is_some() {
                return Err(Error::new(ErrorKind::AlreadyExists, format!("Can not write a Themisto 2 index to {}: the directory is not empty", path.display())));
            }
        } else {
            std::fs::create_dir_all(path).map_err(|e| Error::new(e.kind(), format!("Could not create directory {}: {}", path.display(), e)))?;
        }
        Ok(Self { path: path.to_path_buf() })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    fn write_component(&self, component: Component, write: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<()>) -> std::io::Result<()> {
        let path = self.path.join(component.file_name());
        log::info!("Writing the {} to {}", component.description(), path.display());
        let mut out = BufWriter::new(File::create(&path).map_err(|e| Error::new(e.kind(), format!("Could not create {}: {}", path.display(), e)))?);
        write(&mut out)?;
        out.flush()
    }

    /// Writes the SBWT in the format of [SbwtIndexVariant::serialize].
    pub fn write_sbwt(&self, sbwt: &SbwtIndex<SubsetMatrix>) -> std::io::Result<()> {
        self.write_component(Component::Sbwt, |out| {
            // The same bytes as SbwtIndexVariant::serialize, which needs an owned SbwtIndexVariant
            out.write_all(&(SUBSET_MATRIX_ID.len() as u64).to_le_bytes())?;
            out.write_all(SUBSET_MATRIX_ID.as_bytes())?;
            sbwt.serialize(out)?;
            Ok(())
        })
    }

    pub fn write_lcs(&self, lcs: &LcsArray) -> std::io::Result<()> {
        self.write_component(Component::Lcs, |out| lcs.serialize(out).map(|_| ()))
    }

    pub fn write_sufgroups(&self, sufgroups: &bitvec::vec::BitVec) -> std::io::Result<()> {
        self.write_component(Component::Sufgroups, |out| crate::util::write_bitvec(sufgroups, out))
    }

    pub fn write_dummy_marks(&self, dummy_marks: &bitvec::vec::BitVec) -> std::io::Result<()> {
        self.write_component(Component::DummyMarks, |out| crate::util::write_bitvec(dummy_marks, out))
    }

    /// Writes the color names, the color set sampling and the color set storage.
    pub fn write_coloring<CSS: ColorSetStorage>(&self, coloring: &Coloring<CSS>) -> std::io::Result<()> {
        self.write_component(Component::ColorNames, |out| Coloring::<CSS>::serialize_color_names(coloring.get_color_names(), out))?;
        self.write_component(Component::ColorSampling, |out| { coloring.get_map().serialize(out); Ok(()) })?;
        self.write_component(Component::ColorStorage, |out| { coloring.get_set_storage().serialize(out); Ok(()) })
    }

    /// Writes every component of `index`, the LCS array `lcs` of its SBWT, and the metadata.
    pub fn write_all<CSS: ColorSetStorage>(self, index: &DbgIndex<CSS>, lcs: &LcsArray, both_orientations: Option<bool>) -> std::io::Result<()> {
        self.write_sbwt(&index.sbwt)?;
        self.write_lcs(lcs)?;
        self.write_sufgroups(&index.sufgroups)?;
        self.write_dummy_marks(&index.dummy_marks)?;
        self.write_coloring(&index.coloring)?;
        self.finish(&Metadata::new(&index.sbwt, &index.coloring, both_orientations))
    }

    /// Writes the metadata, which completes the index.
    pub fn finish(self, metadata: &Metadata) -> std::io::Result<()> {
        let path = self.path.join(METADATA_FILE);
        let mut out = BufWriter::new(File::create(&path)?);
        serde_json::to_writer_pretty(&mut out, metadata).map_err(|e| Error::other(e.to_string()))?;
        out.write_all(b"\n")?;
        out.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ColorSetView;

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("themisto2_index_dir_test_{}_{}", name, std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    #[test]
    fn write_and_read_roundtrip() {
        let seqs: Vec<(&[u8], usize)> = vec![(b"ACGTTGCATGCAAGTCCGATAG", 0), (b"TGCATGCAAGTCCGTTTTTT", 1)];
        let (index, lcs) = DbgIndex::<SparseDenseStorage>::new_from_small_input_with_lcs(&seqs, 5, 3, 1);

        let dir = temp_dir("roundtrip");
        IndexWriter::create(&dir).unwrap().write_all(&index, &lcs, Some(true)).unwrap();
        assert!(IndexWriter::create(&dir).is_err(), "Writing into a non-empty directory must fail");

        let opened = IndexDir::open(&dir).unwrap();
        assert_eq!(opened.metadata().k, 5);
        assert_eq!(opened.metadata().n_colors, 2);
        assert_eq!(opened.metadata().both_orientations, Some(true));
        assert!(opened.component_sizes().iter().all(|(_, size)| size.is_some()));
        assert_eq!(&opened.load_color_names().unwrap(), index.get_color_names());
        assert_eq!(opened.load_lcs().unwrap(), lcs);

        let loaded = opened.load_dbg_index().unwrap();
        assert_eq!(loaded.sufgroups, index.sufgroups);
        assert_eq!(loaded.dummy_marks, index.dummy_marks);
        assert!(sbwt::SubsetSeq::has_select_support(loaded.sbwt.sbwt()));
        for colex in 0..loaded.sbwt.n_sets() {
            if loaded.dummy_marks[colex] { continue; }
            let a: Vec<usize> = index.colex_to_set(colex).iter().collect();
            let b: Vec<usize> = loaded.colex_to_set(colex).iter().collect();
            assert_eq!(a, b);
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn missing_component_is_an_error() {
        let seqs: Vec<(&[u8], usize)> = vec![(b"ACGTTGCATGCAAGTCCGATAG", 0)];
        let (index, lcs) = DbgIndex::<SparseDenseStorage>::new_from_small_input_with_lcs(&seqs, 5, 3, 1);
        let dir = temp_dir("missing");
        IndexWriter::create(&dir).unwrap().write_all(&index, &lcs, Some(true)).unwrap();
        std::fs::remove_file(dir.join(Component::Sufgroups.file_name())).unwrap();

        let opened = IndexDir::open(&dir).unwrap();
        let err = opened.load_sufgroups().err().unwrap();
        assert_eq!(err.kind(), ErrorKind::NotFound);
        assert!(err.to_string().contains("kmers.sufgroups"), "{}", err);
        assert!(opened.load_dbg_index().is_err());
        assert!(opened.load_sbwt(false).is_ok());

        assert!(IndexDir::open(&dir.join("nonexistent")).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
