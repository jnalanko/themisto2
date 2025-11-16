use std::collections::HashSet;

use simple_sds_sbwt::int_vector::IntVector;
use simple_sds_sbwt::serialize::Serialize;
use simple_sds_sbwt::{ops::{Access, BitVec, Push, Rank, Resize, Vector}, raw_vector::PushRaw};
use bitvec::slice::BitSlice;
use bitvec::bitvec;

use crate::coloring_interface::{ColorSetOwned, ColorSetView};

/*
 *
 * Structs
 * 
 */

/// A data structure for storing arbitary set of sets of integers, such that dense
/// sets are encoded as bitmaps, and sparse sets as lists of integers.
pub struct SparseDenseStorage{
    dense_sets: BitMaps,
    sparse_sets: SortedIntVecs,
    n_colors: usize,
    is_dense_marks: simple_sds_sbwt::bit_vector::BitVector, // Has rank support.
}

// A set of lists of integers, stored in concatenated form. Each
// list is assumed to be sorted.
struct SortedIntVecs {
    // Concatenation of IntVecs. Better keep this private because the sets need
    // to maintained in sorted order for intersections. Since it's private, we only need 
    // to make sure that code in this module respects that. Users can get
    // immutable views to this via IntVecSlice.
    concat: IntVector, 

    // Ends of individual intvecs, such that ends[0] = 0 and ends[i+1] is the
    // exclusive end of the i-th vector.
    ends: Vec<usize>, 
}

#[derive(Copy, Clone)]
pub struct IntVecSlice<'a> {
    vec: &'a IntVector,
    start: usize,
    end: usize, // Exclusive end
}

// A set of sets encoded as bitmaps.
struct BitMaps {
    bitmap_data: bitvec::vec::BitVec, // Concatenation of bit vectors
    individual_length: usize, // Length of each bitmap in bitmap_data
}

#[derive(Clone)]
pub enum SetType {
    Dense(bitvec::vec::BitVec),
    Sparse(IntVector),
}

#[derive(Clone)]
pub struct SparseDenseColorSetOwned {
    set: SetType,
}

// This enum is only for passing references to individual sets around. The actual
// sets are stored in concatenated form somewhere else in memory. 
#[derive(Copy, Clone)]
pub enum SparseDenseColorSetView<'a> {
    Dense(&'a BitSlice),
    Sparse(IntVecSlice<'a>),
}

pub struct ColorSetViewIterator<'a> {
    set: SparseDenseColorSetView<'a>,
    pos: usize, // Interpreted differently depending of whether this is Sparse or Dense
}

/*
 *
 * Trait implementations for above structs 
 * 
 */


impl<'a> Iterator for ColorSetViewIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match &self.set {
            SparseDenseColorSetView::Dense(bit_slice) => {
                // Rewind to the next 1-bit
                if let Some(offset) = bit_slice[self.pos..].first_one(){
                    let ret = self.pos + offset;
                    self.pos = ret + 1; // Starting point for the next iteration
                    Some(ret)
                } else {
                    None
                }
            },
            SparseDenseColorSetView::Sparse(int_vec_slice) => {
                if self.pos == int_vec_slice.end - int_vec_slice.start {
                    None
                } else {
                    let x = int_vec_slice.vec.get(int_vec_slice.start + self.pos);
                    self.pos += 1;
                    Some(x as usize)
                }
            },
        }
    }
}

impl crate::coloring_interface::ColorSetOwned for SparseDenseColorSetOwned {
    type Iter<'a> = ColorSetViewIterator<'a> where Self: 'a;

    fn iter(&self) -> Self::Iter<'_> {
        match &self.set {
            SetType::Dense(bv) => {
                ColorSetViewIterator{
                    set: SparseDenseColorSetView::Dense(bv.as_bitslice()),
                    pos: 0,
                }
            },
            SetType::Sparse(iv) => {
                ColorSetViewIterator{
                    set: SparseDenseColorSetView::Sparse(IntVecSlice{
                        vec: iv,
                        start: 0,
                        end: iv.len(),
                    }),
                    pos: 0,
                }
            },
        }
    }
}

impl SparseDenseColorSetOwned {

    // n_colors-1 is the maximum color id supported by this set
    pub fn new(elements: impl Iterator<Item = usize>, n_colors: usize) -> Self {
        let elements: Vec<usize> = elements.collect();
        let bits_per_color = n_colors.next_power_of_two().trailing_zeros() as usize;
        if is_dense_set(elements.len(), bits_per_color, n_colors) {
            let mut bv = bitvec![0; n_colors];
            for color in elements.iter() {
                bv.set(color, true);
            }
            SparseDenseColorSetOwned {set: SetType::Dense(bv)}
        } else {
            let mut iv = IntVector::new(bits_per_color).unwrap();
            for color in elements.iter() {
                iv.push(color as u64);
            }
            SparseDenseColorSetOwned {set: SetType::Sparse(iv)}
        }
    }
}

impl<'a> crate::coloring_interface::ColorSetView<'a> for SparseDenseColorSetView<'a> {
    type Iter = ColorSetViewIterator<'a>;

    fn iter(&self) -> Self::Iter {
        ColorSetViewIterator{
            set: *self,
            pos: 0,
        }
    }
    
    fn len(&self) -> usize {
        match self {
            SparseDenseColorSetView::Dense(bv) => {
                bv.count_ones()
            },
            SparseDenseColorSetView::Sparse(iv) => {
                iv.end - iv.start
            },
        }
    }
}

impl crate::coloring_interface::ColorSetStorage for SparseDenseStorage {
    type SetView<'a> = SparseDenseColorSetView<'a>; 
    type OwnedSet = SparseDenseColorSetOwned;

    fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow> {
        self.get(id)
    }

    fn get_empty_set(&self) -> Self::OwnedSet {
        let bit_width = self.n_colors.next_power_of_two().trailing_zeros() as usize;
        Self::OwnedSet {
            set: SetType::Sparse(IntVector::new(bit_width).unwrap()),
        }
    }

    fn get_full_set(&self) -> Self::OwnedSet {
        Self::OwnedSet {
            set: SetType::Dense(bitvec![1; self.n_colors]),
        }
    }
    
    fn new(sets: impl Iterator<Item = impl Iterator<Item = usize>>, n_colors: usize) -> SparseDenseStorage {
        log::info!("Encoding color sets");
        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;
        let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

        let mut sparse_sets = SortedIntVecs::new(color_id_bit_width);
        let mut dense_sets = BitMaps::new(n_colors);

        let mut buf = Vec::<usize>::new();
        let mut n_sets_total = 0_usize;
        for set in sets {
            buf.clear();
            buf.extend(set);
            if is_dense_set(buf.len(), color_id_bit_width, n_colors) {
                let mut bm = bitvec![0; n_colors];
                for color in buf.iter() {
                    bm.set(color, true);
                }
                dense_sets.push(&bm);
                is_dense_marks.push_bit(true);
            } else {
                buf.sort_unstable(); // We need sorted sets for intersections
                sparse_sets.push(buf.iter());
                is_dense_marks.push_bit(false);
            }

            n_sets_total += 1;
        }

        sparse_sets.shrink_to_fit();
        dense_sets.shrink_to_fit();

        log::info!("{}% of the sets are sparse", sparse_sets.n_sets() as f64 / n_sets_total as f64 * 100.0);

        // Add rank support to dense marks
        log::info!("Building rank support for dense marks");
        let mut is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::from(is_dense_marks);
        is_dense_marks.enable_rank();

        SparseDenseStorage {
            is_dense_marks, 
            sparse_sets,
            dense_sets,
            n_colors
        }
    }

    fn serialize<W: std::io::Write>(&self, out: &mut W) {
        bincode::serialize_into(out.by_ref(), &self.n_colors).unwrap();
        self.is_dense_marks.serialize(out).unwrap();
        self.sparse_sets.serialize(out);
        self.dense_sets.serialize(out);
    }

    fn load<R: std::io::Read>(input: &mut R) -> Self {
        let n_colors: usize = bincode::deserialize_from(input.by_ref()).unwrap();
        let is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::load(input).unwrap();
        let sparse_sets = SortedIntVecs::load(input);
        let dense_sets = BitMaps::load(input);

        assert_eq!(is_dense_marks.len(), sparse_sets.n_sets() + dense_sets.n_sets());
        assert_eq!(n_colors, dense_sets.individual_length);
        assert!(sparse_sets.concat.width() >= n_colors.next_power_of_two().trailing_zeros() as usize);

        Self {is_dense_marks, sparse_sets, dense_sets, n_colors}
    }
    
    fn view_to_owned(&self, view: &Self::SetView<'_>) -> Self::OwnedSet {
        SparseDenseColorSetOwned {
            set: match view {
                SparseDenseColorSetView::Dense(bv) => {
                    SetType::Dense(bv.to_bitvec())
                },
                SparseDenseColorSetView::Sparse(iv_slice) => {
                    let len = iv_slice.end - iv_slice.start;
                    let mut new_iv = IntVector::with_capacity(len, iv_slice.vec.width()).unwrap();
                    for i in iv_slice.start..iv_slice.end {
                        new_iv.push(iv_slice.vec.get(i));
                    }
                    SetType::Sparse(new_iv)
                },
            }
        }
    }
    
    fn owned_to_view<'a>(&self, owned: &'a Self::OwnedSet) -> Self::SetView<'a> {
        match &owned.set {
            SetType::Dense(bv) => {
                SparseDenseColorSetView::Dense(bv.as_bitslice())
            },
            SetType::Sparse(iv) => {
                SparseDenseColorSetView::Sparse(IntVecSlice{
                    vec: iv,
                    start: 0,
                    end: iv.len(),
                })
            },
        }
    }
    
    // This is where we need the assumption that the IntVecs are sorted.
    // After this call, `a` has the intersection of `a` and `b`, in sorted
    // order when it's in IntVec form.
    fn intersect(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
        match (&mut a.set, b) {
            (SetType::Dense(a_bv), SparseDenseColorSetView::Dense(b_bv)) => {
                // Both dense -> bitwise AND
                *a_bv &= *b_bv;
                // Todo: re-encode into sparse if the intersection is small? 
            },
            (SetType::Dense(a_bv), SparseDenseColorSetView::Sparse(b_iv_slice)) => {
                // Intersection of Sparse and Dense will be Sparse   
                let s = b_iv_slice.start;
                let e = b_iv_slice.end;
                let v = &b_iv_slice.vec;
                let mut new_elements = IntVector::with_capacity(e-s, v.width()).unwrap();
                for v_idx in s..e {
                    let x = v.get(v_idx) as usize;
                    if a_bv[x] {
                        new_elements.push(x as u64);
                    }
                }
                *a = SparseDenseColorSetOwned { set: SetType::Sparse(new_elements) };
            },
            (SetType::Sparse(a_iv), SparseDenseColorSetView::Dense(b_bv)) => {
                // Intersection of Sparse and Dense will be Sparse.
                // Remove the elements of a that do not occur in b.
                // We do this in-place.
                let mut new_set_end = 0_usize;
                for i in 0..a_iv.len() {
                    if b_bv[a_iv.get(i) as usize] {
                        a_iv.set(new_set_end, a_iv.get(i));
                        new_set_end += 1; 
                    }
                }
                a_iv.resize(new_set_end, 0);
            },
            (SetType::Sparse(a_iv), SparseDenseColorSetView::Sparse(b_iv_slice)) => {
                // Intersection of sparse and sparse will be sparse.
                // We assume that the sets are sorted.
                // Let's do a scanning in-place intersection into a.
                let mut new_set_end = 0_usize;
                let mut b_pos = b_iv_slice.start;
                for i in 0..a_iv.len() {
                    let a_elem = a_iv.get(i) as usize;
                    while b_pos < b_iv_slice.end && (b_iv_slice.vec.get(b_pos) as usize) < a_elem {
                        b_pos += 1;
                    }
                    if b_pos == b_iv_slice.end {
                        break; // No more elements in b
                    }
                    if b_iv_slice.vec.get(b_pos) as usize == a_elem {
                        // Element is in both sets -> keep it
                        a_iv.set(new_set_end, a_iv.get(i));
                        new_set_end += 1;
                    }
                }
                a_iv.resize(new_set_end, 0);
            },
        }
    }
    
    // This should maintain the sorted order of a
    fn union(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
        // Really slow dummy implementation
        let mut elements: Vec<usize> = a.iter().chain(b.iter()).collect();
        elements.sort_unstable();
        elements.dedup();
        *a = SparseDenseColorSetOwned::new(elements.into_iter(), self.n_colors);
    }

     

}

/*
 *
 * Implementation of non-trait functions
 *
 */

impl SparseDenseStorage {

    pub fn get(&self, id: usize) -> SparseDenseColorSetView<'_> {
        if self.is_dense_marks.get(id) {
            let set_idx = self.is_dense_marks.rank(id);
            SparseDenseColorSetView::Dense(self.dense_sets.get(set_idx))
        } else {
            let set_idx = self.is_dense_marks.rank_zero(id);
            SparseDenseColorSetView::Sparse(self.sparse_sets.get(set_idx))
        }
    }
}

impl SortedIntVecs {
    fn new(bit_width: usize) -> Self {
        SortedIntVecs{concat: IntVector::new(bit_width).unwrap(), ends: vec![0]}
    }

    fn push(&mut self, set: impl IntoIterator<Item = usize>) { // Pushes a new set of integers
        for x in set {
            self.concat.push(x as u64);
        }
        self.ends.push(self.concat.len());
    }

    fn shrink_to_fit(&mut self) {
        self.concat.resize(self.concat.len(), 0);
    }

    fn get(&self, vec_idx: usize) -> IntVecSlice<'_> {
        IntVecSlice{vec: &self.concat, start: self.ends[vec_idx], end: self.ends[vec_idx+1]}
    }

    fn n_sets(&self) -> usize {
        self.ends.len() - 1 // Minus 1 because there is a 0 at the start of ends
    }

    fn serialize(&self, out: &mut impl std::io::Write) {
        // Serialize using bincode
        self.concat.serialize(out).unwrap();
        bincode::serialize_into(out, &self.ends).unwrap();
    }

    fn load(input: &mut impl std::io::Read) -> Self {
        // Deserialize using bincode
        let intvec_data = IntVector::load(input).unwrap();
        let ends: Vec<usize> = bincode::deserialize_from(input).unwrap();
        assert!(!ends.is_empty() && ends[0] == 0); // The first end must be 0
        SortedIntVecs{concat: intvec_data, ends}
    }

}

impl BitMaps {
    fn new(individual_length: usize) -> Self {
        BitMaps{bitmap_data: bitvec::vec::BitVec::new(), individual_length}
    }

    fn push(&mut self, bv: &bitvec::slice::BitSlice) {
        assert_eq!(bv.len(), self.individual_length);
        self.bitmap_data.extend_from_bitslice(bv);
    }

    fn shrink_to_fit(&mut self) {
        self.bitmap_data.shrink_to_fit();
    }

    fn get(&self, bitmap_idx: usize) -> &BitSlice {
        &self.bitmap_data[bitmap_idx*self.individual_length .. (bitmap_idx + 1) * self.individual_length]
    }

    #[allow(dead_code)]
    fn n_sets(&self) -> usize {
        self.bitmap_data.len() / self.individual_length
    }

    pub fn serialize(&self, out: &mut impl std::io::Write) {
        // Serialize using bincode
        bincode::serialize_into(out.by_ref(), &self.bitmap_data).unwrap();
        bincode::serialize_into(out.by_ref(), &self.individual_length).unwrap();
    }

    pub fn load(input: &mut impl std::io::Read) -> Self {
        // Deserialize using bincode
        let bitmap_data: bitvec::vec::BitVec = bincode::deserialize_from(input.by_ref()).unwrap();
        let individual_length: usize = bincode::deserialize_from(input.by_ref()).unwrap();
        assert!(individual_length > 0);
        BitMaps{bitmap_data, individual_length}
    }
}

fn is_dense_set(n_elements: usize, bits_per_color: usize, n_colors: usize) -> bool {
    let intvec_size = n_elements * bits_per_color;
    let bitmap_size = n_colors;
    bitmap_size <= intvec_size
}

#[cfg(test)]
mod tests {
    use crate::coloring_interface::*;
    use super::*;

    #[test]
    fn test_sparse_dense_storage() {
        let n_colors = 1000;
        let sets = [
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10], // Sparse
            (100..900).collect(), // Dense
            (0..n_colors).collect(), // Dense
            vec![0, 400, 700, 999], // Sparse
            vec![], // Sparse
            vec![42], // Sparse
        ];

        let should_be_sparse = [
            true,
            false,
            false,
            true,
            true,
            true,
        ];

        let iter_of_iter = sets.iter().map(|s| s.iter());
        let storage = SparseDenseStorage::new(iter_of_iter, n_colors);

        // Serialize and load
        let mut buf: Vec<u8> = vec![];
        storage.serialize(&mut buf);
        let storage = SparseDenseStorage::load(&mut buf.as_slice());

        // Check that we can retrieve the sets correctly
        for (i, true_set) in sets.iter().enumerate() {
            let view = storage.get(i);
            let owned = storage.view_to_owned(&view);
            let view_of_owned = storage.owned_to_view(&owned);

            eprintln!("{:?}", view.iter().collect::<Vec<usize>>());
            eprintln!("{:?}", owned.iter().collect::<Vec<usize>>());
            eprintln!("{:?}", view_of_owned.iter().collect::<Vec<usize>>());

            if should_be_sparse[i] {
                assert!(matches!(view, SparseDenseColorSetView::Sparse(_)));
                assert!(matches!(owned.set, SetType::Sparse(_)));
                assert!(matches!(view_of_owned, SparseDenseColorSetView::Sparse(_)));
            } else {
                assert!(matches!(view, SparseDenseColorSetView::Dense(_)));
                assert!(matches!(owned.set, SetType::Dense(_)));
                assert!(matches!(view_of_owned, SparseDenseColorSetView::Dense(_)));
            }

            assert_eq!(view.iter().collect::<Vec::<usize>>(), *true_set);
            assert_eq!(owned.iter().collect::<Vec::<usize>>(), *true_set);
            assert_eq!(view_of_owned.iter().collect::<Vec::<usize>>(), *true_set);
        }

        // Check all pairwise unions and intersections
        for i in 0..sets.len() {
            for j in 0..sets.len() {

                // Compute true union 
                let mut true_union: Vec<usize> = vec![];
                true_union.extend(sets[i].iter());
                true_union.extend(sets[j].iter());
                true_union.sort();
                true_union.dedup();

                // Compute true intersection
                let mut i_hash: HashSet<usize> = HashSet::new();
                i_hash.extend(sets[i].iter());
                let mut true_intersection: Vec<usize> = vec![];
                for x in sets[j].iter() {
                    if i_hash.contains(&x) {
                        true_intersection.push(x);
                    }
                }
                true_intersection.sort();

                // Compare to our union
                let mut owned_i = storage.view_to_owned(&storage.get(i));
                let view_j = storage.get(j);
                storage.union(&mut owned_i, &view_j);
                let mut our_union = owned_i.iter().collect::<Vec::<usize>>();
                our_union.sort();

                assert_eq!(our_union, true_union);

                // Compare to our intersection
                let mut owned_i = storage.view_to_owned(&storage.get(i));
                let view_j = storage.get(j);
                storage.intersect(&mut owned_i, &view_j);
                let mut our_intersection = owned_i.iter().collect::<Vec::<usize>>();
                our_intersection.sort();
                assert_eq!(our_intersection, true_intersection);
            }
        }
    }
}