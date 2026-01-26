use std::ops::Range;
use std::path::Path;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering::{Acquire, Release, SeqCst};

use simple_sds_sbwt::raw_vector::AccessRaw;
use simple_sds_sbwt::serialize::Serialize;
use simple_sds_sbwt::{ops::{BitVec, Rank}, raw_vector::PushRaw};
use bitvec::slice::BitSlice;
use bitvec::bitvec;

use crate::atomic_bitmap::AtomicBitmap;
use crate::int_vec::{AtomicCompactIntVec, CompactIntVec};
use crate::coloring_interface::{ColorSetOwned, ColorSetView};
use crate::iterators::{USizeIterator, USizeIteratorGenerator};
use crate::set_of_sets_construction;

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
#[derive(Clone, Debug)]
struct SortedIntVecs {
    // Concatenation of IntVecs. Better keep this private because the sets need
    // to maintained in sorted order for intersections. Since it's private, we only need 
    // to make sure that code in this module respects that. Users can get
    // immutable views to this via IntVecSlice.
    concat: CompactIntVec, 

    // Starts of individual sets, with a final sentinel value at the exclusive end of the last one.
    starts: Vec<usize>, 
}

#[derive(Copy, Clone, Debug)]
pub struct IntVecSlice<'a> {
    vec: &'a CompactIntVec,
    start: usize,
    end: usize, // Exclusive end
}

// A set of sets encoded as bitmaps.
#[derive(Clone, Debug)]
struct BitMaps {
    bitmap_data: bitvec::vec::BitVec, // Concatenation of bit vectors
    individual_length: usize, // Length of each bitmap in bitmap_data
}

#[derive(Clone)]
pub enum OwnedSparseDense {
    Dense(bitvec::vec::BitVec),
    Sparse(CompactIntVec),
}

// This enum is only for passing references to individual sets around. The actual
// sets are stored in concatenated form somewhere else in memory. 
#[derive(Copy, Clone, Hash, Eq, PartialEq, Debug)]
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

impl PartialEq for IntVecSlice<'_> {
    fn eq(&self, other: &Self) -> bool {
        if self.end - self.start != other.end - other.start {
            return false; // Different length
        }

        // Todo: could do word-by-word comparison
        let len = self.end - self.start;
        (0..len).into_iter().all(|i| self.vec.get(self.start + i) == other.vec.get(other.start + i))
    }
}
impl Eq for IntVecSlice<'_> {}

impl std::hash::Hash for IntVecSlice<'_> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        for i in self.start..self.end {
            self.vec.get(i).hash(state);
        }
    }
}

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
                    Some(x)
                }
            },
        }
    }
}

impl crate::coloring_interface::ColorSetOwned for OwnedSparseDense {
    type Iter<'a> = ColorSetViewIterator<'a> where Self: 'a;

    fn iter(&self) -> Self::Iter<'_> {
        match &self {
            OwnedSparseDense::Dense(bv) => {
                ColorSetViewIterator{
                    set: SparseDenseColorSetView::Dense(bv.as_bitslice()),
                    pos: 0,
                }
            },
            OwnedSparseDense::Sparse(iv) => {
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
    
    fn len(&self) -> usize {
        match &self {
            OwnedSparseDense::Dense(bit_vec) => bit_vec.count_ones(),
            OwnedSparseDense::Sparse(compact_int_vec) => compact_int_vec.len(),
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
    type OwnedSet = OwnedSparseDense;
    
    fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow> {
        self.get(id)
    }

    fn get_empty_set(&self) -> Self::OwnedSet {
        let bit_width = self.n_colors.next_power_of_two().trailing_zeros() as usize;
        OwnedSparseDense::Sparse(CompactIntVec::new(0, bit_width))
    }

    fn get_full_set(&self) -> Self::OwnedSet {
        OwnedSparseDense::Dense(bitvec![1; self.n_colors])
    }
    
    fn new(mut sets: impl USizeIteratorGenerator, n_colors: usize) -> Box<Self> {
        log::info!("Encoding color sets");
        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;
        let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

        // Since CompactIntVec does not support push, we need to simulate it here.
        let mut sparse_concat_capacity = 1; // Current capacity of sparse_concat
        let mut sparse_concat_len = 0; // Current len of sparse_concat
        let mut sparse_concat = CompactIntVec::new(sparse_concat_capacity, color_id_bit_width);
        let mut sparse_concat_starts: Vec<usize> = vec![];

        let mut dense_sets = BitMaps::new(n_colors);

        let mut buf = Vec::<usize>::new();
        let mut n_sets_total = 0_usize;
        while let Some(mut set) = sets.next() {
            buf.clear();
            while let Some(x) = set.next() {
                buf.push(x);
            }
            if is_dense_formula(buf.len(), color_id_bit_width, n_colors) {
                let mut bm = bitvec![0; n_colors];
                for color in buf.iter() {
                    bm.set(color, true);
                }
                dense_sets.push(&bm);
                is_dense_marks.push_bit(true);
            } else {
                buf.sort_unstable(); // We need sorted sets for intersections

                // Push the sorted set to the sparse concat
                sparse_concat_starts.push(sparse_concat_len);
                for x in buf.iter() {
                    if sparse_concat_len == sparse_concat_capacity {
                        // Make more space
                        sparse_concat_capacity *= 2;
                        sparse_concat.resize(sparse_concat_capacity);
                    }
                    sparse_concat.set(sparse_concat_len, x);
                    sparse_concat_len += 1;
                }

                is_dense_marks.push_bit(false);
            }

            n_sets_total += 1;
        }

        sparse_concat_starts.push(sparse_concat_len); // Final sentinel start, as required

        sparse_concat.resize(sparse_concat_len); // Shrink to fit
        sparse_concat_starts.shrink_to_fit();
        dense_sets.shrink_to_fit();

        let sparse_sets = SortedIntVecs {
            concat: sparse_concat,
            starts: sparse_concat_starts,
        };

        log::info!("{}% of the sets are sparse", sparse_sets.n_sets() as f64 / n_sets_total as f64 * 100.0);

        // Add rank support to dense marks
        log::info!("Building rank support for dense marks");
          let mut is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::from(is_dense_marks);
        is_dense_marks.enable_rank();

        Box::new(SparseDenseStorage {
            is_dense_marks, 
            sparse_sets,
            dense_sets,
            n_colors
        })
    }

    fn new_parallel(mut element_gen: impl crate::set_of_sets_construction::ParallelElementGenerator, n_colors: usize, set_sizes: &[usize], n_threads: usize) -> Box<Self> {

        // TODO This function has a ton of repetition. Should clean this up

        log::info!("Encoding color sets");

        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;
        let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

        let mut n_dense_sets = 0;
        let mut n_sparse_sets = 0;
        let mut total_sparse_size = 0;
        for &size in set_sizes.iter() {
            if is_dense_formula(size, color_id_bit_width, n_colors) {
                is_dense_marks.push_bit(true);
                n_dense_sets += 1;
            } else {
                is_dense_marks.push_bit(false);
                n_sparse_sets += 1;
                total_sparse_size += size;
            }
        }
        // Todo: shrink to fit dense marks

        dbg!(&n_sparse_sets, &n_dense_sets, &color_id_bit_width, &n_colors);

        // Zero-initialized data with atomic updates 
        let sparse_set_concat = AtomicCompactIntVec::new(total_sparse_size, color_id_bit_width);
        let dense_sets = AtomicBitmap::new(n_colors * n_dense_sets);

        let sparse_set_insertion_points: Vec<AtomicU64> = (0..n_sparse_sets).map(|_| std::sync::atomic::AtomicU64::new(0)).collect();
        let mut total_sparse_set_size = 0_usize;
        let mut ugly_sparse_id_again = 0_usize;
        for (i, s) in set_sizes.iter().enumerate() {
            if !is_dense_marks.bit(i) {
                sparse_set_insertion_points[ugly_sparse_id_again].store(total_sparse_set_size as u64, Release);
                total_sparse_set_size += s;
                ugly_sparse_id_again += 1;
            }
        }

        log::info!("Building rank support for dense marks");
        let mut is_dense_marks = simple_sds_sbwt::bit_vector::BitVector::from(is_dense_marks);
        is_dense_marks.enable_rank();

        // Fill in the set elements to the zero-initialized data
        log::info!("Encoding sets");
        let callback = |element: crate::set_of_sets_construction::SetElement| {
            let set_id = element.set_id;
            let color_id = element.color;
            if is_dense_marks.get(set_id) {
                let dense_id = is_dense_marks.rank(set_id);
                dense_sets.set(dense_id * n_colors + color_id, true);
            } else {
                let sparse_id = is_dense_marks.rank_zero(set_id);
                
                // Increment the insertion point. If I understand correctly, fetch_update will
                // fetch the value, apply a function to it, and try to compare-and-swap it in.
                // If the value has changed, it will fetch and try again until it succeeds.
                // When succesful, it will return the old value that was updated.
                // TODO: fetch_add instead?

                let insertion_point = sparse_set_insertion_points[sparse_id].fetch_update(SeqCst, SeqCst, |x| Some(x+1)).unwrap();

                // Hooray, we have succeeded in incrementing the insertion point. This means
                // that no other thread can write at the point, so we can safely write it. 

                sparse_set_concat.set(insertion_point as usize, color_id);
            }
        };
        element_gen.run(callback, n_threads);

        // Atomicity no longer needed
        let sparse_set_concat = CompactIntVec::from_atomic(sparse_set_concat);
        let mut sparse_set_starts: Vec<usize> = sparse_set_insertion_points.into_iter().map(
            |x| x.load(Acquire) as usize
        ).collect();

        // After inserting everything, the insertion points are now that the starts of
        // the next set. Let's add the start of the first set to the beginning:
        sparse_set_starts.insert(0, 0);

        // Make a SortedIntVecs out of the sparse set concatenation
        sparse_set_starts.shrink_to_fit();
        let mut sorted_sparse_sets = SortedIntVecs {
            concat: sparse_set_concat,
            starts: sparse_set_starts,
        };

        log::info!("Sorting sparse sets");
        let mut sort_buf = Vec::<usize>::new();
        for sparse_id in 0..sorted_sparse_sets.n_sets() {
            sort_buf.clear();

            let piece = sorted_sparse_sets.get(sparse_id);
            for i in piece.start..piece.end {
                sort_buf.push(piece.vec.get(i));
            }
            sort_buf.sort();
            for offset in 0..sort_buf.len() {
                sorted_sparse_sets.assign_element(sparse_id, offset, sort_buf[offset]);
            }
        }
        log::info!("Sparse sets sorted");

        Box::new(Self {
            dense_sets: BitMaps { 
                bitmap_data: dense_sets.into_bitvec(), 
                individual_length: n_colors
            },
            sparse_sets: sorted_sparse_sets,
            n_colors,
            is_dense_marks,
        })
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
        assert!(sparse_sets.concat.bit_width() >= n_colors.next_power_of_two().trailing_zeros() as usize);

        Self {is_dense_marks, sparse_sets, dense_sets, n_colors}
    }
    
    fn view_to_owned(&self, view: &Self::SetView<'_>) -> Self::OwnedSet {
        match view {
            SparseDenseColorSetView::Dense(bv) => {
                OwnedSparseDense::Dense(bv.to_bitvec())
            },
            SparseDenseColorSetView::Sparse(iv_slice) => {
                let len = iv_slice.end - iv_slice.start;
                let mut new_iv = CompactIntVec::new(len, iv_slice.vec.bit_width());
                let mut push_idx = 0_usize;
                #[allow(clippy::explicit_counter_loop)] // I don't want to touch it
                for i in iv_slice.start..iv_slice.end {
                    new_iv.set(push_idx, iv_slice.vec.get(i));
                    push_idx += 1;
                }
                OwnedSparseDense::Sparse(new_iv)
            },
        }
    }
    
    fn owned_to_view<'a>(&self, owned: &'a Self::OwnedSet) -> Self::SetView<'a> {
        match &owned {
            OwnedSparseDense::Dense(bv) => {
                SparseDenseColorSetView::Dense(bv.as_bitslice())
            },
            OwnedSparseDense::Sparse(iv) => {
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
    fn intersect(&self, mut a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
        match (&mut a, b) {
            (OwnedSparseDense::Dense(a_bv), SparseDenseColorSetView::Dense(b_bv)) => {
                // Both dense -> bitwise AND
                *a_bv &= *b_bv;

                // Check if we need to re-encode as sparse
                if self.is_dense(a_bv.count_ones()) {
                    // Re-encode
                    *a = self.new_owned_set(a_bv.iter_ones());
                }
            },
            (OwnedSparseDense::Dense(a_bv), SparseDenseColorSetView::Sparse(b_iv_slice)) => {
                // Intersection of Sparse and Dense will be Sparse   
                let s = b_iv_slice.start;
                let e = b_iv_slice.end;
                let v = &b_iv_slice.vec;

                // Let's make space for the biggest possible intersection. We'll
                // truncante at the end.
                let mut new_elements = CompactIntVec::new(e-s, v.bit_width());
                let mut push_idx = 0_usize;
                for v_idx in s..e {
                    let x = v.get(v_idx);
                    if a_bv[x] {
                        new_elements.set(push_idx, x);
                        push_idx += 1;
                    }
                }
                new_elements.resize(push_idx); // Truncate to actual size
                *a = OwnedSparseDense::Sparse(new_elements);
            },
            (OwnedSparseDense::Sparse(a_iv), SparseDenseColorSetView::Dense(b_bv)) => {
                // Intersection of Sparse and Dense will be Sparse.
                // Remove the elements of a that do not occur in b.
                // We do this in-place.
                let mut new_set_end = 0_usize;
                for i in 0..a_iv.len() {
                    if b_bv[a_iv.get(i)] {
                        a_iv.set(new_set_end, a_iv.get(i));
                        new_set_end += 1; 
                    }
                }
                a_iv.resize(new_set_end); // Truncate the vector
            },
            (OwnedSparseDense::Sparse(a_iv), SparseDenseColorSetView::Sparse(b_iv_slice)) => {
                // Intersection of sparse and sparse will be sparse.
                // We assume that the sets are sorted.
                // Let's do a scanning in-place intersection into a.
                let mut new_set_end = 0_usize;
                let mut b_pos = b_iv_slice.start;
                for i in 0..a_iv.len() {
                    let a_elem = a_iv.get(i);
                    while b_pos < b_iv_slice.end && (b_iv_slice.vec.get(b_pos)) < a_elem {
                        b_pos += 1;
                    }
                    if b_pos == b_iv_slice.end {
                        break; // No more elements in b
                    }
                    if b_iv_slice.vec.get(b_pos) == a_elem {
                        // Element is in both sets -> keep it
                        a_iv.set(new_set_end, a_iv.get(i));
                        new_set_end += 1;
                    }
                }
                a_iv.resize(new_set_end); // Truncate the vector
            },
        }
    }
    
    // This should maintain the sorted order of a
    fn union(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
        // Really slow dummy implementation
        let mut elements: Vec<usize> = a.iter().chain(b.iter()).collect();
        elements.sort_unstable();
        elements.dedup();
        *a = self.new_owned_set(elements.into_iter());
    }
    
    fn n_sets(&self) -> usize {
        self.is_dense_marks.len()
    }
    
    fn n_colors(&self) -> usize {
        self.n_colors
    }
    
}

/*
 *
 * Implementation of non-trait functions
 *
 */

 // Returns true if a set with n_elements elements should be encoded as a sparse
 // set with integer bit width sparse_bit_width, or as a dense set supporting up
 // to n_colors.
fn is_dense_formula(n_elements: usize, sparse_bit_width: usize, n_colors: usize) -> bool {

    // sparse_bit_width bits per element, plus 64 bits to encode the offset in the concatenation
    let intvec_size = n_elements * sparse_bit_width + 64;

    // One bit per color
    let bitmap_size = n_colors;

    bitmap_size <= intvec_size
}

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

    fn sparse_bit_width(&self) -> usize {
        self.sparse_sets.concat.bit_width()
    }

    fn is_dense(&self, n_elements: usize) -> bool {
        is_dense_formula(n_elements, self.sparse_bit_width(), self.n_colors)
    }

    // n_colors-1 is the maximum color id supported by this set
    pub fn new_owned_set(&self, elements: impl Iterator<Item = usize>) -> OwnedSparseDense {
        let n_colors = self.n_colors;
        let elements: Vec<usize> = elements.collect();
        let bits_per_color = n_colors.next_power_of_two().trailing_zeros() as usize;
        if self.is_dense(elements.len()) {
            let mut bv = bitvec![0; n_colors];
            for color in elements.iter() {
                bv.set(color, true);
            }
            OwnedSparseDense::Dense(bv)
        } else {
            let mut push_idx = 0_usize;
            let mut iv = CompactIntVec::new(elements.len(), bits_per_color);
            #[allow(clippy::explicit_counter_loop)] // I don't want to touch it
            for color in elements.iter() {
                iv.set(push_idx, color);
                push_idx += 1;
            }
            OwnedSparseDense::Sparse(iv)
        }
    }

    // The Range<usize> associates to each element generator gives the color id range of that generator
    fn new_parallel_to_disk(mut element_gens: Vec<(impl crate::set_of_sets_construction::ParallelElementGenerator, Range<usize>)>, set_sizes: Vec<usize>, output_prefix: &Path, n_threads: usize) {
        log::info!("Encoding color sets to disk in {} pieces", element_gens.len());

        assert!(element_gens.len() > 0);
        assert!(element_gens.first().unwrap().1.start == 0);
        for i in 1..element_gens.len() {
            assert!(element_gens[i-1].1.end == element_gens[i].1.start);
        }
        let n_colors = element_gens.last().unwrap().1.end;

        let color_id_bit_width = n_colors.next_power_of_two().trailing_zeros() as usize;
        let mut is_dense_marks = simple_sds_sbwt::raw_vector::RawVector::new();

        let mut n_dense_sets = 0;
        let mut n_sparse_sets = 0;
        for size in set_sizes.iter() {
            if is_dense_formula(size, color_id_bit_width, n_colors) {
                is_dense_marks.push_bit(true);
                n_dense_sets += 1;
            } else { // Sparse
                is_dense_marks.push_bit(false);
                n_sparse_sets += 1;
            }
        }

        // Todo: shrink to fit dense marks
        let sparse_set_insertion_points: Vec<usize> = vec![0; n_sparse_sets];
        let mut n_sparse_sets_again = 0;
        let mut total_sparse_size = 0;
        for (i, size) in set_sizes.iter().enumerate() {
            if !is_dense_marks.bit(i) {
                sparse_set_insertion_points[n_sparse_sets_again] = total_sparse_size;
                total_sparse_size += size;
                n_sparse_sets_again += 1;
            } 
        }
        drop(set_sizes); // Free memory

        dbg!(&n_sparse_sets, &n_dense_sets, &color_id_bit_width, &n_colors);

        // Create three files:
        /*
            * A file for dense sets: `[data_n_words: usize][bitvec_len: usize][n_colors: usize][data]`
            * A file for sparse sets: `[data_n_words: usize][bit_width: usize][n_colors: usize][n_sets: usize][data][starts]`
            * A file for is_dense_marks: [data n_words: usize][bitvec_len: usize][data]
            (No rank support serialized for now)
        */
        let sparse_filename = crate::util::append_to_filename(output_prefix, ".sparse");
        let dense_filename = crate::util::append_to_filename(output_prefix, ".dense");
        let marks_filename = crate::util::append_to_filename(output_prefix, ".marks");

        let sparse_file = std::fs::File::create(sparse_filename).unwrap();
        let total_sparse_bits = 64 + 64 + 64 + 64 + total_sparse_size*color_id_bit_width + (n_sparse_sets + 1) * 64;
        sparse_file.set_len(total_sparse_bits.next_multiple_of(8) as u64).unwrap();

        let dense_file = std::fs::File::create(dense_filename).unwrap();
        let total_dense_bits: usize = 64 + 64 + 64 + n_dense_sets * n_colors;
        dense_file.set_len(total_dense_bits.next_multiple_of(8) as u64).unwrap();

        let marks_file = std::fs::File::create(marks_filename).unwrap();
        let total_marks_bits = 64 + 64 + is_dense_marks.len();
        marks_file.set_len(total_marks_bits.next_multiple_of(8) as u64).unwrap();

        for (element_gen, color_id_range) in element_gens {
            let slice_set_sizes: Vec<usize> = set_of_sets_construction::compute_set_sizes_assuming_no_duplicate_elements(&mut element_gen, is_dense_marks.len(), n_threads);
            element_gen.rewind();
            let piece = <Self as crate::coloring_interface::ColorSetStorage>::new_parallel(element_gen, n_colors, &slice_set_sizes, n_threads);
            Self::write_piece(piece, color_id_range, &mut sparse_file, &mut dense_file, sparse_write_pointers, n_threads);
        }
    }

}

impl SortedIntVecs {
    #[allow(dead_code)]
    fn new(bit_width: usize) -> Self {
        SortedIntVecs{concat: CompactIntVec::new(0, bit_width), starts: vec![0]}
    }

    #[allow(dead_code)]
    fn new_with_sizes(set_sizes: impl Iterator<Item = usize>, n_sets: usize, bit_width: usize) -> Self {
        let mut starts: Vec<usize> = vec![0; n_sets + 1];
        let mut total_set_size = 0_usize;
        for (i, s) in set_sizes.enumerate() {
            starts[i+1] = starts[i] + s;
            total_set_size += s;
        }

        let concat = CompactIntVec::new(total_set_size, bit_width);

        Self { concat, starts }
    }

    fn assign_element(&mut self, set_id: usize, offset_in_set: usize, element: usize) { // Pushes a new set of integers
        let size = self.starts[set_id+1] - self.starts[set_id];
        assert!(offset_in_set < size);
        self.concat.set(self.starts[set_id] + offset_in_set, element);
    }

    /*
    fn shrink_to_fit(&mut self) {
        self.concat.resize(self.concat.len(), 0);
    }
    */

    fn get(&self, set_id: usize) -> IntVecSlice<'_> {
        IntVecSlice{vec: &self.concat, start: self.starts[set_id], end: self.starts[set_id+1]}
    }

    fn n_sets(&self) -> usize {
        self.starts.len() - 1 // Minus 1 because there is a sentinel start at the end
    }

    fn serialize(&self, out: &mut impl std::io::Write) {
        // Serialize using bincode
        self.concat.serialize(out);
        bincode::serialize_into(out, &self.starts).unwrap();
    }

    fn load(input: &mut impl std::io::Read) -> Self {
        // Deserialize using bincode
        let concat = CompactIntVec::load(input);
        let starts: Vec<usize> = bincode::deserialize_from(input).unwrap();
        assert!(!starts.is_empty() && *starts.last().unwrap() == concat.len());
        SortedIntVecs{concat, starts}
    }

}

impl BitMaps {
    fn new(individual_length: usize) -> Self {
        BitMaps{bitmap_data: bitvec::vec::BitVec::new(), individual_length}
    }

    #[allow(dead_code)]
    fn new_with_zero_init(individual_length: usize, n_sets: usize) -> Self {
        BitMaps{bitmap_data: bitvec![0; n_sets*individual_length], individual_length}
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
    fn get_mut(&mut self, bitmap_idx: usize) -> &mut BitSlice {
        &mut self.bitmap_data[bitmap_idx*self.individual_length .. (bitmap_idx + 1) * self.individual_length]
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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use crate::{coloring_interface::*, iterators::VecVecUsizeIteratorGenerator};
    use super::*;

    #[test]
    fn test_sparse_dense_storage() {
        let n_colors = 1000;
        let sets = vec![
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

        // convert(...) turns an interator into a streaming iterator
        let iter_of_iter = VecVecUsizeIteratorGenerator::new(sets.clone());
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
                assert!(matches!(owned, OwnedSparseDense::Sparse(_)));
                assert!(matches!(view_of_owned, SparseDenseColorSetView::Sparse(_)));
            } else {
                assert!(matches!(view, SparseDenseColorSetView::Dense(_)));
                assert!(matches!(owned, OwnedSparseDense::Dense(_)));
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