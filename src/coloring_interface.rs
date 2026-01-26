// This file describes an interface to a color set storage.

// There are two traits for color sets: one where the underlying data is borrows, and one
// where it is owned. They are separate because abstracting over ownership in Rust is 
// currently super complicated and full of landmines and depends on obscure details of generic associated
// types. Don't do it. See: https://lucumr.pocoo.org/2022/9/11/abstracting-over-ownership/

use crate::iterators::{USizeIteratorGenerator, VecIterator};

// This trait represents a read-only storage struct that stores many color sets.
// The sets are viewed through returned structs implementing the associated color set
// view class. 
pub trait ColorSetStorage {

    // A generic associated color set view type. We could have e.g.
    // ColorSetStorage<BitMapColorSet<'a>> and ColorSetStorage<VecColorSet<'a>>.
    type SetView<'a>: ColorSetView<'a> + Clone + Send + Sync + Eq + std::hash::Hash + std::fmt::Debug where Self: 'a;

    // An owned version of SetView
    type OwnedSet: ColorSetOwned + Clone + Send + Sync;

    // Gives a set with a lifetime linked to the lifetime of the &self borrow.
    fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow>;

    // A fully generic way to construct the storage from a stream of color sets.
    // Returning box so that I can provide the method new_from_iter_of_iters
    // because there the size of the return value must be known at compile time.
    fn new(sets: impl USizeIteratorGenerator, n_colors: usize) -> Box<Self>; // Todo: rename to match better with `new_from_transpose`?

    // Construct in parallel from a generator of set_of_sets_construction::SetElement
    fn new_parallel(element_gen: impl crate::set_of_sets_construction::ParallelElementGenerator, n_colors: usize, set_sizes: &[usize], n_threads: usize) -> Box<Self>;

    // Construct in parallel in chunks, one element generator for each, and write the data
    // directly to disk in chunks.
    fn new_parallel_to_disk(element_gens: Vec<(impl crate::set_of_sets_construction::ParallelElementGenerator, std::ops::Range<usize>)>, set_sizes: Vec<usize>, output_prefix: &std::path::Path, n_threads: usize);

    fn n_sets(&self) -> usize;
    fn n_colors(&self) -> usize;

    fn get_empty_set(&self) -> Self::OwnedSet;
    fn get_full_set(&self) -> Self::OwnedSet;

    fn serialize<W: std::io::Write>(&self, out: &mut W);
    fn load<R: std::io::Read>(input: &mut R) -> Self;

    // Functions to convert between views and owned sets.
    // One would think that these should be methods of SetView and Ownedset
    // called "to_owned" and "as_view". But those types do not know what is their
    // corresponding view or owned type. I did not want to add those as associated
    // types to the view and owned types because then they are not linked to the
    // associated types here, or to each other, so for example if we had a view,
    // and made it owned, and again a view, the type would have been 
    // Storage::View::Owned::View even though Storage::View and 
    // Storage::View::Owned::View are the same type, but the compiler does not
    // see that. So the solution is to put the conversion functions here at the
    // Storage trait, and now the types do not nest like that.
    #[allow(dead_code)] // There will be useful later
    fn view_to_owned(&self, view: &Self::SetView<'_>) -> Self::OwnedSet;
    #[allow(dead_code)] // These will be useful later
    fn owned_to_view<'a>(&self, owned: &'a Self::OwnedSet) -> Self::SetView<'a>;

    // Set intersection: a := a ∩ b
    fn intersect(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>);

    // Set union: a := a ∪ b
    fn union(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>);

    // Provided method: new from an iterator that gives iterators to color sets.
    // While iterating, stores the color sets in a reused local buffer.
    #[allow(dead_code)] // Could be useful later
    fn new_from_iter_of_iters<
        InnerIter: Iterator<Item = usize>, 
        OuterIter: Iterator<Item = InnerIter>>
    (it: OuterIter, n_colors: usize) -> Box::<Self> {
        let stream = ColorSetStreamFromIters{iters: it, cur_set: vec![]};
        Self::new(stream, n_colors)
    }
}

// A color set view that does not own the data, but can return an
// iterator into it. The lifetime 'a is not referred to in the methods here,
// but we need it so that implementors have a lifetime parameter to work with. 
// ColorSetStorage uses this 'a to link it to the lifetime of the storage.
pub trait ColorSetView<'a> {

    // This associated iterator type may have lifetime parameters even though they
    // are not listed here. The iterator must iterate the elements in increasing order!
    // This assumption is useful e.g. for fast set operations on the sets. Todo:
    // I should make Iter also a subtype of a marker trait that represents that the
    // iterator gives its values in order. This would prevent accidentally using
    // a non-sorted iterator here. 
    type Iter: Iterator<Item = usize>;

    // The returned iterator may have lifetime parameters even though they are 
    // not listed here. It is just a generic type that implements Iterator<usize>. 
    // If it does have lifetime parameters, they are *not* linked to the lifetime 
    // of the &self borrow, but implementors can choose to link them to 'a.
    // Consequences:
    // * The iterator may not borrow from &self because the compiler does not 
    //   understand that the lifetimes are related, so it would see that we
    //   are returning something with lifetime of &self, but that has no relation
    //   to 'a, so it can not do required the lifetime subtyping.
    // * The iterator may outlive the set because it does not borrow from the
    //   set itself but rather from some external structure with lifetime 'a.
    fn iter(&self) -> Self::Iter;
    
    fn len(&self) -> usize;
}

pub trait ColorSetOwned {
    type Iter<'a>: Iterator<Item = usize> where Self: 'a;

    // This is different from ColorSetView because here the borrow in the
    // iterator is tied to the &self borrow, allowing us to return values
    // that borrow from &self.
    fn iter(&self) -> Self::Iter<'_>; 

    fn len(&self) -> usize;
}

impl ColorSetOwned for Vec<usize> {

    //type Iter = std::vec::IntoIter<usize>;
    type Iter<'a> = std::iter::Copied<std::slice::Iter<'a, usize>>;

    fn iter(&self) -> Self::Iter<'_> {
        self.as_slice().iter().copied()
    }

    fn len(&self) -> usize {
        self.len()
    }

}

// Convenient struct to get a color set stream from an iterator
// of iterators.
struct ColorSetStreamFromIters<T : IterOfIters> {
    iters: T,
    cur_set: Vec<usize>, 
}

impl<T: IterOfIters> USizeIteratorGenerator for ColorSetStreamFromIters<T> {
    type Iter<'a> = VecIterator<'a> where Self: 'a;

    fn next<'b>(&'b mut self) -> Option<Self::Iter<'b>> {
        if let Some(set_iter) = self.iters.next() {
            self.cur_set.clear();
            self.cur_set.extend(set_iter);
            Some(VecIterator::new(&self.cur_set))
        } else {
            None
        }
    }
    
}


pub trait IterOfIters {
    type Iter: Iterator<Item = usize>;
    fn next(&mut self) -> Option<Self::Iter>;
}

impl <InnerIter: Iterator<Item = usize>, OuterIter: Iterator<Item = InnerIter>> IterOfIters for OuterIter {
    type Iter = InnerIter;

    fn next(&mut self) -> Option<Self::Iter> {
        self.next()
    }
}



#[cfg(test)]
mod tests {

    use crate::iterators::USizeIterator;

    use super::*;

    pub struct ColorSetStorageVec {
        v: Vec<Vec<usize>>,
    }

    impl ColorSetStorage for ColorSetStorageVec {
        type SetView<'a> = SliceColorSet<'a> where Self: 'a;
        type OwnedSet = Vec<usize>;

        fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow> {
            SliceColorSet { // Dummy implementation
                slice: &self.v[id],
            }
        }
        
        fn n_sets(&self) -> usize {
            self.v.len()
        }

        fn n_colors(&self) -> usize {
            unimplemented!();
        }
        
        fn get_empty_set(&self) -> Self::OwnedSet {
            vec![]
        }
        
        fn get_full_set(&self) -> Self::OwnedSet {
            todo!()
        }

        
        #[allow(unused_variables)] // It's a dummy anyway
        fn serialize<W: std::io::Write>(&self, out: &mut W) {
            todo!()
        }
        
        #[allow(unused_variables)] // It's a dummy anyway
        fn load<R: std::io::Read>(input: &mut R) -> Self {
            todo!()
        }

        fn view_to_owned(&self, view: &Self::SetView<'_>) -> Self::OwnedSet {
            view.slice.to_vec()
        } 

        fn owned_to_view<'a>(&self, owned: &'a Self::OwnedSet) -> Self::SetView<'a> {
            SliceColorSet { // Dummy implementation
                slice: owned.as_slice()
            }
        }
        
        #[allow(unused_variables)] // It's a dummy anyway
        fn intersect(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
            todo!()
        }
        
        #[allow(unused_variables)] // It's a dummy anyway
        fn union(&self, a: &mut Self::OwnedSet, b: &Self::SetView<'_>) {
            todo!()
        }
        
        #[allow(unused_variables)] // It's a dummy anyway
        fn new_parallel(element_gen: impl crate::set_of_sets_construction::ParallelElementGenerator, n_colors: usize, set_sizes: &[usize], n_threads: usize) -> Box<Self> {
            todo!()
        }
        
        fn new(mut sets: impl USizeIteratorGenerator, _n_colors: usize) -> Box<Self> {
            while let Some(mut set) = sets.next() {
                while let Some(elem) = set.next() {
                    println!("{}", elem);
                }
            }
            todo!();
        }

    }


    #[derive(Debug, Clone, Hash, Eq, PartialEq)]
    pub struct SliceColorSet<'storage> {
        slice: &'storage [usize],
    }

    impl<'storage> ColorSetView<'storage> for SliceColorSet<'storage> {

        // The iterator type depends on the same lifetime as the ColorSet
        type Iter = SliceColorSetIter<'storage>;

        // The lifetime in the returned iter is NOT linked to the lifetime of the
        // &self borrow. So it is allowed to last longer than the borrow and in fact
        // even longer than Self.
        fn iter(&self) -> Self::Iter {
            SliceColorSetIter {
                slice: self.slice,
                pos: 0,
            }
        }
        
        fn len(&self) -> usize {
            todo!()
        }
    }

    #[derive(Debug)]
    pub struct SliceColorSetIter<'storage> {
        slice: &'storage [usize],
        pos: usize,
    }

    // We can keep the lifetime parameter anonymous because we do
    // not return references color color ids, but owned color ids instead.
    impl Iterator for SliceColorSetIter<'_> {
        type Item = usize;

        fn next(&mut self) -> Option<Self::Item> {
            if self.pos >= self.slice.len() {
                None
            } else {
                let item = self.slice[self.pos];
                self.pos += 1;
                Some(item)
            }
        }
    }

    fn generic<CSS: ColorSetStorage>(storage: CSS, true_sets: &[Vec<usize>]) {

        // Calling get_set returns a color set with the same lifetime as the borrow
        let set0 = storage.get_set_view(0);
        let iter0 = set0.iter();
        let set1 = storage.get_set_view(1);
        let iter1 = set1.iter();

        // Can drop a set but still keep and use the iterator since
        // the iterator does not borrow from the set, but depends on the
        // lifetime of the storage instead
        drop(set0);
        iter0.for_each(|x| println!("{}", x));
        iter1.for_each(|x| println!("{}", x));

        let owned = storage.get_empty_set();
        owned.iter().for_each(|x| println!("{}", x));
        owned.iter().for_each(|x| println!("{}", x)); // Can iterate multiple times

        for id in 0..true_sets.len() {
            let view = storage.get_set_view(id);
            assert_eq!(view.iter().collect::<Vec::<usize>>(), true_sets[id]);
            let owned = storage.view_to_owned(&view);
            assert_eq!(owned.iter().collect::<Vec::<usize>>(), true_sets[id]);
            let owned_view = storage.owned_to_view(&owned);
            assert_eq!(owned_view.iter().collect::<Vec::<usize>>(), true_sets[id]);
        }


    }

    #[test]
    fn color_set_traits() {
        let true_sets = vec![vec![1,2,3], vec![4,5,6], vec![7,8,9]];
        let storage = ColorSetStorageVec { v: true_sets.clone()};
        generic(storage, &true_sets);
    }

}