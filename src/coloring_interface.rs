// There are two traits for color sets: one where the underlying data is borrows, and one
// where it is owned. They are separate because abstracting over ownership in Rust is 
// currently super complicated and full of landmines and depends on obscure details of generic associated
// types. Don't do it. See: https://lucumr.pocoo.org/2022/9/11/abstracting-over-ownership/

struct ColorSetStorageVec {
    v: Vec<usize>,
}

// This trait represents a read-only storage struct that stores many color sets.
// The sets are viewed through returned structs implementing the associated color set
// view class. 
pub trait ColorSetStorage {

    // A generic associated color set view type. We could have e.g.
    // ColorSetStorage<BitMapColorSet<'a>> and ColorSetStorage<VecColorSet<'a>>.
    type SetView<'a>: ColorSetView<'a> where Self: 'a;

    // An owned version of SetView
    type OwnedSet: ColorSetOwned;

    // Gives a set with a lifetime linked to the lifetime of the &self borrow.
    fn get_set_view<'borrow>(&'borrow self, id: usize) -> Self::SetView<'borrow>;
    fn get_owned_set(&self, id: usize) -> Self::OwnedSet;

    // Takes an iterator of iterators: Each inner iterator iterates the elements of one color set.
    // The color ids are in the range 0..n_colors.
    fn new(sets: impl Iterator<Item = impl Iterator<Item = usize>>, n_colors: usize) -> Self;

    fn get_empty_set(&self) -> Self::OwnedSet;
    fn get_full_set(&self) -> Self::OwnedSet;

    fn serialize<W: std::io::Write>(&self, out: &mut W);
    fn load<R: std::io::Read>(input: &mut R) -> Self;
}

// A color set view that does not own the data, but can return an
// iterator into it. The lifetime 'a is not referred to in the methods here,
// but we need it so that implementors have a lifetime parameter to work with. 
// ColorSetStorage uses this 'a to link it to the lifetime of the storage.
pub trait ColorSetView<'a> {

    // This associated iterator type may have lifetime parameters even though they
    // are not listed here.
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
    fn iter<'me>(&'me self) -> Self::Iter;
}

pub trait ColorSetOwned {
    type Iter<'a>: Iterator<Item = usize> where Self: 'a;

    fn intersect(&mut self, other: &impl ColorSetOwned);
    fn union(&mut self, other: &impl ColorSetOwned);

    // This is different from ColorSetView because here the borrow in the
    // iterator is tied to the &self borrow, allowing us to return values
    // that borrow from &self.
    fn iter(&self) -> Self::Iter<'_>; 
}

impl ColorSetOwned for Vec<usize> {

    //type Iter = std::vec::IntoIter<usize>;
    type Iter<'a> = std::iter::Copied<std::slice::Iter<'a, usize>>;

    fn intersect(&mut self, other: &impl ColorSetOwned) {
        todo!()
    }

    fn union(&mut self, other: &impl ColorSetOwned) {
        todo!()
    }

    fn iter(&self) -> Self::Iter<'_> {
        self.as_slice().iter().copied()
    }
}

impl ColorSetStorage for ColorSetStorageVec {
    type SetView<'a> = SliceColorSet<'a> where Self: 'a;
    type OwnedSet = Vec<usize>;

    fn get_set_view<'borrow>(&'borrow self, index: usize) -> Self::SetView<'borrow> {
        SliceColorSet { // Dummy implementation
            slice: &self.v[index..index + 5],
        }
    }
    
    fn get_owned_set(&self, id: usize) -> Self::OwnedSet {
        todo!()
    }
    
    fn get_empty_set(&self) -> Self::OwnedSet {
        vec![]
    }
    
    fn get_full_set(&self) -> Self::OwnedSet {
        todo!()
    }
    
    fn new(sets: impl Iterator<Item = impl Iterator<Item = usize>>, n_colors: usize) -> Self {
        for set in sets {
            for elem in set {
                println!("{}", elem);
            }
        }
        todo!();
    }
    
    fn serialize<W: std::io::Write>(&self, out: &mut W) {
        todo!()
    }
    
    fn load<R: std::io::Read>(input: &mut R) -> Self {
        todo!()
    }


}


#[derive(Debug)]
struct SliceColorSet<'storage> {
    slice: &'storage [usize],
}

impl<'storage> ColorSetView<'storage> for SliceColorSet<'storage> {

    // The iterator type depends on the same lifetime as the ColorSet
    type Iter = SliceColorSetIter<'storage>;

    // The lifetime in the returned iter is NOT linked to the lifetime of the
    // &self borrow. So it is allowed to last longer than the borrow and in fact
    // even longer than Self.
    fn iter<'me>(&'me self) -> Self::Iter {
        SliceColorSetIter {
            slice: self.slice,
            pos: 0,
        }
    }
}

#[derive(Debug)]
struct SliceColorSetIter<'storage> {
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

#[derive(Debug)]
struct OwnedColorSet {
    set: Vec<usize>
}

struct OwnedColorSetIter<'a> {
    pos: usize,
    set: &'a [usize],
}

impl<'a> Iterator for OwnedColorSetIter<'a> {
    type Item = usize;
    fn next(&mut self) -> Option<usize> {
        if self.pos >= self.set.len() {
            None
        } else {
            let item = self.set[self.pos];
            self.pos += 1;
            Some(item)
        }
    }
}

mod tests {

    use super::*;

    fn generic<CSS: ColorSetStorage>(storage: CSS) {
        // 'a is is the generic lifetime associated with color set objects from CSS 

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
        owned.iter().for_each(|x| println!("{}", x));

    }

    #[test]
    fn color_set_traits() {
        let storage = ColorSetStorageVec { v: vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]};
        generic(storage);
    }

}