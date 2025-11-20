// This file has some iterator-like traits that are used in this crate.
// We need to define our own iterator traits because the lifetime constraints
// in the standard iterator trait do not allow iterators creating iterators
// that stream over data from the parent iterator. At least I could not make it work.

// This is different from Iterator<Item = usize> because this has a lifetime
// parameter attached to it. This will be needed in the USizeIteratorGenerator trait
pub trait USizeIterator<'a> {
    fn next(&mut self) -> Option<usize>;
}

pub trait USizeIteratorGenerator {
    type Iter<'a>: USizeIterator<'a> where Self: 'a;

    // Here we link the lifetime in the iterator to the
    // lifetime of the self-borrow.
    fn next<'b>(&'b mut self) -> Option<Self::Iter<'b>>;
}

pub struct VecVecUsizeIteratorGenerator {
    pub(crate) sets: Vec<Vec<usize>>,
    pub(crate) pos: usize,
}

impl VecVecUsizeIteratorGenerator {
    pub fn new(vecs: Vec<Vec<usize>>) -> Self {
        Self {sets: vecs, pos: 0}
    }
}

pub struct VecIterator<'a> {
    inner: std::slice::Iter<'a, usize>,
}

impl<'a> VecIterator<'a> {
    pub fn new(vec: &'a Vec<usize>) -> Self {
        Self { inner: vec.as_slice().iter() }
    }
}

impl<'a> USizeIterator<'a> for VecIterator<'a> {
    fn next(&mut self) -> Option<usize> {
        self.inner.next().copied()
    }
}

impl USizeIteratorGenerator for VecVecUsizeIteratorGenerator {
    type Iter<'a> = VecIterator<'a>;
    
    fn next<'a>(&'a mut self) -> Option<Self::Iter<'a>> {
        if self.pos == self.sets.len() {
            None
        } else {
            let iter = VecIterator{inner: self.sets[self.pos].iter()};
            self.pos += 1;
            Some(iter)
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_iterator_traits() {
        let mut storage = VecVecUsizeIteratorGenerator{sets: vec![vec![1,2,3], vec![2,3,4], vec![3,4,5]], pos: 0};
        while let Some(mut iter) = storage.next() {
            while let Some(x) = iter.next() {
                println!("{}", x); 
            }
            println!("=="); 
        }
    }

    //fn print_all_generic<'a, 'b, Inner: USizeIterator<'a>, IterGenerator: USizeIteratorGenerator<Iter<'a> = Inner> + 'a>(mut gen: IterGenerator){
    fn print_all_generic(mut gen: impl USizeIteratorGenerator){
        while let Some(mut iter) = gen.next() {
            while let Some(x) = iter.next() {
                println!("{}", x); 
            }
            println!("=="); 
        }

    }
}
