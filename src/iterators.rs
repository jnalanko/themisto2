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


#[cfg(test)]
mod tests {

    struct VecVecStorage {
        sets: Vec<Vec<usize>>,
        pos: usize,
    }

    struct MyIterator<'a> {
        inner: std::slice::Iter<'a, usize>,
    }

    impl<'a> USizeIterator<'a> for MyIterator<'a> {
        fn next(&mut self) -> Option<usize> {
            self.inner.next().copied()
        }
    }

    impl USizeIteratorGenerator for VecVecStorage {
        type Iter<'a> = MyIterator<'a>;
        
        fn next<'a>(&'a mut self) -> Option<Self::Iter<'a>> {
            if self.pos == self.sets.len() {
                None
            } else {
                let iter = MyIterator{inner: self.sets[self.pos].iter()};
                self.pos += 1;
                Some(iter)
            }
        }
    }

    use super::*;

    #[test]
    fn test_iterator_traits() {
        let mut storage = VecVecStorage{sets: vec![vec![1,2,3], vec![2,3,4], vec![3,4,5]], pos: 0};
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
