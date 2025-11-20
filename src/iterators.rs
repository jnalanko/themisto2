trait MyIteratorTrait<'a> {
    fn next(&mut self) -> Option<usize>;
}

struct MyIterator<'a> {
    inner: std::slice::Iter<'a, usize>,

    owned: Vec<usize>,
    owned_pos: usize,    
}

impl<'a> MyIteratorTrait<'a> for MyIterator<'a> {
    fn next(&mut self) -> Option<usize> {
        self.inner.next().copied()

        /*
        if self.owned_pos == self.owned.len() { None }
        else {
            let x = self.owned[self.owned_pos];
            self.owned_pos += 1;
            Some(x)
        }
        */
    }
}

trait MyIteratorGeneratorTrait {
    type Iter<'a>: MyIteratorTrait<'a> where Self: 'a;

    fn next<'a>(&'a mut self) -> Option<Self::Iter<'a>>;
}

struct VecVecStorage {
    sets: Vec<Vec<usize>>,
    pos: usize,
}

impl MyIteratorGeneratorTrait for VecVecStorage {
    type Iter<'a> = MyIterator<'a>;
    
    fn next<'a>(&'a mut self) -> Option<Self::Iter<'a>> {
        if self.pos == self.sets.len() {
            None
        } else {
            let iter = MyIterator{inner: self.sets[self.pos].iter(), owned: self.sets[self.pos].clone(), owned_pos: 0};
            self.pos += 1;
            Some(iter)
        }
    }
}

fn main() {
    let mut storage = VecVecStorage{sets: vec![vec![1,2,3], vec![2,3,4], vec![3,4,5]], pos: 0};
    while let Some(mut iter) = storage.next() {
        while let Some(x) = iter.next() {
            println!("{}", x); 
        }
        println!("=="); 
    }

}
