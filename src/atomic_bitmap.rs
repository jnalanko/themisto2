use bitvec::vec::BitVec;

pub struct AtomicBitmap {
    pub data: Vec<std::sync::atomic::AtomicU64>,
    pub len: usize, // Number of bits stored. The last word may be only partially used.
} 

impl AtomicBitmap {
    pub fn new(len: usize) -> Self {
        let n_words = len.div_ceil(64);
        let data = (0..n_words).map(|_| std::sync::atomic::AtomicU64::new(0)).collect();
        Self {data, len}
    }

    pub fn set(&self, index: usize, value: bool) {
        assert!(index < self.len);
        let word_idx = index / 64;
        let bit_idx = index % 64;
        let mask = 1_u64 << bit_idx;
        if value {
            self.data[word_idx].fetch_or(mask, std::sync::atomic::Ordering::Release);
        } else {
            self.data[word_idx].fetch_and(!mask, std::sync::atomic::Ordering::Release);
        }
    }

    pub fn get(&self, index: usize) -> bool {
        assert!(index < self.len);
        let word_idx = index / 64;
        let bit_idx = index % 64;
        let mask = 1_u64 << bit_idx;
        (self.data[word_idx].load(std::sync::atomic::Ordering::Acquire) & mask) != 0
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn into_bitvec(self) -> BitVec<usize, bitvec::order::Lsb0> {
        // This assumes that usize is an u64 in Lsb0 byte order. Let's hope for the best.
        // At least the answer will be really wrong if that is not the case, so it
        // should be easy to notice. TODO: check for this somehow?

        // Make sizes not atomic
        let non_atomic: Vec<usize> = self.data.into_iter().map(
            |x| x.load(std::sync::atomic::Ordering::Acquire) as usize
        ).collect();
        let mut bv = BitVec::from_vec(non_atomic);
        bv.truncate(self.len); // Truncate the leftover bits in the last word
        bv
    }
}