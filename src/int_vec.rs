use std::sync::atomic::Ordering::{Acquire, Relaxed, Release};

#[derive(Debug, Clone)]
pub struct CompactIntVec {
    pub data: Vec<u64>,
    pub len: usize, // Number of stored integers. The last word may be only partially used.
    pub bit_width: usize,
}

impl CompactIntVec {
    pub fn new(len: usize, bit_width: usize) -> Self {
        assert!(bit_width <= 64);
        let n_words = (len * bit_width + 63) / 64;
        let data = vec![0; n_words];
        Self {data, len, bit_width}
    }

    pub fn get(&self, i: usize) -> usize {
        assert!(i < self.len);
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1;
            let bits = (self.data[word_idx] >> word_offset) & mask;
            bits as usize
        } else { // Combine bits from two words

            let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
            let n_bits2 = self.bit_width - n_bits1; // Rest of the bits from the start of the second word
            debug_assert!(n_bits1 + n_bits2 == self.bit_width);

            let x1 = self.data[word_idx] >> word_offset; // Tail of the first word
            let x2 = self.data[word_idx + 1] & ((1_u64 << n_bits2) - 1); // Head of the second word

            (x1 | (x2 << n_bits1)) as usize // Piece together
        }
    }

    pub fn set(&mut self, i: usize, x: usize) {
        assert!(i < self.len);
        debug_assert!((x as u64) < (1_u64 << self.bit_width));
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1; // Hopefully computed at compile time
            self.data[word_idx] &= !(mask << word_offset); // Clear the bits
            self.data[word_idx] |= (x as u64) << word_offset ; // Set new bits
        } else { // Combine bits from two words
            let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
            let n_bits2 = self.bit_width - n_bits1; // Rest of the bits from the start of the second word

            let mask1 = (1_u64 << n_bits1) - 1;
            let clearmask1 = !(mask1 << word_offset);
            let setmask1 = (x as u64 & mask1) << word_offset;

            self.data[word_idx] &= clearmask1; // Clear the bits
            self.data[word_idx] |= setmask1; // Set the bits

            let mask2 = (1_u64 << n_bits2) - 1;
            let clearmask2 = !mask2;
            let setmask2 = x as u64 >> n_bits1;

            self.data[word_idx + 1] &= clearmask2; // Clear the bits
            self.data[word_idx + 1] |= setmask2; // Set the bits
        }
    }

    pub fn from_vec(v: Vec::<usize>, bit_width: usize) -> Self {
        if let Some(v_max) = v.iter().max() {
            assert!(*v_max < (1_usize << bit_width));
        }
        let mut ret = Self::new(v.len(), bit_width);
        for (i, x) in v.into_iter().enumerate() {
            ret.set(i, x);
        }
        ret
    }

    pub fn from_atomic(other: AtomicCompactIntVec) -> Self {
        let len = other.len;
        let bit_width = other.bit_width;
        let data: Vec<u64> = other.data.into_iter().map(|x| x.load(Relaxed)).collect();
        Self{ data, len, bit_width }
    }

    pub fn serialize(&self, writer: &mut impl std::io::Write) {
        // Write length and bit width
        writer.write_all(&(self.len as u64).to_le_bytes()).unwrap();
        writer.write_all(&(self.bit_width as u64).to_le_bytes()).unwrap();

        // Write data
        bincode::serialize_into(writer, &self.data).unwrap();
    }

    pub fn load(reader: &mut impl std::io::Read) -> Self {
        let mut buf8 = [0_u8; 8];

        // Read length
        reader.read_exact(&mut buf8).unwrap();
        let len = u64::from_le_bytes(buf8) as usize;

        // Read bit width
        reader.read_exact(&mut buf8).unwrap();
        let bit_width = u64::from_le_bytes(buf8) as usize;

        // Read data
        let data = bincode::deserialize_from(reader).unwrap();

        Self{ data, len, bit_width }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn bit_width(&self) -> usize {
        self.bit_width
    }

    pub fn into_parts(self) -> (Vec<u64>, usize, usize) {
        (self.data, self.len, self.bit_width)
    }

    /// If new_len > current length, new elements are zero-initialized.
    /// If new_len < current length, the vector is truncated.
    pub fn resize(&mut self, new_len: usize) {
        let new_n_words = (new_len * self.bit_width).div_ceil(64);
        self.data.resize(new_n_words, 0);
        self.len = new_len;
    }

}


/// This is atomic only in a limited sense! See comments at member functions get and set.
pub struct AtomicCompactIntVec {
    pub data: Vec<std::sync::atomic::AtomicU64>,
    pub len: usize, // Number of stored integers. The last word may be only partially used.
    pub bit_width: usize,
}

impl AtomicCompactIntVec {
    pub fn new(len: usize, bit_width: usize) -> Self {
        assert!(bit_width <= 64);
        let n_words = (len * bit_width + 63) / 64;
        let data = (0..n_words).map(|_| std::sync::atomic::AtomicU64::new(0)).collect();
        Self {data, len, bit_width}
    }

    pub fn new_with_universe_size(len: usize, universe_size: usize) -> Self {
        let bit_width = universe_size.next_power_of_two().trailing_zeros() as usize;
        Self::new(len, bit_width)
    }

    /// This operation is not atomic in the sense that if some thread is writing to a
    /// value that is being read, the result could be mixed! It's atomic in the sense that
    /// it is safe to load a value even if another thread modifies a nearby value that is
    /// stored in the same word.
    pub fn get(&self, i: usize) -> usize {
        assert!(i < self.len);
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1;
            let bits = (self.data[word_idx].load(Acquire) >> word_offset) & mask;
            bits as usize
        } else { // Combine bits from two words

            let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
            let n_bits2 = self.bit_width - n_bits1; // Rest of the bits from the start of the second word
            debug_assert!(n_bits1 + n_bits2 == self.bit_width);

            let x1 = self.data[word_idx].load(Acquire) >> word_offset; // Tail of the first word
            let x2 = self.data[word_idx + 1].load(Acquire) & ((1_u64 << n_bits2) - 1); // Head of the second word

            (x1 | (x2 << n_bits1)) as usize // Piece together
        }
    }

    /// This operation is not atomic in the sense that if two threads try to modify
    /// the same value, then the result could be a mix of the two updates! This is
    /// atomic in the sense that modifying elements at nearby indices is ok even if
    /// they are in the same word.
    pub fn set(&self, i: usize, x: usize) {
        assert!(i < self.len);
        debug_assert!((x as u64) < (1_u64 << self.bit_width));
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1; // Hopefully computed at compile time
            self.data[word_idx].fetch_and(!(mask << word_offset), Release); // Clear the bits
            self.data[word_idx].fetch_or((x as u64) << word_offset, Release) ; // Set new bits
        } else { // Combine bits from two words
            let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
            let n_bits2 = self.bit_width - n_bits1; // Rest of the bits from the start of the second word

            let mask1 = (1_u64 << n_bits1) - 1;
            let clearmask1 = !(mask1 << word_offset);
            let setmask1 = (x as u64 & mask1) << word_offset;

            self.data[word_idx].fetch_and(clearmask1, Release); // Clear the bits
            self.data[word_idx].fetch_or(setmask1, Release); // Set the bits

            let mask2 = (1_u64 << n_bits2) - 1;
            let clearmask2 = !mask2;
            let setmask2 = x as u64 >> n_bits1;

            self.data[word_idx + 1].fetch_and(clearmask2, Release); // Clear the bits
            self.data[word_idx + 1].fetch_or(setmask2, Release); // Set the bits
        }
    }

    // Returns the underlying data, length and bit width
    pub fn into_parts(self) -> (Vec<usize>, usize, usize) {
        let data = self.data.into_iter().map(
            |x| x.load(Acquire) as usize
        ).collect();
        (data, self.len, self.bit_width) 
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn serialize_and_load() {
        use super::CompactIntVec;
        let mut vec = CompactIntVec::new(105, 30);
        for i in 0..100 {
            vec.set(i, i*i);
        }

        let mut buf: Vec<u8> = vec![];
        vec.serialize(&mut buf);

        let vec2 = CompactIntVec::load(&mut buf.as_slice());

        for i in 0..100 {
            assert_eq!(vec2.get(i), i*i);
        }
        for i in 100..vec2.len() {
            assert_eq!(vec2.get(i), 0);
        }

        assert_eq!(vec.data, vec2.data);
        assert_eq!(vec.len, vec2.len);
        assert_eq!(vec.bit_width, vec2.bit_width);
    }
}