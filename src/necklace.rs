use core::cmp::min;
use num_traits::int::PrimInt;
use std::collections::VecDeque;

#[derive(Debug)]
pub struct LexMinQueue<const W: usize, T: Ord + Copy> {
    deq: VecDeque<(T, usize)>,
    min_pos: VecDeque<usize>,
    pos: usize,
}

impl<const W: usize, T: Ord + Copy> LexMinQueue<W, T> {
    pub fn new() -> Self {
        Self {
            deq: VecDeque::with_capacity(W),
            min_pos: VecDeque::with_capacity(W),
            pos: 0,
        }
    }

    pub fn iter_min_pos(&self) -> impl Iterator<Item = usize> + '_ {
        self.min_pos.iter().map(|&pos| (pos + W - self.pos) % W)
    }

    pub fn insert_full<I: DoubleEndedIterator<Item = T>>(&mut self, vals: I) {
        self.deq.clear();
        self.min_pos.clear();
        let mut vals = vals.rev();
        let mut min = vals.next().unwrap();
        let mut pos = (self.pos + W - 1) % W;
        self.deq.push_front((min, pos));
        for u in vals.take(W - 1) {
            pos = (pos + W - 1) % W;
            if u <= min {
                min = u;
                self.deq.push_front((min, pos));
            }
        }
        while self.min_pos.len() < self.deq.len() && self.deq[self.min_pos.len()].0 == self.deq[0].0
        {
            self.min_pos.push_back(self.deq[self.min_pos.len()].1);
        }
    }

    pub fn insert(&mut self, u: T) {
        if !self.deq.is_empty() && self.deq[0].1 == self.pos {
            self.deq.pop_front();
            self.min_pos.pop_front();
        }
        let mut i = self.deq.len();
        while i > 0 && self.deq[i - 1].0 > u {
            i -= 1;
        }
        self.deq.truncate(i);
        self.min_pos.truncate(i);
        self.deq.push_back((u, self.pos));
        while self.min_pos.len() < self.deq.len() && self.deq[self.min_pos.len()].0 == self.deq[0].0
        {
            self.min_pos.push_back(self.deq[self.min_pos.len()].1);
        }
        self.pos = (self.pos + 1) % W;
    }

    pub fn insert2(&mut self, u: T, v: T) {
        let next_pos = (self.pos + 1) % W;
        if !self.deq.is_empty() && self.deq[0].1 == self.pos {
            self.deq.pop_front();
            self.min_pos.pop_front();
        }
        if !self.deq.is_empty() && self.deq[0].1 == next_pos {
            self.deq.pop_front();
            self.min_pos.pop_front();
        }
        let w = min(u, v);
        let mut i = self.deq.len();
        while i > 0 && self.deq[i - 1].0 > w {
            i -= 1;
        }
        self.deq.truncate(i);
        self.min_pos.truncate(i);
        if u <= v {
            self.deq.push_back((u, self.pos));
        }
        self.deq.push_back((v, next_pos));
        while self.min_pos.len() < self.deq.len() && self.deq[self.min_pos.len()].0 == self.deq[0].0
        {
            self.min_pos.push_back(self.deq[self.min_pos.len()].1);
        }
        self.pos = (next_pos + 1) % W;
    }
}

impl<const W: usize, T: Ord + Copy> Default for LexMinQueue<W, T> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug)]
pub struct NecklaceQueue<const N: usize, const W: usize, T: PrimInt> {
    word: T,
    min_queue: LexMinQueue<W, T>,
}

macro_rules! impl_t {
($($T:ty),+) => {$(
impl<const N: usize, const W: usize> NecklaceQueue<N, W, $T> {
    const M: usize = N - W + 1;
    const MASK: $T = (1 << N) - 1;
    const MIN_MASK: $T = (1 << Self::M) - 1;

    pub fn new() -> Self {
        Self {
            word: 0,
            min_queue: LexMinQueue::new(),
        }
    }

    pub fn new_from_word(word: $T) -> Self {
        let mut res = Self::new();
        res.insert_full(word);
        res
    }

    #[inline]
    fn rotation(&self, p: usize) -> $T {
        ((self.word << p) & Self::MASK) | (self.word >> (N - p))
    }

    pub fn get_necklace_pos(&self) -> ($T, usize) {
        min(
            self.min_queue
                .iter_min_pos()
                .map(|p| (self.rotation(p), p))
                .min()
                .unwrap(),
            (W..N)
                .map(|p| (self.rotation(p), p))
                .min()
                .unwrap(),
        )
    }

    pub fn insert_full(&mut self, word: $T) {
        self.word = word & Self::MASK;
        let vals = (0..W).map(|p|
            (word >> (N - p - Self::M)) & Self::MIN_MASK
        );
        self.min_queue.insert_full(vals);
    }

    pub fn insert(&mut self, x: $T) {
        self.word = ((self.word << 1) & Self::MASK) | (x & 0b1);
        self.min_queue.insert(self.word & Self::MIN_MASK);
    }

    pub fn insert2(&mut self, x: $T) {
        self.word = ((self.word << 2) & Self::MASK) | (x & 0b11);
        self.min_queue.insert2((self.word >> 1) & Self::MIN_MASK, self.word & Self::MIN_MASK);
    }
}

impl<const N: usize, const W: usize> Default for NecklaceQueue<N, W, $T> {
    fn default() -> Self {
    Self::new()
    }
}
)*}}

impl_t!(u8, u16, u32, u64, u128);

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    const N: usize = 8;
    const M: usize = 5;
    const W: usize = N - M + 1;

    #[test]
    fn test_lex_min_queue_insert_full() {
        let mut min_queue = LexMinQueue::<W, _>::new();
        min_queue.insert_full([2, 1, 2, 1].iter());
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 3, W - 1],
            "{:?}",
            min_queue
        );
    }

    #[test]
    fn test_lex_min_queue_insert() {
        let mut min_queue = LexMinQueue::<W, _>::new();
        min_queue.insert(3);
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 1],
            "{:?}",
            min_queue
        );
        min_queue.insert(1);
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 1],
            "{:?}",
            min_queue
        );
        min_queue.insert(2);
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 2],
            "{:?}",
            min_queue
        );
        min_queue.insert(3);
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 3],
            "{:?}",
            min_queue
        );
        min_queue.insert(1);
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 4, W - 1],
            "{:?}",
            min_queue
        );
        min_queue.insert(2);
        assert_eq!(
            min_queue.iter_min_pos().collect_vec(),
            vec![W - 2],
            "{:?}",
            min_queue
        );
    }

    #[test]
    fn test_necklace_queue() {
        let mut necklace_queue = NecklaceQueue::<N, W, u64>::new_from_word(0b10010110);
        assert_eq!(necklace_queue.get_necklace_pos(), (0b00101101, 1));
        necklace_queue.insert(0);
        assert_eq!(necklace_queue.get_necklace_pos(), (0b00001011, N - 2));
    }
}
