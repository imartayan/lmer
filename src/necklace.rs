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
        self.min_pos
            .iter()
            .map(|&pos| pos.wrapping_sub(self.pos) % W)
    }

    pub fn insert(&mut self, u: T) {
        let mut i = self.deq.len();
        while i > 0 {
            if self.deq[i - 1].0 <= u {
                break;
            }
            i -= 1;
        }
        self.deq.truncate(i);
        self.min_pos.truncate(i);
        self.deq.push_back((u, self.pos));
        if i == 0 {
            self.min_pos.push_back(self.deq[0].1);
        } else {
            if self.pos == self.deq[0].1 {
                self.deq.pop_front();
                self.min_pos.pop_front();
            }
            while self.min_pos.len() < self.deq.len()
                && self.deq[self.min_pos.len()].0 == self.deq[0].0
            {
                self.min_pos.push_back(self.deq[self.min_pos.len()].1);
            }
        }
        self.pos = self.pos.wrapping_add(1) % W;
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

    pub fn new(word: $T) -> Self {
        Self {
            word: word & Self::MASK,
            min_queue: LexMinQueue::new(),
        }
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
            ((N - Self::M)..(N - 1))
                .map(|p| (self.rotation(p), p))
                .min()
                .unwrap(),
        )
    }

    pub fn insert(&mut self, x: $T) {
        self.word = ((self.word << 1) & Self::MASK) | x;
        self.min_queue.insert(self.word & Self::MIN_MASK);
    }
}
)*}}

impl_t!(u8, u16, u32, u64, u128);

#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 8;
    const M: usize = 4;
    const W: usize = N - M + 1;

    #[test]
    fn test_lex_min_queue() {
        let mut queue = LexMinQueue::<W, _>::new();
        for x in [1, 4, 1, 2, 3] {
            queue.insert(x);
            println!("{:?}", queue);
        }
    }

    #[test]
    fn test_necklace() {
        let mut queue = NecklaceQueue::<N, W, u64>::new(u64::MAX);
        for x in [0u64, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1] {
            queue.insert(x);
            println!("{:?}", queue.get_necklace_pos());
        }
    }
}
