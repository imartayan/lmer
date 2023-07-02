use crate::kmer::{Base, Kmer};
use std::cmp::min;

pub trait Lyndon<const K: usize, T: Base>: Kmer<K, T> {
    #[inline]
    fn rot_right(x: T) -> T {
        ((x & T::one()) << (2 * K - 2)) | (x >> 1)
    }
    fn lmer(self) -> T {
        let mut res = self.canonical().to_int() >> 1;
        let mut rot = res;
        for _ in 1..(2 * K - 1) {
            rot = Self::rot_right(rot);
            res = min(res, rot);
        }
        res
    }
}

impl<const K: usize, T: Base, KT: Kmer<K, T>> Lyndon<K, T> for KT {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::{RawKmer, RevComp};

    #[test]
    fn test_lmer_rc() {
        for i in 0..1_000_000 {
            let kmer = RawKmer::<15, u32>::from_int(i);
            assert_eq!(kmer.lmer(), kmer.rev_comp().lmer());
        }
    }
}
