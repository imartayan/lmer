use crate::kmer::{Base, Kmer};
use core::cmp::min;

#[inline]
pub fn bin_rot_right<const N: usize, T: Base>(x: T) -> T {
    ((x & T::one()) << (N - 1)) | (x >> 1)
}

pub fn necklace<const N: usize, T: Base>(x: T) -> T {
    let mut res = x;
    let mut rot = x;
    for _ in 1..N {
        rot = bin_rot_right::<N, T>(rot);
        res = min(res, rot);
    }
    res
}

pub fn necklace_index<const N: usize, T: Base>(x: T) -> (T, usize) {
    let mut res = x;
    let mut rot = x;
    let mut idx = 0;
    for i in 1..N {
        rot = bin_rot_right::<N, T>(rot);
        if rot < res {
            res = rot;
            idx = i;
        }
    }
    (res, idx)
}

pub trait Lyndon<const K: usize, T: Base>: Kmer<K, T> {
    #[inline]
    fn bin_rot_right(x: T) -> T {
        ((x & T::one()) << (2 * K - 2)) | (x >> 1)
    }
    fn lmer(self) -> T {
        let mut res = self.canonical().to_int() >> 1;
        let mut rot = res;
        for _ in 1..(2 * K - 1) {
            rot = Self::bin_rot_right(rot);
            res = min(res, rot);
        }
        res
    }
    fn lmer_index(self) -> (T, usize) {
        let mut res = self.canonical().to_int() >> 1;
        let mut idx = 0;
        let mut rot = res;
        for i in 1..(2 * K - 1) {
            rot = Self::bin_rot_right(rot);
            if rot < res {
                res = rot;
                idx = i;
            }
        }
        (res, idx)
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
