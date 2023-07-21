use crate::kmer::{Base, Kmer};
use core::cmp::min;

#[inline]
pub fn bin_rot_right<const K: usize, T: Base>(x: T) -> T {
    ((x & T::one()) << (2 * K - 2)) | (x >> 1)
}

pub fn necklace<const K: usize, T: Base>(x: T) -> T {
    let mut res = x;
    let mut rot = x;
    for _ in 1..(2 * K - 1) {
        rot = bin_rot_right::<K, T>(rot);
        res = min(res, rot);
    }
    res
}

pub fn necklace_index<const K: usize, T: Base>(x: T) -> (T, usize) {
    let mut res = x;
    let mut rot = x;
    let mut idx = 0;
    for i in 1..(2 * K - 1) {
        rot = bin_rot_right::<K, T>(rot);
        if rot < res {
            res = rot;
            idx = i;
        }
    }
    (res, idx)
}

pub trait Lyndon<const K: usize, T: Base>: Kmer<K, T> {
    fn lmer(self) -> T {
        necklace::<K, T>(self.canonical().to_int() >> 1)
    }
    fn lmer_index(self) -> (T, usize) {
        necklace_index::<K, T>(self.canonical().to_int() >> 1)
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
