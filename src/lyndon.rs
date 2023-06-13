use crate::kmer::{Base, Kmer};
use std::cmp::min;

pub trait Lyndon<const K: usize, T: Base> {
    fn rot_right(self) -> Self;
    fn necklace(self) -> Self;
    fn lmer(self) -> Self;
}

macro_rules! impl_Lyndon {
($($t:ty),+) => {$(
    impl<const K: usize> Lyndon<K, $t> for Kmer<K, $t> {
        #[inline]
        fn rot_right(self) -> Self {
            Self::from_int(((self.to_int() & 1) << (2 * K - 1)) | (self.to_int() >> 1))
        }
        fn necklace(self) -> Self {
            let mut res = self;
            let mut rot = self;
            for _ in 1..(2 * K) {
                rot = rot.rot_right();
                res = min(res, rot);
            }
            res
        }
        fn lmer(self) -> Self {
            self.canonical().necklace()
        }
    }
)*}}

impl_Lyndon!(u8, u16, u32, u64, u128);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_necklace() {
        let kmer = Kmer::<3, u32>::from_nucs(b"CGT");
        println!("{:06b}", kmer.to_int());
        println!("{:06b}", kmer.necklace().to_int());
        println!("{:06b}", kmer.rev_comp().to_int());
        println!("{:06b}", kmer.rev_comp().necklace().to_int());
    }

    #[test]
    fn test_necklace_rc() {
        for i in 0..1_000_000 {
            let kmer = Kmer::<15, u32>::from_int(i);
            assert_eq!(kmer.lmer(), kmer.rev_comp().lmer());
        }
    }
}
