#![allow(dead_code)]
use crate::kmer::{Base, Kmer, RawKmer};
use crate::lyndon::Lyndon;
use core::fmt::Binary;
use rand::Rng;
use std::collections::BTreeSet;

pub fn all_lmers<const K: usize>() -> BTreeSet<u32> {
    let mut lmers = BTreeSet::new();
    for i in 0..=RawKmer::<K, u32>::MASK {
        let kmer = RawKmer::<K, u32>::from_int(i);
        lmers.insert(kmer.lmer());
    }
    lmers
}

pub fn random_kmers<const K: usize, T: Base, KT: Kmer<K, T>>(n: usize) -> Vec<KT> {
    let mut res = Vec::new();
    let mut rng = rand::thread_rng();
    let mut kmer = KT::new();
    for i in 0..(n + K - 1) {
        let base: T = T::bases()[rng.gen_range(0..4)];
        if i < K - 1 {
            kmer = kmer.extend(base);
        } else {
            kmer = kmer.append(base);
            res.push(kmer);
        }
    }
    res
}

pub fn print_kmer_int<const K: usize, T: Base, KT: Kmer<K, T>>(kmer: KT) {
    println!("{}", kmer.to_int());
}

pub fn print_bin<T: Binary>(x: T) {
    println!("{:b}", x);
}

pub fn print_kmer_bin<const K: usize, T: Base, KT: Kmer<K, T>>(kmer: KT) {
    print_bin(kmer.to_int());
}

pub fn print_kmer_nucs<const K: usize, T: Base, KT: Kmer<K, T>>(kmer: KT) {
    kmer.to_nucs()
        .map(|x| x as char)
        .iter()
        .for_each(|c| print!("{}", c));
    println!();
}

pub fn print_lmers_ints<const K: usize>() {
    all_lmers::<K>()
        .iter()
        .for_each(|&lmer| println!("{}", lmer));
}
