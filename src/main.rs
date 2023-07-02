#![allow(dead_code)]
mod kmer;
mod lyndon;
mod rank;
mod reads;
mod utils;
// use kmer::{Base, Kmer, RawKmer};
use lyndon::Lyndon;
use rank::Ranker;
use utils::*;

// pub const K: usize = 11;
// pub type T = u32;
// const USELESS: T = 32 - 2 * K as T;
pub const K: usize = 15;
pub type T = u32;
const N: usize = 2 * K - 1;

fn main() {
    // all_lmers().iter().for_each(|&lmer| print_bin(lmer));
    // all_lmers()
    //     .iter()
    //     .for_each(|&lmer| println!("{:09b}", lmer));
    // all_lmers()
    //     .iter()
    //     .for_each(|&lmer| println!("{}", lmer.leading_zeros() - USELESS));

    let kmers = random_kmers::<K>(100000);
    // let lmers: Vec<T> = kmers.iter().map(|&kmer| kmer.lmer()).collect();
    let ranker = Ranker::<N, T>::new();
    let _ranks: Vec<T> = kmers.iter().map(|&kmer| ranker.rank(kmer.lmer())).collect();

    // let seq = b"GCTCACGC";
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_nucs(kmer));
    // println!();
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_nucs(kmer.necklace()));
    // println!();
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_bin(kmer));
    // println!();
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_bin(kmer.necklace()));
}
