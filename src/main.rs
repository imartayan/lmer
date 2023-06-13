#![allow(dead_code)]
mod kmer;
mod lyndon;
mod reads;
mod utils;
use kmer::Kmer;
use lyndon::Lyndon;
use utils::*;

pub const K: usize = 5;
pub type T = u32;

fn main() {
    // all_lmers().iter().for_each(|&lmer| print_kmer_nucs(lmer));
    // let kmers = random_kmers(10);
    // kmers.iter().for_each(|&kmer| print_kmer_nucs(kmer));
    // println!();
    // kmers
    //     .iter()
    //     .for_each(|&kmer| print_kmer_nucs(kmer.necklace()));
    let seq = b"GCTCACGC";
    Kmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_nucs(kmer));
    println!();
    Kmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_nucs(kmer.necklace()));
    println!();
    Kmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_bin(kmer));
    println!();
    Kmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_bin(kmer.necklace()));
}
