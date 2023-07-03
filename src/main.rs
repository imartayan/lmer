#![allow(dead_code)]
#![allow(unused_imports)]
mod kmer;
mod lyndon;
mod rank;
mod reads;
mod utils;

use kmer::Kmer;
use lyndon::Lyndon;
use rank::Ranker;
use utils::*;

pub type T = u32;
const K: usize = 15;
const B: usize = 5;
type RB = roaring::RoaringBitmap;

// pub type T = u64;
// const K: usize = 31;
// const B: usize = 6;
// type RB = roaring::RoaringTreemap;

const N: usize = 2 * K - 1;
type TSet = std::collections::BTreeSet<T>;

use std::fs::File;
use std::io::Write;

use sucds::mii_sequences::EliasFanoBuilder;
use sucds::Serializable;

fn main() {
    let mut set_kmers_cp = TSet::new();
    let mut set_lmers_idx = TSet::new();

    // let mut rb_kmers = RB::new();
    // let mut rb_kmers_cp = RB::new();
    // let mut rb_lmers = RB::new();
    // let mut rb_lmers_idx = RB::new();

    let kmers = random_kmers::<K>(1_000_000);
    let ranker = Ranker::<N, T>::new();
    for kmer in kmers {
        let kmer_cp = kmer.canonical().to_int() >> 1;
        let (lmer, idx) = kmer.lmer_index();
        let rank = ranker.rank(lmer);
        let rank_idx = (rank << B) | idx as T;

        set_kmers_cp.insert(kmer_cp);
        set_lmers_idx.insert(rank_idx);

        // rb_kmers.insert(kmer.to_int());
        // rb_kmers_cp.insert(kmer_cp);
        // rb_lmers.insert(rank);
        // rb_lmers_idx.insert(rank_idx);
    }

    // let mut out_kmers = File::create("sorted_kmers.txt").unwrap();
    // set_kmers_cp
    //     .iter()
    //     .for_each(|x| writeln!(out_kmers, "{}", x).unwrap());

    // let mut out_lmers = File::create("sorted_lmers.txt").unwrap();
    // set_lmers_idx
    //     .iter()
    //     .for_each(|x| writeln!(out_lmers, "{}", x).unwrap());

    // set_lmers_idx.iter().for_each(|x| println!("{}", x))

    // println!("Kmer bitmap: {} bytes", rb_kmers.serialized_size());
    // println!("Kmer cp bitmap: {} bytes", rb_kmers_cp.serialized_size());
    // println!("Lmer bitmap: {} bytes", rb_lmers.serialized_size());
    // println!("Lmer idx bitmap: {} bytes", rb_lmers_idx.serialized_size());

    // let mut efb = EliasFanoBuilder::new(1 << (N + 1), set_kmers_cp.len())
    //     .expect("Failed to create Elias-Fano Builder");
    // efb.extend(set_kmers_cp.iter().map(|&x| x as usize))
    //     .unwrap();
    // let ef = efb.build();
    // println!("Kmer cp EF: {} bytes", ef.size_in_bytes());

    // let mut efb = EliasFanoBuilder::new(1 << (N + 1), set_lmers_idx.len())
    //     .expect("Failed to create Elias-Fano Builder");
    // efb.extend(set_lmers_idx.iter().map(|&x| x as usize))
    //     .unwrap();
    // let ef = efb.build();
    // println!("Lmer idx EF: {} bytes", ef.size_in_bytes());

    // let seq = b"GCTCACGC";
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_nucs(kmer));
    // println!();
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_nucs(kmer.necklace()));
    // println!();
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_bin(kmer));
    // println!();
    // RawKmer::<K, T>::iter_from_nucs(seq.iter()).for_each(|kmer| print_kmer_bin(kmer.necklace()));
}
