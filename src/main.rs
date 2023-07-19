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

// const K: usize = 15;
// pub type T = u32;
// pub type L = u32;
// const B: usize = 5;
// type RB = roaring::RoaringBitmap;

const K: usize = 19;
pub type T = u64;
pub type L = u32;
const B: usize = 6;
type RB = roaring::RoaringTreemap;

// const K: usize = 31;
// pub type T = u64;
// pub type L = u64;
// const B: usize = 6;
// type RB = roaring::RoaringTreemap;

// const K: usize = 35;
// pub type T = u128;
// pub type L = u64;
// const B: usize = 7;

const N: usize = 2 * K - 1;
type TSet = std::collections::BTreeSet<T>;
type LSet = std::collections::BTreeSet<L>;

use std::fs::File;
use std::io::Write;

use sucds::mii_sequences::EliasFanoBuilder;
use sucds::Serializable;

fn main() {
    let mut set_kmers_cp = TSet::new();
    let mut set_lmers_idx = TSet::new();
    let mut set_ranks = LSet::new();

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
        set_ranks.insert(rank as L);

        // rb_kmers.insert(kmer.to_int());
        // rb_kmers_cp.insert(kmer_cp);
        // rb_lmers.insert(rank);
        // rb_lmers_idx.insert(rank_idx);
    }

    let mut out_kmers = File::create(format!("sorted_kmers_{}.txt", K)).unwrap();
    set_kmers_cp
        .iter()
        .for_each(|x| writeln!(out_kmers, "{}", x).unwrap());

    let mut out_lmers = File::create(format!("sorted_lmers_{}.txt", K)).unwrap();
    set_lmers_idx
        .iter()
        .for_each(|x| writeln!(out_lmers, "{}", x).unwrap());

    let mut out_ranks = File::create(format!("sorted_ranks_{}.txt", K)).unwrap();
    set_ranks
        .iter()
        .for_each(|x| writeln!(out_ranks, "{}", x).unwrap());

    // println!("Kmer bitmap: {} bytes", rb_kmers.serialized_size());
    // println!("Kmer cp bitmap: {} bytes", rb_kmers_cp.serialized_size());
    // println!("Lmer bitmap: {} bytes", rb_lmers.serialized_size());
    // println!("Lmer idx bitmap: {} bytes", rb_lmers_idx.serialized_size());
}
