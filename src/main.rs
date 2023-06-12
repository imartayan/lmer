#![allow(dead_code)]
mod kmer;
mod lyndon;
mod reads;
use kmer::Kmer;
use lyndon::Lyndon;
use std::collections::BTreeSet;

const K: usize = 5;
type T = u32;

fn main() {
    let mut lmers = BTreeSet::new();
    for i in 0..=Kmer::<K, T>::MASK {
        let kmer = Kmer::<K, T>::from_int(i);
        lmers.insert(kmer.lmer());
    }
    lmers
        .iter()
        .map(|lmer| lmer.to_nucs().map(|x| x as char))
        .for_each(|nucs| {
            nucs.iter().for_each(|c| print!("{}", c));
            println!();
        });
}
