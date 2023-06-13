use crate::kmer::Kmer;
use crate::lyndon::Lyndon;
use crate::{K, T};
use rand::Rng;
use std::collections::BTreeSet;

pub fn all_lmers() -> BTreeSet<Kmer<K, T>> {
    let mut lmers = BTreeSet::new();
    for i in 0..=Kmer::<K, T>::MASK {
        let kmer = Kmer::<K, T>::from_int(i);
        lmers.insert(kmer.lmer());
    }
    lmers
}

pub fn random_kmers(n: usize) -> Vec<Kmer<K, T>> {
    let mut res = Vec::new();
    let mut rng = rand::thread_rng();
    let mut kmer = Kmer::<K, T>::new();
    for i in 0..n {
        let base: T = rng.gen_range(0..4);
        if i < K - 1 {
            kmer = kmer.extend(base);
        } else {
            kmer = kmer.append(base);
            res.push(kmer);
        }
    }
    res
}

pub fn print_kmer_int(kmer: Kmer<K, T>) {
    println!("{}", kmer.to_int());
}

pub fn print_kmer_bin(kmer: Kmer<K, T>) {
    println!("{:b}", kmer.to_int());
}

pub fn print_kmer_nucs(kmer: Kmer<K, T>) {
    kmer.to_nucs()
        .map(|x| x as char)
        .iter()
        .for_each(|c| print!("{}", c));
    println!();
}

pub fn print_lmers_ints() {
    all_lmers()
        .iter()
        .for_each(|&lmer| println!("{}", lmer.to_int()));
}

pub fn print_lmers_nucs() {
    all_lmers().iter().for_each(|&lmer| print_kmer_nucs(lmer));
}
