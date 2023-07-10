// use lmer::kmer::Kmer;
use lmer::lyndon::Lyndon;
use lmer::rank::Ranker;
use lmer::utils::random_kmers;

type T = u64;
const K: usize = 19;
const N: usize = 2 * K - 1;
type R = u32;
type TSet = std::collections::BTreeSet<R>;

use std::fs::File;
use std::io::Write;

// use sucds::mii_sequences::EliasFanoBuilder;
// use sucds::Serializable;

fn main() {
    let mut set_ranks = TSet::new();

    let kmers = random_kmers::<K>(1_000_000);
    let ranker = Ranker::<N, T>::new();
    for kmer in kmers {
        let lmer = kmer.lmer();
        let rank = ranker.rank(lmer) as R;
        set_ranks.insert(rank);
    }

    let mut out_ranks = File::create(format!("sorted_ranks_{}.txt", K)).unwrap();
    set_ranks
        .iter()
        .for_each(|x| writeln!(out_ranks, "{}", x).unwrap());

    // let rank_max = ranker.rank(1 << N - 1) as usize;

    // let mut efb = EliasFanoBuilder::new(rank_max + 1, set_ranks.len())
    //     .expect("Failed to create Elias-Fano Builder");
    // efb.extend(set_ranks.iter().map(|&x| x as usize)).unwrap();
    // let ef = efb.build();
    // println!("Ranks EF: {} bytes", ef.size_in_bytes());
}
