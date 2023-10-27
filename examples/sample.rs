use lmer::constants::{CANON_BITS, K, KT, LT, ROT_BITS};
use lmer::kmer::{Kmer, RawKmer};
use lmer::lyndon::Lyndon;
use lmer::rank::Ranker;
use lmer::utils::*;
use std::fs::File;
use std::io::Write;

type TSet = std::collections::BTreeSet<KT>;
type LSet = std::collections::BTreeSet<LT>;

fn main() {
    let mut set_kmers_cp = TSet::new();
    let mut set_lmers_idx = TSet::new();
    let mut set_ranks = LSet::new();

    let kmers = random_kmers::<K, KT, RawKmer<K, KT>>(1_000_000);
    let ranker = Ranker::<CANON_BITS, KT>::new();
    for kmer in kmers {
        let kmer_cp = kmer.canonical().to_int() >> 1;
        let (lmer, idx) = kmer.lmer_index();
        let rank = ranker.rank(lmer);
        let rank_idx = (rank << ROT_BITS) | idx as KT;

        set_kmers_cp.insert(kmer_cp);
        set_lmers_idx.insert(rank_idx);
        set_ranks.insert(rank as LT);
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
}
