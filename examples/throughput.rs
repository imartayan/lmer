use lmer::kmer::RawKmer;
use lmer::lyndon::Lyndon;
use lmer::rank::Ranker;
use lmer::utils::*;
use std::time::Instant;

const N: usize = 100_000;
type T = u64;

fn main() {
    bench_throughput::<11, { 2 * 11 - 1 }>();
    if T::BITS >= 64 {
        bench_throughput::<21, { 2 * 21 - 1 }>();
        bench_throughput::<31, { 2 * 31 - 1 }>();
    }
    if T::BITS >= 128 {
        bench_throughput::<41, { 2 * 41 - 1 }>();
        bench_throughput::<51, { 2 * 51 - 1 }>();
        bench_throughput::<61, { 2 * 61 - 1 }>();
    }
}

fn bench_throughput<const K: usize, const B: usize>() {
    let kmers = random_kmers::<K, T, RawKmer<K, T>>(N);
    let mut lmers = Vec::new();
    lmers.reserve(N);
    let mut ranks = Vec::new();
    ranks.reserve(N);
    println!("Throughput for k={}:", K);

    let now = Instant::now();
    for kmer in kmers {
        lmers.push(kmer.lmer());
    }
    let elapsed = now.elapsed().as_nanos();
    println!("{} ns/kmer to compute necklace", elapsed / N as u128);

    let now = Instant::now();
    let ranker = Ranker::<B, T>::new();
    for lmer in lmers {
        ranks.push(ranker.rank(lmer));
    }
    let elapsed = now.elapsed().as_nanos();
    println!("{} ns/lmer to compute rank", elapsed / N as u128);
}
