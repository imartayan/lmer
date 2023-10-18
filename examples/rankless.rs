use lmer::kmer::RawKmer;
use lmer::lyndon::Lyndon;
use lmer::rank::Ranker;
use lmer::utils::*;
use num_format::{Locale, ToFormattedString};
use roaring::RoaringBitmap;
use std::collections::BTreeSet;
use std::time::Instant;
use sucds::mii_sequences::EliasFanoBuilder;
use sucds::Serializable;

// type T = u32;
type T = u64;

fn main() {
    bench_rankless::<1_000_000, 13, { 2 * 13 - 1 }>();
    bench_rankless::<1_000_000, 15, { 2 * 15 - 1 }>();
    if T::BITS >= 64 {
        bench_rankless_long::<1_000_000, 17, { 2 * 17 - 1 }>();
        bench_rankless_long::<1_000_000, 19, { 2 * 19 - 1 }>();
        bench_rankless_long::<1_000_000, 31, { 2 * 31 - 1 }>();
    }
}

fn bench_rankless<const N: usize, const K: usize, const B: usize>() {
    println!("K = {}, {} k-mers", K, N.to_formatted_string(&Locale::fr));
    println!();
    let kmers = random_kmers::<K, T, RawKmer<K, T>>(N);
    println!("Roaring without ranks");
    bench_roaring::<N, K, B, false>(&kmers);
    println!("Roaring with ranks");
    bench_roaring::<N, K, B, true>(&kmers);
    println!("EF without ranks");
    bench_ef::<N, K, B, false>(&kmers);
    println!("EF with ranks");
    bench_ef::<N, K, B, true>(&kmers);
}

fn bench_rankless_long<const N: usize, const K: usize, const B: usize>() {
    println!("K = {}, {} k-mers", K, N.to_formatted_string(&Locale::fr));
    println!();
    let kmers = random_kmers::<K, T, RawKmer<K, T>>(N);
    println!("EF without ranks");
    bench_ef::<N, K, B, false>(&kmers);
    println!("EF with ranks");
    bench_ef::<N, K, B, true>(&kmers);
}

fn gcd(a: usize, b: usize) -> usize {
    let mut a = a;
    let mut b = b;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn phi(n: usize) -> usize {
    (1..=n).filter(|&i| gcd(n, i) == 1).count()
}

fn num_necklaces(n: usize) -> usize {
    (1..=n)
        .filter(|d| n % d == 0)
        .map(|d| phi(n / d) << d)
        .sum::<usize>()
        / n
}

fn bench_roaring<const N: usize, const K: usize, const B: usize, const RANK: bool>(
    kmers: &Vec<RawKmer<K, T>>,
) {
    let u = if RANK { num_necklaces(B) } else { 1 << B };
    let mut roaring = RoaringBitmap::new();
    let ranker = Ranker::<B, T>::new();

    let values: BTreeSet<u32> = kmers
        .iter()
        .map(|kmer| {
            let lmer = kmer.lmer();
            if RANK {
                ranker.rank(lmer) as u32
            } else {
                lmer as u32
            }
        })
        .collect();
    let n_values = values.len();
    roaring.append(values).unwrap();
    println!("density = {:.2e}", roaring.len() as f64 / u as f64);
    println!(
        "{:.2} bits / necklace",
        roaring.serialized_size() as f64 * 8.0 / n_values as f64
    );

    let mut _count = 0u64;
    let now = Instant::now();
    for kmer in kmers {
        let lmer = kmer.lmer();
        if RANK {
            let rank = ranker.rank(lmer) as u32;
            if roaring.contains(rank) {
                _count += 1;
            }
        } else {
            let lmer = lmer as u32;
            if roaring.contains(lmer) {
                _count += 1;
            }
        }
    }
    let elapsed = now.elapsed().as_nanos();
    println!("{} ns / query", elapsed / N as u128);
    println!();
}

fn bench_ef<const N: usize, const K: usize, const B: usize, const RANK: bool>(
    kmers: &Vec<RawKmer<K, T>>,
) {
    let u = if RANK { num_necklaces(B) } else { 1 << B };
    let mut efb = EliasFanoBuilder::new(u, N).expect("Failed to create Elias-Fano Builder");
    let ranker = Ranker::<B, T>::new();

    let values: BTreeSet<usize> = kmers
        .iter()
        .map(|kmer| {
            let lmer = kmer.lmer();
            if RANK {
                ranker.rank(lmer) as usize
            } else {
                lmer as usize
            }
        })
        .collect();
    let n_values = values.len();
    efb.extend(values).unwrap();
    let ef = efb.build().enable_rank();
    println!("density = {:.2e}", ef.len() as f64 / u as f64);
    println!(
        "{:.2} bits / necklace",
        ef.size_in_bytes() as f64 * 8.0 / n_values as f64
    );

    let mut _count = 0u64;
    let now = Instant::now();
    for kmer in kmers {
        let lmer = kmer.lmer();
        if RANK {
            let rank = ranker.rank(lmer) as usize;
            if ef.successor(rank) == Some(rank) {
                _count += 1;
            }
        } else {
            let lmer = lmer as usize;
            if ef.successor(lmer) == Some(lmer) {
                _count += 1;
            }
        }
    }
    let elapsed = now.elapsed().as_nanos();
    println!("{} ns / query", elapsed / N as u128);
    println!();
}
