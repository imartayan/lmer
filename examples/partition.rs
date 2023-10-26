use core::cmp;
use itertools::Itertools;
use lmer::kmer::RawKmer;
use lmer::lyndon::Lyndon;
use lmer::rank::Ranker;
use lmer::utils::random_kmers;
use std::collections::BTreeSet;

const N: usize = 1_000_000;
const RANK: bool = false;
const K: usize = 31;
const B: usize = 2 * K - 1;
type T = u64;

/// cost of bound + number of lmers
const FIXED_COST: f64 = (T::BITS * 2) as f64;
/// approximation factor
const EPS: f64 = 0.3;

/// cost to encode [`val[i]`, `val[j-1]`]
fn weight(val: &Vec<usize>, i: usize, j: usize) -> f64 {
    debug_assert!(i < j, "{i} >= {j}");
    debug_assert!(j <= val.len(), "{j} >= {}", val.len());
    let n = j - i;
    let u = val[j - 1] - val[i] + 1;
    if n == u {
        return FIXED_COST;
    }
    let n = n as f64;
    let u = u as f64;
    let l = (u / n).log2().ceil();
    let ef = 2.0 * n + n * l;
    let cost = FIXED_COST + if ef < u { ef } else { u };
    (cost / 8.0).ceil() * 8.0
}

/// compute a (1+Ɛ) approximation of the optimal partition
fn partition(val: &Vec<usize>) -> (Vec<usize>, f64) {
    let n = val.len();
    let plain_cost = weight(&val, 0, n);
    let q = (plain_cost / FIXED_COST).log(1.0 + EPS).floor() as usize;
    // let q = (1.0 / EPS1).log(1.0 + EPS).floor() as usize;

    let mut dist = vec![f64::INFINITY; n + 1];
    dist[0] = 0.0;
    let mut pred = vec![0usize; n + 1];
    let mut windows = vec![Some(1usize); q + 1];
    windows.push(None);
    let thresholds: Vec<f64> = (0..=q)
        .map(|h| FIXED_COST * (1.0 + EPS).powi(h as i32))
        .collect();

    for i in 0..n {
        for h in 0..=q {
            if let Some(w) = windows[h] {
                let mut j = cmp::max(w, i + 1);
                while j < n && weight(&val, i, j) <= thresholds[h] {
                    j += 1;
                }
                windows[h] = if weight(&val, i, j) > thresholds[h] {
                    Some(j)
                } else {
                    None
                };
            }
        }
        for &w in windows.iter().dedup() {
            let j = if let Some(j) = w { j - 1 } else { n };
            let d = dist[i] + weight(&val, i, j);
            if dist[j] > d {
                dist[j] = d;
                pred[j] = i;
            }
        }
    }

    let mut bounds = Vec::new();
    let mut j = n;
    while j > 0 {
        j = pred[j];
        bounds.push(j);
    }
    bounds.reverse();
    (bounds, dist[n])
}

fn main() {
    println!("K={K}, Ɛ={EPS}, {N} k-mers");
    let ranker = Ranker::<B, T>::new();
    let kmers = random_kmers::<K, T, RawKmer<K, T>>(2 * N);
    let mut lmers = BTreeSet::new();
    lmers.extend(kmers[..N].iter().map(|&kmer| {
        let lmer = kmer.lmer();
        if RANK {
            ranker.rank(lmer) as usize
        } else {
            lmer as usize
        }
    }));
    let val: Vec<_> = lmers.iter().map(|&x| x).collect();
    let n = val.len();

    println!("Plain cost/entry: {:.2}", weight(&val, 0, n) / n as f64);
    let (bounds, cost) = partition(&val);
    println!("Partition cost/entry: {:.2}", cost / n as f64);
    println!("using {} partition(s)", bounds.len());
    println!();

    lmers.extend(kmers[N..].iter().map(|&kmer| {
        let lmer = kmer.lmer();
        if RANK {
            ranker.rank(lmer) as usize
        } else {
            lmer as usize
        }
    }));
    let val2: Vec<_> = lmers.iter().map(|&x| x).collect();
    let n2 = val2.len();

    println!("After inserting {} new lmers:", n2 - n);
    println!("Plain cost/entry: {:.2}", weight(&val2, 0, n2) / n2 as f64);

    let bounds_val: Vec<_> = bounds.iter().map(|&i| val[i]).collect();
    let bounds2: Vec<_> = bounds_val
        .iter()
        .map(|v| val2.binary_search(v).unwrap())
        .collect();
    let cost2 = bounds2
        .iter()
        .zip(bounds2.iter().skip(1))
        .map(|(&i, &j)| weight(&val2, i, j))
        .sum::<f64>()
        + weight(&val2, *bounds2.last().unwrap(), n2);
    println!("Old partition cost/entry: {:.2}", cost2 / n2 as f64);
    println!("using {} partition(s)", bounds2.len());

    let (bounds2_new, cost2_new) = partition(&val2);
    println!("New partition cost/entry: {:.2}", cost2_new / n2 as f64);
    println!("using {} partition(s)", bounds2_new.len());
}
