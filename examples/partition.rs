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

fn main() {
    println!("K={K}, Ɛ={EPS}, {N} k-mers");
    let kmers = random_kmers::<K, T, RawKmer<K, T>>(N);
    let ranker = Ranker::<B, T>::new();
    let lmers: BTreeSet<usize> = kmers
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
    let val: Vec<usize> = lmers.iter().map(|&x| x).collect();

    let n = val.len();
    let init_cost = weight(&val, 0, n);
    let q = (init_cost / FIXED_COST).log(1.0 + EPS).floor() as usize;
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

    println!("Initial cost/entry: {:.2}", init_cost / n as f64);
    println!("Partitioned cost/entry: {:.2}", dist[n] / n as f64);
    let mut bounds = vec![n];
    let mut j = n;
    while j > 0 {
        j = pred[j];
        bounds.push(j);
    }
    bounds.reverse();
    // println!("using {} partition(s)", bounds.len() - 1);
    // let stop = cmp::min(bounds.len() - 1, 10);
    // println!("{stop} first partition(s):");
    // for i in 0..stop {
    //     let (i, j) = (bounds[i], bounds[i + 1]);
    //     println!(
    //         "width={:.1e},\tcost={:.2}",
    //         val[j - 1] - val[i] + 1,
    //         weight(&val, i, j) / (j - i) as f64
    //     );
    // }
}
