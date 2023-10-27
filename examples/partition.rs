use lmer::constants::{CANON_BITS, K, KT};
use lmer::kmer::RawKmer;
use lmer::lyndon::Lyndon;
use lmer::partition::Partitioner;
use lmer::rank::Ranker;
use lmer::utils::random_kmers;
use std::collections::BTreeSet;

/// approximation factor
const EPS: f64 = 0.3;
const N: usize = 1_000_000;
const RANK: bool = false;

fn main() {
    println!("K={K}, Ɛ={EPS}, {N} k-mers");
    let partitioner = Partitioner::new();
    let ranker = Ranker::<CANON_BITS, KT>::new();
    let kmers = random_kmers::<K, KT, RawKmer<K, KT>>(2 * N);
    let mut lmers = BTreeSet::new();

    lmers.extend(kmers[..N].iter().map(|&kmer| {
        let lmer = kmer.lmer();
        if RANK {
            ranker.rank(lmer)
        } else {
            lmer
        }
    }));
    let val: Vec<_> = lmers.iter().copied().collect();
    let n = val.len();

    println!(
        "Plain cost/entry: {} bits",
        partitioner.cost(&val, 0, n) / n
    );
    let (partition, cost) = partitioner.partition(&val, EPS);
    println!("Partition cost/entry: {:.2} bits", cost as f64 / n as f64);
    println!("using {} partition(s)", partition.len());
    println!();

    lmers.extend(kmers[N..].iter().map(|&kmer| {
        let lmer = kmer.lmer();
        if RANK {
            ranker.rank(lmer)
        } else {
            lmer
        }
    }));
    let val2: Vec<_> = lmers.iter().copied().collect();
    let n2 = val2.len();

    println!("After inserting {} new lmers:", n2 - n);
    println!(
        "Plain cost/entry: {} bits",
        partitioner.cost(&val2, 0, n2) / n2
    );
    let cost_alt = partitioner.cost_with_partition(&val2, &partition);
    println!(
        "Old partition cost/entry: {:.2} bits",
        cost_alt as f64 / n2 as f64
    );
    println!("using {} partition(s)", partition.len());
    let (partition2, cost2) = partitioner.partition(&val2, EPS);
    println!(
        "New partition cost/entry: {:.2} bits",
        cost2 as f64 / n2 as f64
    );
    println!("using {} partition(s)", partition2.len());
}
