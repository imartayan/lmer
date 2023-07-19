use lmer::lyndon::{bin_rot_right, necklace};
use rand::Rng;
use std::collections::HashSet;

type T = u64;
const K: usize = 19;
const N: usize = 2 * K - 1;
const MASK: T = (1 << N) - 1;

use std::fs::File;
use std::io::Write;

fn main() {
    let mut rng = rand::thread_rng();
    let mut kmer = 0u64;
    for _ in 0..(N - 1) {
        kmer = (kmer << 1) | rng.gen_range(0..2);
    }
    let mut set_necks = HashSet::new();
    for _ in 0..1000 {
        kmer = ((kmer << 1) | rng.gen_range(0..2)) & MASK;
        set_necks.insert(necklace::<N, _>(kmer));
    }
    println!("{} necklaces", set_necks.len());

    let mut dot = File::create("sample.dot").unwrap();
    writeln!(dot, "digraph {{ node[shape=point]").unwrap();

    for &neck in set_necks.iter() {
        let mut rot = neck;
        for _ in 0..N {
            let s = (rot >> 1) | ((1 - (rot & 1)) << (N - 1));
            let ns = necklace::<N, _>(s);
            if set_necks.contains(&ns) {
                writeln!(dot, "{:b} -> {:b}", neck, ns).unwrap();
            }
            rot = bin_rot_right::<N, _>(rot);
        }
    }

    writeln!(dot, "}}").unwrap();
}
