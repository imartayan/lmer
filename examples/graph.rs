use lmer::constants::{CANON_BITS, K, KT};
use lmer::lyndon::{bin_rot_right, necklace};
use rand::Rng;
use std::collections::HashSet;

const MASK: KT = (1 << CANON_BITS) - 1;

use std::fs::File;
use std::io::Write;

fn main() {
    let mut rng = rand::thread_rng();
    let mut kmer = 0u64;
    for _ in 0..(CANON_BITS - 1) {
        kmer = (kmer << 1) | rng.gen_range(0..2);
    }
    let mut set_necks = HashSet::new();
    for _ in 0..1000 {
        kmer = ((kmer << 1) | rng.gen_range(0..2)) & MASK;
        set_necks.insert(necklace::<K, _>(kmer));
    }
    println!("{} necklaces", set_necks.len());

    let mut dot = File::create("graph_necklaces_bin.dot").unwrap();
    writeln!(dot, "digraph {{ node[shape=point]").unwrap();

    for &neck in set_necks.iter() {
        let mut rot = neck;
        for _ in 0..CANON_BITS {
            let s = (rot >> 1) | ((1 - (rot & 1)) << (CANON_BITS - 1));
            let ns = necklace::<K, _>(s);
            if set_necks.contains(&ns) {
                writeln!(dot, "{:b} -> {:b}", neck, ns).unwrap();
            }
            rot = bin_rot_right::<K, _>(rot);
        }
    }

    writeln!(dot, "}}").unwrap();
}
