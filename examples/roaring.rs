use std::env;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

// pub type T = u32;
// type RB = roaring::RoaringBitmap;

pub type T = u64;
type RB = roaring::RoaringTreemap;

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = args.get(1).expect("No filename given").as_str();

    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);

    let mut rb = RB::new();
    for line in reader.lines() {
        let val: T = line.unwrap().parse().unwrap();
        rb.insert(val);
    }

    println!("Roaring size in bits: {}", rb.serialized_size() * 8);
}
