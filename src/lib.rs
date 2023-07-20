pub mod kmer;
pub mod lyndon;
pub mod rank;
pub mod reads;
pub mod utils;

// Loads runtime-provided constants for which declarations
// will be generated at `$OUT_DIR/constants.rs`.
pub mod constants {
    include!(concat!(env!("OUT_DIR"), "/constants.rs"));
}
