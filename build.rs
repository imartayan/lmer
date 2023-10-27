/// This is a hack to support dynamic K values.
/// K values are implemented as a const generic in our code
/// as we expect it to remain constant across executions
/// and benefit from compile-time optimizations.
/// This build script will set the value of K at compile-time
/// from an environment variable, so one can easily build
/// the project "just in time" with the desired K value.
/// This will not re-build if the K value does not change.
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=K");

    let out_dir: std::path::PathBuf = std::env::var("OUT_DIR")
        .expect("Failed to obtain OUT_DIR")
        .into();
    let mut code = Vec::new();

    let k: usize = std::env::var("K")
        .unwrap_or_else(|_| "31".into())
        .parse()
        .expect("Failed to parse K");
    assert!(k >= 1, "K must be ≥ 1");
    assert!(k < 64, "K must be < 64");
    assert!(k % 2 == 1, "K must be odd");
    code.push(format!("pub const K: usize = {k};"));

    let kmer_bits = 2 * k;
    code.push(format!("pub const KMER_BITS: usize = {kmer_bits};"));

    let canon_bits = 2 * k - 1;
    code.push(format!("pub const CANON_BITS: usize = {canon_bits};"));

    let t = select_type(kmer_bits);
    code.push(format!("pub type T = {t};"));

    let rot_bits = canon_bits.next_power_of_two().ilog2() as usize;
    code.push(format!("pub const ROT_BITS: usize = {rot_bits};"));

    let lmer_bits = kmer_bits - rot_bits;
    code.push(format!("pub const LMER_BITS: usize = {lmer_bits};"));

    let lmer_t = select_type(lmer_bits);
    code.push(format!("pub type LmerT = {lmer_t};"));

    std::fs::write(out_dir.join("constants.rs"), code.join("\n"))
        .expect("Failed to write const file");
}

fn select_type(n_bits: usize) -> &'static str {
    match n_bits.next_power_of_two() {
        8 => "u8",
        16 => "u16",
        32 => "u32",
        64 => "u64",
        _ => "u128",
    }
}
