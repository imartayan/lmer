# Lmer

Sampling 1 million consecutive Kmers and Lmers, and writing them to a file:
```sh
cargo r -r --example sample
```
By default `K=31` but you can specify it as follows:
```sh
K=15 cargo r -r --example sample
```

Compressing sorted integers with Roaring Bitmaps:
```sh
cargo r -r --example roaring -- sorted_values.txt
```
