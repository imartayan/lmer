use core::cmp;
use itertools::Itertools;
use num_traits::cast::AsPrimitive;
use num_traits::int::PrimInt;
use std::marker::PhantomData;
use std::mem::size_of;

pub struct Partitioner<T: PrimInt + AsPrimitive<usize>> {
    fixed_cost: usize,
    _phantom: PhantomData<T>,
}

impl<T: PrimInt + AsPrimitive<usize>> Partitioner<T> {
    pub fn new() -> Self {
        Self {
            fixed_cost: size_of::<T>() * 8 * 2,
            _phantom: PhantomData,
        }
    }

    /// cost to encode [`val[i]`, `val[j-1]`]
    pub fn cost(&self, val: &Vec<T>, i: usize, j: usize) -> usize {
        debug_assert!(i < j, "{i} >= {j}");
        debug_assert!(j <= val.len(), "{j} >= {}", val.len());
        let n = j - i;
        let u = (val[j - 1] - val[i]).as_() + 1;
        if n == u {
            return self.fixed_cost;
        }
        let l = u.div_ceil(n).next_power_of_two().ilog2() as usize; // lg(u/n)
        let ef = 2 * n + n * l;
        let cost = self.fixed_cost + if ef < u { ef } else { u };
        cost.next_multiple_of(8)
    }

    /// compute a (1+Ɛ) approximation of the optimal partition
    pub fn partition(&self, val: &Vec<T>, eps: f64) -> (Vec<T>, usize) {
        debug_assert!(!val.is_empty(), "empty vector");
        let n = val.len();
        let plain_cost = self.cost(val, 0, n);
        let q = (plain_cost as f64 / self.fixed_cost as f64).log(1.0 + eps) as usize;
        // let q = (1.0 / EPS1).log(1.0 + eps) as usize;

        let mut dist = vec![usize::MAX; n + 1]; // distance from 0
        dist[0] = 0;
        let mut pred = vec![0usize; n + 1];
        let mut windows = vec![Some(1usize); q + 1];
        windows.push(None);
        let mut t = self.fixed_cost as f64;
        let mut thresholds = vec![self.fixed_cost];
        for _ in 1..=q {
            t *= 1.0 + eps;
            thresholds.push(t as usize);
        }

        for i in 0..n {
            for h in 0..=q {
                if let Some(w) = windows[h] {
                    let mut j = cmp::max(w, i + 1);
                    while j < n && self.cost(val, i, j) <= thresholds[h] {
                        j += 1;
                    }
                    windows[h] = if self.cost(val, i, j) > thresholds[h] {
                        Some(j)
                    } else {
                        None
                    };
                }
            }
            for &w in windows.iter().dedup() {
                let j = if let Some(j) = w { j - 1 } else { n };
                let d = dist[i].saturating_add(self.cost(val, i, j));
                if dist[j] > d {
                    dist[j] = d;
                    pred[j] = i;
                }
            }
        }

        let mut partition = Vec::new();
        let mut i = n;
        while i > 0 {
            i = pred[i];
            partition.push(val[i]);
        }
        partition.reverse();
        (partition, dist[n])
    }

    pub fn cost_with_partition(&self, val: &Vec<T>, partition: &Vec<T>) -> usize {
        debug_assert!(!val.is_empty(), "empty vector");
        let n = val.len();
        let bounds: Vec<_> = partition
            .iter()
            .map(|v| val.binary_search(v).unwrap())
            .collect();
        bounds
            .iter()
            .zip(bounds.iter().skip(1))
            .map(|(&i, &j)| self.cost(val, i, j))
            .sum::<usize>()
            + self.cost(val, *bounds.last().unwrap(), n)
    }
}
