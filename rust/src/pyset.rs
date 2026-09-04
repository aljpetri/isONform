//! CPython's `set` iteration order, for integer keys.
//!
//! # Why this exists
//!
//! `compute_equal_reads` takes `list(current_node_support)[0]` as an isoform's id
//! and `list(current_node_support)` as its member order, and that order decides
//! which sequence seeds spoa --- which is order-sensitive. Finding 28 diverged
//! here, ordering members ascending instead, on the grounds that reproducing
//! CPython's order would mean modelling its probing and resize rule. This is that
//! model, so the port can be byte-identical to the reference rather than merely
//! equivalent.
//!
//! # What has to be reproduced
//!
//! A CPython set is an open-addressing table of `mask + 1` slots, and iteration
//! is **slot order**, so the answer is a function of the values and of the order
//! they were inserted (which decides who wins a collision). For `int` keys
//! `hash(x) == x`, except `hash(-1) == -2`; read ids are positive, so the hash is
//! the value.
//!
//! The probe sequence (CPython 3.10 `Objects/setobject.c`): start at
//! `hash & mask`, scan up to `LINEAR_PROBES` **consecutive** slots when they fit
//! below the mask, then jump with `i = i * 5 + 1 + perturb` and
//! `perturb >>= PERTURB_SHIFT`. The table grows when `fill * 5 >= mask * 3`, to
//! `used * 4` (or `used * 2` above 50 000) rounded up to a power of two, and a
//! resize **reinserts in old slot order** --- so growth history is part of the
//! answer, not just the final size.
//!
//! Only insertion and iteration are modelled. That is all `compute_equal_reads`
//! needs on the sets whose order reaches the output: every one of them is built
//! fresh by `set.intersection(edge_supp)`, and `edge_supp` is a **list**, so
//! CPython iterates that list and inserts matches in list order. Deletion
//! (`pop`, `-=`) only touches `node_support_left`, whose order does not reach the
//! result.

const MINSIZE: usize = 8;
const LINEAR_PROBES: usize = 9;
const PERTURB_SHIFT: u32 = 5;

#[derive(Clone, Copy, PartialEq)]
enum Slot {
    Empty,
    Used(u64),
}

/// A CPython set of non-negative integer keys, insert and iterate only.
pub struct PySet {
    table: Vec<Slot>,
    mask: usize,
    fill: usize,
    used: usize,
}

impl Default for PySet {
    fn default() -> Self {
        Self::new()
    }
}

impl PySet {
    pub fn new() -> Self {
        Self {
            table: vec![Slot::Empty; MINSIZE],
            mask: MINSIZE - 1,
            fill: 0,
            used: 0,
        }
    }

    /// `set.add`. No dummy handling: this models insert-only sets.
    pub fn add(&mut self, key: u64) {
        if self.insert(key) {
            // `fill * 5 >= mask * 3` --- the load factor CPython grows at.
            if self.fill * 5 >= self.mask * 3 {
                let minused = if self.used > 50_000 {
                    self.used * 2
                } else {
                    self.used * 4
                };
                self.resize(minused);
            }
        }
    }

    /// Places `key` if absent. Returns whether a slot was newly occupied.
    fn insert(&mut self, key: u64) -> bool {
        let mut perturb = key;
        let mut i = (key as usize) & self.mask;
        loop {
            // The linear run is taken only when it fits below the mask, exactly
            // as `probes = (i + LINEAR_PROBES <= mask) ? LINEAR_PROBES : 0`.
            let probes = if i + LINEAR_PROBES <= self.mask {
                LINEAR_PROBES
            } else {
                0
            };
            for j in 0..=probes {
                let at = i + j;
                match self.table[at] {
                    Slot::Empty => {
                        self.table[at] = Slot::Used(key);
                        self.fill += 1;
                        self.used += 1;
                        return true;
                    }
                    Slot::Used(k) if k == key => return false,
                    Slot::Used(_) => {}
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i
                .wrapping_mul(5)
                .wrapping_add(1)
                .wrapping_add(perturb as usize))
                & self.mask;
        }
    }

    /// `set_table_resize`: grow to a power of two above `minused` and reinsert
    /// every live entry **in old slot order**.
    fn resize(&mut self, minused: usize) {
        let mut newsize = MINSIZE;
        while newsize <= minused {
            newsize <<= 1;
        }
        let old = std::mem::replace(&mut self.table, vec![Slot::Empty; newsize]);
        self.mask = newsize - 1;
        self.fill = 0;
        self.used = 0;
        for slot in old {
            if let Slot::Used(k) = slot {
                self.insert(k);
            }
        }
    }

    /// `list(s)` --- slot order.
    pub fn iter_order(&self) -> Vec<u64> {
        self.table
            .iter()
            .filter_map(|s| match s {
                Slot::Used(k) => Some(*k),
                Slot::Empty => None,
            })
            .collect()
    }
}

/// `list(set(values))` for a known insertion order.
pub fn order_of<I: IntoIterator<Item = u64>>(values: I) -> Vec<u64> {
    let mut s = PySet::new();
    for v in values {
        s.add(v);
    }
    s.iter_order()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Replays vectors generated by the real interpreter. Without
    /// `PYSET_VECTORS` this is a no-op, like the other differential oracles.
    #[test]
    fn matches_cpython_on_recorded_vectors() {
        let Ok(path) = std::env::var("PYSET_VECTORS") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("PYSET_VECTORS unreadable");
        let (mut n, mut bad) = (0usize, 0usize);
        for line in data.lines().filter(|l| !l.is_empty()) {
            let (ins, want) = line.split_once('\t').expect("insertion\twant");
            let ins: Vec<u64> = ins.split(',').map(|x| x.parse().unwrap()).collect();
            let want: Vec<u64> = want.split(',').map(|x| x.parse().unwrap()).collect();
            n += 1;
            let got = order_of(ins.iter().copied());
            if got != want {
                if bad < 3 {
                    eprintln!(
                        "MISMATCH n={} \n  want {:?}\n  got  {:?}",
                        ins.len(),
                        &want[..want.len().min(16)],
                        &got[..got.len().min(16)]
                    );
                }
                bad += 1;
            }
        }
        eprintln!("pyset oracle: {n} vectors, {bad} mismatches");
        assert_eq!(bad, 0, "{bad} of {n} vectors disagree with CPython");
    }
}
