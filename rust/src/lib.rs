//! isONform — reference-free isoform reconstruction from long-read data.
//!
//! A Rust port of <https://github.com/aljpetri/isONform>. The Python reference
//! is normative for observable behaviour; see `PORTING.md` at the repository
//! root for the method, the reconnaissance, and the reference defects found so
//! far.
//!
//! [`cli`] came first: the plan's method point 1 is CLI parity locked by unit
//! tests, before any algorithm. [`graph`] is the first algorithm stage —
//! the interval graph and its construction, ported from
//! `modules/GraphGeneration.py`.

pub mod align;
pub mod anchors;
pub mod batch_merge;
pub mod cli;
pub mod driver;
pub mod fastq;
pub mod graph;
pub mod graph_build;
pub mod intervals;
pub mod isoforms;
pub mod minimizers;
pub mod parasail;
pub mod poa;
pub mod reads;
pub mod simplify;
pub mod wis;
