//! isONform — reference-free isoform reconstruction from long-read data.
//!
//! A Rust port of <https://github.com/aljpetri/isONform>. The Python reference
//! is normative for observable behaviour; see `PORTING.md` at the repository
//! root for the method, the reconnaissance, and the reference defects found so
//! far.
//!
//! Only [`cli`] exists yet. It is deliberately first: the plan's method point 1
//! is CLI parity locked by unit tests, before any algorithm.

pub mod cli;
