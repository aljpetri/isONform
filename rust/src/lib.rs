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
pub mod parallel;
pub mod parasail;
pub mod poa;
pub mod pyset;
pub mod reads;
pub mod simd;
pub mod simplify;
pub mod sketch;
pub mod wfa;
/// Which configuration the port runs in. `ISONFORM_FAITHFUL` selects it.
///
/// Three points on one axis, because "how close to the reference" is the axis
/// every measurement in `PORTING.md` is taken along:
///
/// * [`Mode::Faithful`] (`ISONFORM_FAITHFUL=1`) --- **byte-identical to the
///   reference**, verified with `cmp` on `sirv_real` at 1 000/2 000/5 000/10 000/
///   20 000/50 000 reads and on `droso`. Every deliberate divergence off,
///   including CPython's set-iteration order ([`crate::pyset`]).
/// * [`Mode::Stable`] (unset --- **the default**) --- the faithful baseline plus
///   WFA2, which is the one optimisation measured never to cost accuracy at
///   depth: transcript parity or better at five of six depths, and 1.7--2.9x
///   faster.
/// * [`Mode::AllOpts`] (`ISONFORM_FAITHFUL=0`) --- every optimisation. Fastest,
///   and it loses transcripts; the arms interact, so it is worse than the sum of
///   its good parts.
///
/// Any single optimisation can still be switched on or off by its own variable,
/// which always wins over the mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Faithful,
    Stable,
    AllOpts,
}

pub fn mode() -> Mode {
    static M: std::sync::OnceLock<Mode> = std::sync::OnceLock::new();
    *M.get_or_init(
        || match std::env::var("ISONFORM_FAITHFUL").ok().as_deref() {
            Some("1") | Some("on") => Mode::Faithful,
            Some("0") | Some("off") => Mode::AllOpts,
            _ => Mode::Stable,
        },
    )
}

/// Reproduce the reference exactly, WFA2 included. Only [`Mode::Faithful`].
pub fn faithful() -> bool {
    mode() == Mode::Faithful
}

/// Keep the reference's semantics everywhere WFA2 is not involved --- node
/// identity, the POA schedule, the reference's own bugs, and CPython's set
/// order. True for both [`Mode::Faithful`] and [`Mode::Stable`].
pub fn reference_semantics() -> bool {
    mode() != Mode::AllOpts
}

pub mod weights;
pub mod wis;
