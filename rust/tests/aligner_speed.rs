//! Wall-clock and agreement for the three alignment backends, on recorded calls.
//!
//! The question is not which aligner is fastest in the abstract --- it is which
//! one is fastest *on the alignments isONform actually performs*, and what it
//! costs in agreement with the reference. Both halves matter: WFA2 is faster
//! than the scalar port and disagrees; the linked C library is faster still and
//! cannot disagree, because it is what the reference calls.
//!
//! ```text
//! PARASAIL_CASES=<recorded calls> cargo test --release --test aligner_speed -- --nocapture
//! ```
//!
//! Record the cases with `bench/dump_reference.py --record-parasail`. Use a
//! corpus of *corrected* reads: the sequences being aligned are consensuses, and
//! their length is what decides the ratio. `bench/corpus/sirv_small` yields
//! ~200 bp; a real corpus yields ~600 bp with a 1 kb tail, and the backends do
//! not rank the same way on both.

use isonform::parasail::Scoring;
use std::time::Instant;

struct Case {
    s1: Vec<u8>,
    s2: Vec<u8>,
    cigar: String,
    score: i32,
    sc: Scoring,
}

fn load(path: &str) -> Vec<Case> {
    let text = std::fs::read_to_string(path).expect("read PARASAIL_CASES");
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let (s1, s2) = (f[0].as_bytes().to_vec(), f[1].as_bytes().to_vec());
        if s1.is_empty() || s2.is_empty() {
            continue;
        }
        out.push(Case {
            s1,
            s2,
            cigar: f[2].to_string(),
            score: f[3].parse().expect("score"),
            sc: Scoring {
                match_score: f[4].parse().expect("match"),
                mismatch: f[5].parse().expect("mismatch"),
                open: f[6].parse().expect("open"),
                ext: f[7].parse().expect("ext"),
            },
        });
    }
    out
}

/// Best of `repeats`, after one warm-up pass. Returns (millis, cigar hits,
/// score hits) against the *recorded reference* answer.
fn time_backend(
    cases: &[Case],
    repeats: usize,
    f: &dyn Fn(&Case) -> (i32, String),
) -> (f64, usize, usize) {
    for c in cases {
        std::hint::black_box(f(c));
    }
    let (mut cigar_hits, mut score_hits) = (0usize, 0usize);
    let mut best = f64::MAX;
    for r in 0..repeats {
        let t = Instant::now();
        let (mut ch, mut sh) = (0usize, 0usize);
        for c in cases {
            let (score, cigar) = f(c);
            if cigar == c.cigar {
                ch += 1;
            }
            if score == c.score {
                sh += 1;
            }
            std::hint::black_box(score);
        }
        best = best.min(t.elapsed().as_secs_f64() * 1000.0);
        if r == 0 {
            cigar_hits = ch;
            score_hits = sh;
        }
    }
    (best, cigar_hits, score_hits)
}

#[test]
fn the_three_backends_measured_on_recorded_calls() {
    let Ok(path) = std::env::var("PARASAIL_CASES") else {
        eprintln!("SKIPPED: set PARASAIL_CASES to the recorded calls");
        return;
    };
    let repeats: usize = std::env::var("ALIGNER_BENCH_REPEATS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let cases = load(&path);
    assert!(!cases.is_empty(), "no usable cases in {path}");

    let mut lens: Vec<usize> = cases.iter().map(|c| c.s1.len().max(c.s2.len())).collect();
    lens.sort_unstable();
    let at = |q: f64| lens[((lens.len() - 1) as f64 * q) as usize];
    println!(
        "{} cases  maxlen median={} p90={} max={}",
        cases.len(),
        at(0.5),
        at(0.9),
        at(1.0)
    );

    let scalar = time_backend(&cases, repeats, &|c| {
        let a = isonform::parasail::semiglobal(&c.s1, &c.s2, c.sc);
        (a.score, a.cigar)
    });
    let wfa = time_backend(&cases, repeats, &|c| {
        match isonform::wfa::semiglobal(&c.s1, &c.s2, c.sc) {
            Some(a) => (a.score, a.cigar),
            // Exactly what the call site does when WFA2 declines.
            None => {
                let a = isonform::parasail::semiglobal(&c.s1, &c.s2, c.sc);
                (a.score, a.cigar)
            }
        }
    });

    println!(
        "\n{:<38}{:>10}{:>10}{:>14}{:>14}",
        "backend", "wall(ms)", "speedup", "cigar==ref", "score==ref"
    );
    let row = |name: &str, (ms, ch, sh): (f64, usize, usize)| {
        println!(
            "{:<38}{:>10.1}{:>9.2}x{:>9} {:>3.1}%{:>9} {:>3.1}%",
            name,
            ms,
            scalar.0 / ms,
            ch,
            100.0 * ch as f64 / cases.len() as f64,
            sh,
            100.0 * sh as f64 / cases.len() as f64
        )
    };
    row("parasail (scalar port)", scalar);
    row("WFA2 (+ parasail fallback)", wfa);
    #[cfg(feature = "parasail-ffi")]
    {
        let ffi = time_backend(&cases, repeats, &|c| {
            let a = isonform::parasail_ffi::semiglobal(&c.s1, &c.s2, c.sc);
            (a.score, a.cigar)
        });
        row("parasail C (libparasail-sys)", ffi);
        // Not a benchmark assertion --- an identity one. If this backend ever
        // stops reproducing the recorded reference exactly, it is not calling
        // what we think it is, and the speed is beside the point.
        assert_eq!(
            ffi.1,
            cases.len(),
            "the linked library must reproduce every recorded CIGAR"
        );
        assert_eq!(ffi.2, cases.len(), "...and every recorded score");
    }
}
