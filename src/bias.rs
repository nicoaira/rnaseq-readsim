//! Library biases: which candidate fragments make it into the library, and
//! how many different molecules the library holds.
//!
//! Every bias is a probability of keeping a candidate fragment. A candidate
//! is a transcript, a fragment length, a position and a strand, drawn as
//! without biases; it is kept with the product of the probabilities below,
//! otherwise a new candidate is drawn from the start (transcript included).
//! A bias therefore changes how many fragments each transcript gives, as it
//! does in a real library, and not only where the fragments lie.
//!
//! - **Fragment GC content** (`--gc-bias-center`, `--gc-bias-width`):
//!   `exp(-((gc - center) / width)² / 2)`, `gc` being the fraction of G and C
//!   in the fragment. PCR and fragmentation favour fragments of intermediate
//!   GC content; how strongly differs between protocols and kits.
//! - **Degradation** (`--degradation`): poly(A) selection of partly degraded
//!   RNA keeps only the 3′ piece of each molecule. With pieces of
//!   exponentially distributed length (mean `m`), a fragment that starts
//!   `d` bases before the 3′ end of its transcript is kept with probability
//!   `exp(-d / m)`. Coverage then rises towards the 3′ end, and long
//!   transcripts give fewer fragments.
//! - **Random-hexamer priming** (`--hexamer-bias`): the bases at the start of
//!   each read are those where a random primer annealed; primers do not
//!   anneal equally well everywhere. A position weight matrix over the first
//!   [`HEXAMER_POSITIONS`] bases of each read gives the probability of
//!   keeping the fragment; `strength` raises every weight to that power (0 =
//!   no bias, 1 = the built-in matrix). Read 2 uses half the strength.
//! - **Library complexity** (`--library-molecules`): the library holds `M`
//!   different molecules and every read pair is drawn from them at random,
//!   with replacement, so that reads of the same molecule (duplicates) arise
//!   as in a real library amplified by PCR. Each molecule is generated from
//!   its own random stream, so a duplicate is the same fragment; its
//!   sequencing errors and qualities are drawn anew.

use crate::rng::Rng;

/// Number of read positions covered by the random-hexamer weights.
pub const HEXAMER_POSITIONS: usize = 9;

/// Built-in random-hexamer weights for the first bases of a read, columns
/// A, C, G, T. Made up for the simulator: the weights give position-dependent
/// base compositions over the first nine bases of the kind random priming
/// produces (FastQC's "Per base sequence content" fails), not the profile of
/// any particular kit.
const HEXAMER_WEIGHTS: [[f64; 4]; HEXAMER_POSITIONS] = [
    [1.45, 0.75, 1.25, 0.55],
    [0.80, 1.30, 0.70, 1.20],
    [1.30, 0.70, 1.20, 0.80],
    [0.85, 1.20, 0.90, 1.05],
    [1.15, 0.90, 1.10, 0.85],
    [0.90, 1.12, 0.92, 1.06],
    [1.08, 0.94, 1.06, 0.92],
    [0.95, 1.05, 0.96, 1.04],
    [1.03, 0.98, 1.02, 0.97],
];

/// Index of a base in the columns A, C, G, T; `None` for other letters.
#[inline]
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn complement_index(i: usize) -> usize {
    3 - i
}

/// Position weights of one read, each scaled so that its largest value is 1.
#[derive(Clone, Debug)]
struct ReadWeights([[f64; 4]; HEXAMER_POSITIONS]);

impl ReadWeights {
    fn new(strength: f64) -> Self {
        let mut w = [[0.0; 4]; HEXAMER_POSITIONS];
        for (row, base) in w.iter_mut().zip(HEXAMER_WEIGHTS.iter()) {
            let max = base.iter().map(|x| x.powf(strength)).fold(0.0, f64::max);
            for (o, x) in row.iter_mut().zip(base.iter()) {
                *o = x.powf(strength) / max;
            }
        }
        ReadWeights(w)
    }
}

/// The biases of one library. All of them are optional.
#[derive(Clone, Debug, Default)]
pub struct Biases {
    gc: Option<(f64, f64)>,
    degradation: Option<f64>,
    hexamer: Option<[ReadWeights; 2]>,
}

/// Settings of the biases, as given on the command line.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BiasSettings {
    /// GC fraction most favoured, between 0 and 1.
    pub gc_center: Option<f64>,
    /// Width of the GC preference (standard deviation, in GC fraction).
    pub gc_width: f64,
    /// Mean length of the 3′ pieces kept from degraded molecules, in bases.
    pub degradation: Option<f64>,
    /// Strength of the random-hexamer bias (0 = none, 1 = built-in matrix).
    pub hexamer: f64,
}

impl BiasSettings {
    pub fn any(&self) -> bool {
        self.gc_center.is_some() || self.degradation.is_some() || self.hexamer > 0.0
    }
}

impl Biases {
    pub fn new(s: &BiasSettings) -> Result<Self, String> {
        let gc = match s.gc_center {
            None => None,
            Some(c) => {
                if !(0.0..=1.0).contains(&c) {
                    return Err(format!(
                        "the GC bias centre must be between 0 and 1 (got {c})"
                    ));
                }
                if !s.gc_width.is_finite() || s.gc_width <= 0.0 {
                    return Err(format!(
                        "the GC bias width must be a positive number (got {})",
                        s.gc_width
                    ));
                }
                Some((c, s.gc_width))
            }
        };
        let degradation = match s.degradation {
            None => None,
            Some(m) if m.is_finite() && m > 0.0 => Some(m),
            Some(m) => {
                return Err(format!(
                    "the mean length of the kept 3' pieces must be a positive number of bases (got {m})"
                ));
            }
        };
        if !s.hexamer.is_finite() || s.hexamer < 0.0 || s.hexamer > 3.0 {
            return Err(format!(
                "the random-hexamer bias strength must be between 0 and 3 (got {})",
                s.hexamer
            ));
        }
        let hexamer = (s.hexamer > 0.0).then(|| {
            [
                ReadWeights::new(s.hexamer),
                ReadWeights::new(s.hexamer / 2.0),
            ]
        });
        Ok(Biases {
            gc,
            degradation,
            hexamer,
        })
    }

    pub fn is_active(&self) -> bool {
        self.gc.is_some() || self.degradation.is_some() || self.hexamer.is_some()
    }

    /// Probability of keeping a fragment that starts at `start` in a
    /// transcript of length `tlen`, given its sequence and the strand of
    /// read 1.
    pub fn keep_probability(
        &self,
        fragment: &[u8],
        start: usize,
        tlen: usize,
        read1_on_transcript_strand: bool,
    ) -> f64 {
        let mut p = 1.0;
        if let Some(m) = self.degradation {
            // distance from the start of the fragment to the 3' end
            let d = (tlen - start) as f64;
            p *= (-d / m).exp();
        }
        if let Some([w1, w2]) = &self.hexamer {
            // read 1 starts at one end of the fragment, read 2 at the other
            let (r1, r2) = if read1_on_transcript_strand {
                (Start::Forward, Start::Reverse)
            } else {
                (Start::Reverse, Start::Forward)
            };
            p *= read_start_weight(fragment, r1, w1) * read_start_weight(fragment, r2, w2);
        }
        if p == 0.0 {
            return 0.0;
        }
        if let Some((center, width)) = self.gc {
            let gc = gc_fraction(fragment);
            let z = (gc - center) / width;
            p *= (-0.5 * z * z).exp();
        }
        p
    }
}

#[derive(Clone, Copy)]
enum Start {
    /// The read is on the fragment's strand: it starts at its first base.
    Forward,
    /// The read is the reverse complement: it starts at its last base.
    Reverse,
}

#[inline]
fn read_start_weight(fragment: &[u8], start: Start, w: &ReadWeights) -> f64 {
    let n = fragment.len().min(HEXAMER_POSITIONS);
    let mut p = 1.0;
    for (pos, row) in w.0.iter().enumerate().take(n) {
        let i = match start {
            Start::Forward => base_index(fragment[pos]),
            Start::Reverse => base_index(fragment[fragment.len() - 1 - pos]).map(complement_index),
        };
        if let Some(i) = i {
            p *= row[i];
        }
    }
    p
}

/// Fraction of G and C among the A, C, G and T of a sequence (0.5 if it has
/// none of them).
pub fn gc_fraction(seq: &[u8]) -> f64 {
    let (mut gc, mut acgt) = (0u32, 0u32);
    for &b in seq {
        match base_index(b) {
            Some(1) | Some(2) => {
                gc += 1;
                acgt += 1
            }
            Some(_) => acgt += 1,
            None => {}
        }
    }
    if acgt == 0 {
        0.5
    } else {
        gc as f64 / acgt as f64
    }
}

/// Keeps a candidate with probability `p`.
#[inline]
pub fn accept(rng: &mut Rng, p: f64) -> bool {
    p >= 1.0 || rng.next_f64() < p
}

/// Expected fraction of read pairs that are the first read of their molecule
/// when `n` reads are drawn with replacement from `m` molecules.
pub fn expected_unique_fraction(n: u64, m: u64) -> f64 {
    if n == 0 || m == 0 {
        return 1.0;
    }
    let (n, m) = (n as f64, m as f64);
    m * (1.0 - (-n / m).exp()) / n
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings() -> BiasSettings {
        BiasSettings {
            gc_center: None,
            gc_width: 0.15,
            degradation: None,
            hexamer: 0.0,
        }
    }

    #[test]
    fn no_bias_keeps_everything() {
        let b = Biases::new(&settings()).unwrap();
        assert!(!b.is_active());
        assert_eq!(b.keep_probability(b"ACGTACGTAC", 0, 100, false), 1.0);
    }

    #[test]
    fn gc_bias_prefers_the_centre() {
        let mut s = settings();
        s.gc_center = Some(0.5);
        let b = Biases::new(&s).unwrap();
        let mid = b.keep_probability(b"ACGTACGTACGTACGTACGT", 0, 20, true);
        let at = b.keep_probability(b"ATATATATATATATATATAT", 0, 20, true);
        let gc = b.keep_probability(b"GCGCGCGCGCGCGCGCGCGC", 0, 20, true);
        assert!((mid - 1.0).abs() < 1e-12);
        assert!(at < 0.01 && gc < 0.01);
    }

    #[test]
    fn degradation_favours_the_three_prime_end() {
        let mut s = settings();
        s.degradation = Some(100.0);
        let b = Biases::new(&s).unwrap();
        let near = b.keep_probability(b"ACGT", 996, 1000, true);
        let far = b.keep_probability(b"ACGT", 0, 1000, true);
        assert!((near - (-0.04f64).exp()).abs() < 1e-12);
        assert!((far - (-10.0f64).exp()).abs() < 1e-12);
    }

    #[test]
    fn hexamer_weights_depend_on_the_read_start_and_strand() {
        let mut s = settings();
        s.hexamer = 1.0;
        let b = Biases::new(&s).unwrap();
        // read 1 on the transcript strand reads the fragment from its first
        // base; read 1 on the other strand from the complement of its last
        let f = b"AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCC";
        let fwd = b.keep_probability(f, 0, 40, true);
        let rev = b.keep_probability(f, 0, 40, false);
        assert!(fwd > 0.0 && rev > 0.0);
        assert!((fwd - rev).abs() > 1e-6);
        // the largest weight of every position is 1
        let w = ReadWeights::new(1.0);
        for row in &w.0 {
            assert!((row.iter().cloned().fold(0.0, f64::max) - 1.0).abs() < 1e-12);
        }
        // strength 0 is no bias
        let w0 = ReadWeights::new(0.0);
        assert!(
            w0.0.iter()
                .all(|row| row.iter().all(|x| (*x - 1.0).abs() < 1e-12))
        );
    }

    #[test]
    fn rejects_bad_settings() {
        let mut s = settings();
        s.gc_center = Some(1.5);
        assert!(Biases::new(&s).is_err());
        let mut s = settings();
        s.gc_center = Some(0.5);
        s.gc_width = 0.0;
        assert!(Biases::new(&s).is_err());
        let mut s = settings();
        s.degradation = Some(-1.0);
        assert!(Biases::new(&s).is_err());
        let mut s = settings();
        s.hexamer = f64::NAN;
        assert!(Biases::new(&s).is_err());
    }

    #[test]
    fn gc_fraction_counts_only_bases() {
        assert_eq!(gc_fraction(b"GGCC"), 1.0);
        assert_eq!(gc_fraction(b"ATNN"), 0.0);
        assert_eq!(gc_fraction(b"NNNN"), 0.5);
        assert_eq!(gc_fraction(b"ACGT"), 0.5);
    }

    #[test]
    fn unique_fraction_of_a_library() {
        assert!((expected_unique_fraction(1000, 1_000_000_000) - 1.0).abs() < 1e-5);
        let f = expected_unique_fraction(1_000_000, 1_000_000);
        assert!((f - (1.0 - (-1.0f64).exp())).abs() < 1e-12);
    }
}
