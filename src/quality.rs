//! Base qualities and substitution errors.
//!
//! A quality profile gives, for read 1 and read 2 and for every cycle, the
//! distribution of the Phred quality value, optionally conditional on the
//! quality of the previous cycle (a first-order Markov chain, so that a read
//! that has started to go bad tends to stay bad). Profiles are either the
//! built-in one or learned from a FASTQ file and stored as JSON.

use std::fmt::Write as _;
use std::io::BufRead;
use std::path::Path;

use serde::Deserialize;

use crate::alias::vose;
use crate::error::{Context, Error, Result};
use crate::fasta::open_text;
use crate::rng::{Rng, threshold_u64};

/// Largest Phred value that can be written with the Phred+33 encoding.
pub const MAX_QUALITY: u8 = 93;
/// Quality given to `N` calls (`#`).
pub const N_QUALITY: u8 = 2;
/// Value of the `format` field of a profile file.
pub const PROFILE_FORMAT: &str = "rnaseq-readsim-quality-profile";
/// Version of the profile file layout.
pub const PROFILE_VERSION: u32 = 1;
/// Longest read that `profile` accepts.
pub const MAX_PROFILE_CYCLES: usize = 10_000;
/// Number of cycles the built-in profile is defined for; longer reads repeat
/// the last cycle.
pub const BUILTIN_CYCLES: usize = 300;

/// The probability that a base with Phred quality `q` is wrong:
/// `10^(-q/10)`.
pub fn error_probability(q: u8) -> f64 {
    10f64.powf(-(q as f64) / 10.0)
}

/// The quality distributions of one read (read 1 or read 2).
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ReadProfile {
    /// The Phred values that occur, in increasing order ("states").
    pub qualities: Vec<u8>,
    /// `cycles[c][s]`: weight (count or frequency) of state `s` at cycle `c`.
    pub cycles: Vec<Vec<f64>>,
    /// `transitions[c - 1][p][s]`: weight of state `s` at cycle `c` given
    /// state `p` at cycle `c - 1`. Absent for independent cycles.
    #[serde(default)]
    pub transitions: Option<Vec<Vec<Vec<f64>>>>,
}

/// A quality profile: one entry for read 1 and, optionally, one for read 2.
#[derive(Clone, Debug, Deserialize, PartialEq)]
pub struct Profile {
    pub format: String,
    pub version: u32,
    #[serde(default)]
    pub description: String,
    pub reads: Vec<ReadProfile>,
}

impl ReadProfile {
    fn validate(&self, what: &str) -> std::result::Result<(), String> {
        let s = self.qualities.len();
        if s == 0 {
            return Err(format!("{what}: 'qualities' is empty"));
        }
        for w in self.qualities.windows(2) {
            if w[0] >= w[1] {
                return Err(format!("{what}: 'qualities' must be strictly increasing"));
            }
        }
        if *self.qualities.last().expect("non-empty") > MAX_QUALITY {
            return Err(format!(
                "{what}: quality values above {MAX_QUALITY} cannot be written as Phred+33"
            ));
        }
        if self.cycles.is_empty() {
            return Err(format!("{what}: 'cycles' is empty"));
        }
        let check_row = |row: &[f64], place: &str, need_mass: bool| {
            if row.len() != s {
                return Err(format!(
                    "{what}: {place} has {} values, expected {s} (one per quality)",
                    row.len()
                ));
            }
            if row.iter().any(|w| !w.is_finite() || *w < 0.0) {
                return Err(format!(
                    "{what}: {place} contains a value that is not a finite, non-negative number"
                ));
            }
            if need_mass && row.iter().sum::<f64>() <= 0.0 {
                return Err(format!("{what}: {place} has no positive weight"));
            }
            Ok(())
        };
        for (c, row) in self.cycles.iter().enumerate() {
            check_row(row, &format!("cycle {}", c + 1), true)?;
        }
        if let Some(t) = &self.transitions {
            if t.len() != self.cycles.len() - 1 {
                return Err(format!(
                    "{what}: 'transitions' has {} entries, expected {} (cycles - 1)",
                    t.len(),
                    self.cycles.len() - 1
                ));
            }
            for (c, matrix) in t.iter().enumerate() {
                if matrix.len() != s {
                    return Err(format!(
                        "{what}: transition into cycle {} has {} rows, expected {s}",
                        c + 2,
                        matrix.len()
                    ));
                }
                for row in matrix {
                    check_row(row, &format!("transition into cycle {}", c + 2), false)?;
                }
            }
        }
        Ok(())
    }
}

impl Profile {
    pub fn validate(&self) -> std::result::Result<(), String> {
        if self.format != PROFILE_FORMAT {
            return Err(format!(
                "not a quality profile ('format' is '{}', expected '{PROFILE_FORMAT}')",
                self.format
            ));
        }
        if self.version != PROFILE_VERSION {
            return Err(format!(
                "profile version {} is not supported (this program reads version {PROFILE_VERSION})",
                self.version
            ));
        }
        if self.reads.is_empty() || self.reads.len() > 2 {
            return Err(format!(
                "'reads' must have one or two entries, found {}",
                self.reads.len()
            ));
        }
        for (i, r) in self.reads.iter().enumerate() {
            r.validate(&format!("read {}", i + 1))?;
        }
        Ok(())
    }

    /// Reads and checks a profile file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .context(|| format!("cannot read quality profile {}", path.display()))?;
        let profile: Profile = serde_json::from_str(&text)
            .context(|| format!("{} is not a valid quality profile", path.display()))?;
        profile
            .validate()
            .context(|| format!("{}", path.display()))?;
        Ok(profile)
    }

    /// The profile of read 1 (`mate` 0) or read 2 (`mate` 1). A profile with
    /// a single entry is used for both reads.
    pub fn read(&self, mate: usize) -> &ReadProfile {
        self.reads.get(mate).unwrap_or(&self.reads[0])
    }

    /// Serialises the profile as JSON, one cycle per line.
    pub fn to_json(&self) -> String {
        fn number(out: &mut String, v: f64) {
            if v.fract() == 0.0 && v.abs() < 1e15 {
                let _ = write!(out, "{}", v as i64);
            } else {
                let _ = write!(out, "{v}");
            }
        }
        fn row(out: &mut String, values: &[f64]) {
            out.push('[');
            for (i, &v) in values.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                number(out, v);
            }
            out.push(']');
        }
        let mut out = String::new();
        out.push_str("{\n");
        let _ = writeln!(out, "\"format\": \"{PROFILE_FORMAT}\",");
        let _ = writeln!(out, "\"version\": {},", self.version);
        let _ = writeln!(
            out,
            "\"description\": {},",
            serde_json::to_string(&self.description).expect("a string serialises")
        );
        out.push_str("\"reads\": [\n");
        for (r, read) in self.reads.iter().enumerate() {
            out.push_str("{\n\"qualities\": [");
            for (i, q) in read.qualities.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                let _ = write!(out, "{q}");
            }
            out.push_str("],\n\"cycles\": [\n");
            for (c, cycle) in read.cycles.iter().enumerate() {
                row(&mut out, cycle);
                out.push_str(if c + 1 < read.cycles.len() {
                    ",\n"
                } else {
                    "\n"
                });
            }
            out.push(']');
            if let Some(t) = &read.transitions {
                out.push_str(",\n\"transitions\": [\n");
                for (c, matrix) in t.iter().enumerate() {
                    out.push('[');
                    for (p, values) in matrix.iter().enumerate() {
                        if p > 0 {
                            out.push(',');
                        }
                        row(&mut out, values);
                    }
                    out.push(']');
                    out.push_str(if c + 1 < t.len() { ",\n" } else { "\n" });
                }
                out.push(']');
            }
            out.push_str(if r + 1 < self.reads.len() {
                "\n},\n"
            } else {
                "\n}\n"
            });
        }
        out.push_str("]\n}\n");
        out
    }
}

/// One read of the built-in profile.
///
/// Three quality values (11, 25 and 37, the bins of a recent patterned
/// flow-cell instrument). The fraction of Q37 calls starts at `start` and
/// falls by `slope` per cycle (never below one half); the other calls are
/// split 60 % Q25 and 40 % Q11. A Q11 call is followed by another one with
/// extra probability 0.7 and a Q25 call with extra probability 0.5, which
/// gives short runs of low-quality bases while keeping the per-cycle
/// fractions exactly as stated.
fn builtin_read(start: f64, slope: f64) -> ReadProfile {
    const STICK: [f64; 3] = [0.7, 0.5, 0.0];
    let marginal = |c: usize| -> [f64; 3] {
        let q37 = (start - slope * c as f64).max(0.5);
        let rest = 1.0 - q37;
        [0.4 * rest, 0.6 * rest, q37]
    };
    let mut cycles = Vec::with_capacity(BUILTIN_CYCLES);
    let mut transitions = Vec::with_capacity(BUILTIN_CYCLES - 1);
    for c in 0..BUILTIN_CYCLES {
        let now = marginal(c);
        cycles.push(now.to_vec());
        if c > 0 {
            let before = marginal(c - 1);
            let free: f64 = (0..3).map(|i| (1.0 - STICK[i]) * before[i]).sum();
            let fresh: Vec<f64> = (0..3)
                .map(|j| ((now[j] - STICK[j] * before[j]) / free).max(0.0))
                .collect();
            let matrix: Vec<Vec<f64>> = (0..3)
                .map(|i| {
                    (0..3)
                        .map(|j| (1.0 - STICK[i]) * fresh[j] + if i == j { STICK[i] } else { 0.0 })
                        .collect()
                })
                .collect();
            transitions.push(matrix);
        }
    }
    ReadProfile {
        qualities: vec![11, 25, 37],
        cycles,
        transitions: Some(transitions),
    }
}

/// The built-in profile: a recent instrument with binned qualities, high
/// quality with a slow decline along the read, read 2 slightly lower.
pub fn builtin_profile() -> Profile {
    Profile {
        format: PROFILE_FORMAT.to_owned(),
        version: PROFILE_VERSION,
        description: "built-in: binned qualities 11/25/37; Q37 fraction 96.5 % falling by 0.05 \
                      points per cycle in read 1 and 95.5 % falling by 0.07 points per cycle in \
                      read 2"
            .to_owned(),
        reads: vec![builtin_read(0.965, 0.0005), builtin_read(0.955, 0.0007)],
    }
}

#[derive(Clone, Copy, Debug)]
struct Entry {
    threshold: u32,
    alias: u32,
}

/// Number of bases changed in one read.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ReadErrors {
    pub substitutions: u32,
    pub n_calls: u32,
}

/// Draws the qualities of a read and applies the substitution errors.
#[derive(Clone, Debug)]
pub struct QualitySampler {
    states: usize,
    /// `states` for a first-order model, 0 for independent cycles: the
    /// current state times this value is the offset of the next row.
    prev_step: usize,
    /// Entries per cycle.
    cycle_step: usize,
    read_len: usize,
    quality_char: Vec<u8>,
    error_threshold: Vec<u64>,
    /// `[cycle][previous state][column]` alias-table entries.
    table: Vec<Entry>,
    /// Whether the read is longer than the profile.
    extended: bool,
}

fn alias_entries(weights: &[f64]) -> std::result::Result<Vec<Entry>, String> {
    let (prob, alias) = vose(weights)?;
    Ok(prob
        .iter()
        .zip(&alias)
        .map(|(&p, &a)| Entry {
            threshold: if p >= 1.0 {
                u32::MAX
            } else {
                (p * 4_294_967_296.0) as u32
            },
            alias: a,
        })
        .collect())
}

impl QualitySampler {
    /// Prepares the sampler of one read for reads of `read_len` bases.
    ///
    /// A read shorter than the profile uses its first cycles; a longer one
    /// repeats the last cycle of the profile.
    pub fn new(profile: &ReadProfile, read_len: usize) -> std::result::Result<Self, String> {
        profile.validate("quality profile")?;
        if read_len == 0 {
            return Err("the read length must be at least 1".into());
        }
        let s = profile.qualities.len();
        let markov = profile.transitions.is_some();
        let prev_dim = if markov { s } else { 1 };
        let n_cycles = profile.cycles.len();
        let mut table = Vec::with_capacity(read_len * prev_dim * s);
        for c in 0..read_len {
            let src = c.min(n_cycles - 1);
            let marginal = alias_entries(&profile.cycles[src])?;
            match &profile.transitions {
                Some(t) if c > 0 && n_cycles > 1 => {
                    let matrix = &t[src.max(1) - 1];
                    for row in matrix {
                        if row.iter().sum::<f64>() > 0.0 {
                            table.extend(alias_entries(row)?);
                        } else {
                            // previous quality never seen at this cycle
                            table.extend_from_slice(&marginal);
                        }
                    }
                }
                _ => {
                    for _ in 0..prev_dim {
                        table.extend_from_slice(&marginal);
                    }
                }
            }
        }
        Ok(QualitySampler {
            states: s,
            prev_step: if markov { s } else { 0 },
            cycle_step: prev_dim * s,
            read_len,
            quality_char: profile.qualities.iter().map(|q| q + 33).collect(),
            error_threshold: profile
                .qualities
                .iter()
                .map(|&q| threshold_u64(error_probability(q)))
                .collect(),
            table,
            extended: read_len > n_cycles,
        })
    }

    pub fn read_len(&self) -> usize {
        self.read_len
    }

    /// Whether the read is longer than the profile it was built from.
    pub fn extended(&self) -> bool {
        self.extended
    }

    /// Writes the sequenced read and its quality string.
    ///
    /// `template` is the error-free read. For every cycle a quality is drawn,
    /// then the base is replaced by one of the three other bases (each with
    /// probability one third) with probability `10^(-Q/10)`. If
    /// `n_threshold` is not zero, the base is instead called `N` with
    /// quality 2 with the corresponding probability.
    #[inline]
    pub fn apply(
        &self,
        rng: &mut Rng,
        template: &[u8],
        seq: &mut [u8],
        qual: &mut [u8],
        n_threshold: u64,
    ) -> ReadErrors {
        debug_assert_eq!(template.len(), self.read_len);
        let s = self.states;
        let s64 = s as u64;
        let mut errors = ReadErrors::default();
        let mut prev = 0usize;
        let rows = self.table.chunks_exact(self.cycle_step);
        for (((&base, out_base), out_qual), row) in template
            .iter()
            .zip(seq.iter_mut())
            .zip(qual.iter_mut())
            .zip(rows)
        {
            let r = rng.next_u64();
            let col = (((r >> 32) * s64) >> 32) as usize;
            let entry = row[prev + col];
            // Branch-free choice between the column and its alias: which one
            // is taken is random, so a branch here would be mispredicted
            // often.
            let keep = 0usize.wrapping_sub(((r as u32) < entry.threshold) as usize);
            let alias = entry.alias as usize;
            let state = alias ^ ((alias ^ col) & keep);
            prev = state * self.prev_step;
            let mut b = base;
            let mut q = self.quality_char[state];
            let wrong = rng.next_u64() < self.error_threshold[state];
            if n_threshold != 0 && rng.next_u64() < n_threshold {
                b = b'N';
                q = N_QUALITY + 33;
                errors.n_calls += 1;
            } else if wrong {
                let k = rng.below(3) as usize;
                let replacement = match base {
                    b'A' => b"CGT"[k],
                    b'C' => b"AGT"[k],
                    b'G' => b"ACT"[k],
                    b'T' => b"ACG"[k],
                    other => other,
                };
                if replacement != base {
                    b = replacement;
                    errors.substitutions += 1;
                }
            }
            *out_base = b;
            *out_qual = q;
        }
        errors
    }
}

/// Accumulates quality counts from the reads of one FASTQ file.
pub struct Learner {
    with_transitions: bool,
    /// `marginal[c][q]`.
    marginal: Vec<Vec<u64>>,
    /// `transitions[c - 1][p * 94 + q]`.
    transitions: Vec<Vec<u64>>,
    reads: u64,
}

const Q_STATES: usize = MAX_QUALITY as usize + 1;

impl Learner {
    pub fn new(with_transitions: bool) -> Self {
        Learner {
            with_transitions,
            marginal: Vec::new(),
            transitions: Vec::new(),
            reads: 0,
        }
    }

    pub fn reads(&self) -> u64 {
        self.reads
    }

    /// Adds the quality string of one read (Phred+33 characters).
    pub fn add(&mut self, quality: &[u8]) -> std::result::Result<(), String> {
        if quality.len() > MAX_PROFILE_CYCLES {
            return Err(format!(
                "a read has {} bases; reads longer than {MAX_PROFILE_CYCLES} bases are not supported",
                quality.len()
            ));
        }
        if let Some(&bad) = quality.iter().find(|&&c| !(33..=126).contains(&c)) {
            return Err(format!(
                "quality character with code {bad} is outside the Phred+33 range (33 to 126)"
            ));
        }
        while self.marginal.len() < quality.len() {
            self.marginal.push(vec![0; Q_STATES]);
        }
        if self.with_transitions {
            while self.transitions.len() + 1 < quality.len() {
                self.transitions.push(vec![0; Q_STATES * Q_STATES]);
            }
        }
        let mut prev = 0usize;
        for (c, &ch) in quality.iter().enumerate() {
            let q = (ch - 33) as usize;
            self.marginal[c][q] += 1;
            if self.with_transitions && c > 0 {
                self.transitions[c - 1][prev * Q_STATES + q] += 1;
            }
            prev = q;
        }
        self.reads += 1;
        Ok(())
    }

    /// The profile of the reads seen so far.
    pub fn finish(&self) -> std::result::Result<ReadProfile, String> {
        if self.reads == 0 || self.marginal.is_empty() {
            return Err("no reads with at least one base".into());
        }
        let used: Vec<usize> = (0..Q_STATES)
            .filter(|&q| self.marginal.iter().any(|m| m[q] > 0))
            .collect();
        let cycles: Vec<Vec<f64>> = self
            .marginal
            .iter()
            .map(|m| used.iter().map(|&q| m[q] as f64).collect())
            .collect();
        let transitions = self.with_transitions.then(|| {
            self.transitions
                .iter()
                .map(|t| {
                    used.iter()
                        .map(|&p| used.iter().map(|&q| t[p * Q_STATES + q] as f64).collect())
                        .collect()
                })
                .collect()
        });
        Ok(ReadProfile {
            qualities: used.iter().map(|&q| q as u8).collect(),
            cycles,
            transitions,
        })
    }
}

/// Learns the profile of one read from a FASTQ file (plain or gzip).
///
/// Uses the first `max_reads` reads (all of them if `max_reads` is 0).
/// Returns the profile and the number of reads used.
pub fn learn_from_fastq(
    path: &Path,
    max_reads: u64,
    with_transitions: bool,
) -> Result<(ReadProfile, u64)> {
    let reader = open_text(path)?;
    learn_from_reader(
        reader,
        &path.display().to_string(),
        max_reads,
        with_transitions,
    )
}

/// Like [`learn_from_fastq`], from any reader; `name` is used in messages.
pub fn learn_from_reader<R: BufRead>(
    mut reader: R,
    name: &str,
    max_reads: u64,
    with_transitions: bool,
) -> Result<(ReadProfile, u64)> {
    let mut learner = Learner::new(with_transitions);
    let mut lines: [Vec<u8>; 4] = Default::default();
    let mut line_no = 0u64;
    'records: loop {
        if max_reads != 0 && learner.reads() >= max_reads {
            break;
        }
        for (i, line) in lines.iter_mut().enumerate() {
            line.clear();
            let n = reader
                .read_until(b'\n', line)
                .context(|| format!("cannot read {name}"))?;
            if n == 0 {
                if i == 0 {
                    break 'records;
                }
                return Err(Error(format!(
                    "{name}: the file ends in the middle of a FASTQ record (line {})",
                    line_no + 1
                )));
            }
            line_no += 1;
            while matches!(line.last(), Some(b'\n' | b'\r')) {
                line.pop();
            }
        }
        let record_line = line_no - 3;
        if lines[0].first() != Some(&b'@') {
            return Err(Error(format!(
                "{name}: line {record_line}: expected a FASTQ header starting with '@'"
            )));
        }
        if lines[2].first() != Some(&b'+') {
            return Err(Error(format!(
                "{name}: line {}: expected the '+' separator line of a FASTQ record",
                record_line + 2
            )));
        }
        if lines[1].len() != lines[3].len() {
            return Err(Error(format!(
                "{name}: line {record_line}: the sequence has {} bases but the quality string has {} characters",
                lines[1].len(),
                lines[3].len()
            )));
        }
        learner
            .add(&lines[3])
            .context(|| format!("{name}: line {}", record_line + 3))?;
    }
    let reads = learner.reads();
    let profile = learner.finish().context(|| name.to_string())?;
    Ok((profile, reads))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_to_error_probability() {
        assert!((error_probability(0) - 1.0).abs() < 1e-15);
        assert!((error_probability(10) - 0.1).abs() < 1e-15);
        assert!((error_probability(20) - 0.01).abs() < 1e-15);
        assert!((error_probability(30) - 0.001).abs() < 1e-15);
        assert!((error_probability(37) - 1.9952623149688795e-4).abs() < 1e-15);
    }

    fn flat(quality: u8, cycles: usize) -> ReadProfile {
        ReadProfile {
            qualities: vec![quality],
            cycles: vec![vec![1.0]; cycles],
            transitions: None,
        }
    }

    #[test]
    fn substitution_rate_follows_quality() {
        for (q, n_reads) in [(10u8, 2_000usize), (20, 20_000), (30, 100_000)] {
            let sampler = QualitySampler::new(&flat(q, 100), 100).unwrap();
            let mut rng = Rng::for_stream(q as u64, 0);
            let template = b"ACGT".repeat(25);
            let mut seq = vec![0u8; 100];
            let mut qual = vec![0u8; 100];
            let mut subs = 0u64;
            let mut mismatches = 0u64;
            for _ in 0..n_reads {
                let e = sampler.apply(&mut rng, &template, &mut seq, &mut qual, 0);
                subs += e.substitutions as u64;
                assert_eq!(e.n_calls, 0);
                mismatches += seq.iter().zip(&template).filter(|(a, b)| a != b).count() as u64;
                assert!(qual.iter().all(|&c| c == q + 33));
                assert!(seq.iter().all(|b| b"ACGT".contains(b)));
            }
            assert_eq!(subs, mismatches);
            let bases = (n_reads * 100) as f64;
            let p = error_probability(q);
            let expected = bases * p;
            let sd = (bases * p * (1.0 - p)).sqrt();
            assert!(
                (subs as f64 - expected).abs() < 5.0 * sd,
                "Q{q}: {subs} substitutions, expected {expected}"
            );
        }
    }

    #[test]
    fn substitutions_are_uniform_over_the_other_bases() {
        let sampler = QualitySampler::new(&flat(0, 50), 50).unwrap();
        let mut rng = Rng::for_stream(2, 0);
        let template = vec![b'A'; 50];
        let mut seq = vec![0u8; 50];
        let mut qual = vec![0u8; 50];
        let mut counts = std::collections::BTreeMap::new();
        for _ in 0..2_000 {
            sampler.apply(&mut rng, &template, &mut seq, &mut qual, 0);
            for &b in &seq {
                *counts.entry(b).or_insert(0u32) += 1;
            }
        }
        // Q0: every base is substituted
        assert!(!counts.contains_key(&b'A'));
        for b in [b'C', b'G', b'T'] {
            let n = counts[&b] as f64;
            assert!((n - 33_333.3).abs() < 1_000.0, "{}: {n}", b as char);
        }
    }

    #[test]
    fn ambiguous_bases_are_kept() {
        let sampler = QualitySampler::new(&flat(0, 10), 10).unwrap();
        let mut rng = Rng::for_stream(2, 0);
        let template = vec![b'N'; 10];
        let mut seq = vec![0u8; 10];
        let mut qual = vec![0u8; 10];
        let e = sampler.apply(&mut rng, &template, &mut seq, &mut qual, 0);
        assert_eq!(seq, template);
        assert_eq!(e, ReadErrors::default());
    }

    #[test]
    fn n_calls() {
        let sampler = QualitySampler::new(&flat(93, 100), 100).unwrap();
        let mut rng = Rng::for_stream(4, 0);
        let template = b"ACGT".repeat(25);
        let mut seq = vec![0u8; 100];
        let mut qual = vec![0u8; 100];
        let mut ns = 0u64;
        for _ in 0..2_000 {
            let e = sampler.apply(
                &mut rng,
                &template,
                &mut seq,
                &mut qual,
                threshold_u64(0.05),
            );
            ns += e.n_calls as u64;
            for (i, &b) in seq.iter().enumerate() {
                if b == b'N' {
                    assert_eq!(qual[i], b'#');
                } else {
                    assert_eq!(b, template[i]);
                    assert_eq!(qual[i], 93 + 33);
                }
            }
        }
        assert!((ns as f64 - 10_000.0).abs() < 500.0, "{ns}");
    }

    #[test]
    fn independent_cycles_follow_the_marginals() {
        let profile = ReadProfile {
            qualities: vec![2, 20, 40],
            cycles: vec![
                vec![0.0, 1.0, 3.0],
                vec![1.0, 1.0, 2.0],
                vec![0.0, 0.0, 5.0],
            ],
            transitions: None,
        };
        let sampler = QualitySampler::new(&profile, 3).unwrap();
        let mut rng = Rng::for_stream(8, 0);
        let template = b"AAA".to_vec();
        let mut seq = vec![0u8; 3];
        let mut qual = vec![0u8; 3];
        let mut counts = [[0u32; 3]; 3];
        let n = 40_000;
        for _ in 0..n {
            sampler.apply(&mut rng, &template, &mut seq, &mut qual, 0);
            for c in 0..3 {
                let s = match qual[c] - 33 {
                    2 => 0,
                    20 => 1,
                    40 => 2,
                    other => panic!("unexpected quality {other}"),
                };
                counts[c][s] += 1;
            }
        }
        assert_eq!(counts[0][0], 0);
        assert!((counts[0][1] as f64 - 10_000.0).abs() < 500.0);
        assert!((counts[1][0] as f64 - 10_000.0).abs() < 500.0);
        assert!((counts[1][2] as f64 - 20_000.0).abs() < 600.0);
        assert_eq!(counts[2], [0, 0, n]);
    }

    #[test]
    fn markov_chain_follows_the_transitions() {
        // Two states that never change after the first cycle.
        let profile = ReadProfile {
            qualities: vec![10, 40],
            cycles: vec![vec![1.0, 3.0]; 20],
            transitions: Some(vec![vec![vec![1.0, 0.0], vec![0.0, 1.0]]; 19]),
        };
        let sampler = QualitySampler::new(&profile, 20).unwrap();
        let mut rng = Rng::for_stream(6, 0);
        let template = vec![b'C'; 20];
        let mut seq = vec![0u8; 20];
        let mut qual = vec![0u8; 20];
        let mut low = 0u32;
        let n = 20_000;
        for _ in 0..n {
            sampler.apply(&mut rng, &template, &mut seq, &mut qual, 0);
            assert!(
                qual.iter().all(|&c| c == qual[0]),
                "quality changed within a read"
            );
            if qual[0] == 10 + 33 {
                low += 1;
            }
        }
        assert!((low as f64 - 5_000.0).abs() < 350.0, "{low}");
    }

    #[test]
    fn longer_reads_repeat_the_last_cycle_and_shorter_ones_truncate() {
        let profile = ReadProfile {
            qualities: vec![10, 40],
            cycles: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            transitions: None,
        };
        let long = QualitySampler::new(&profile, 5).unwrap();
        assert!(long.extended());
        let mut rng = Rng::for_stream(6, 0);
        let mut seq = vec![0u8; 5];
        let mut qual = vec![0u8; 5];
        long.apply(&mut rng, b"NNNNN", &mut seq, &mut qual, 0);
        assert_eq!(qual, [73, 43, 43, 43, 43]);
        let short = QualitySampler::new(&profile, 1).unwrap();
        assert!(!short.extended());
        let mut qual = vec![0u8; 1];
        short.apply(&mut rng, b"N", &mut seq[..1], &mut qual, 0);
        assert_eq!(qual, [73]);
    }

    #[test]
    fn builtin_profile_is_valid_and_declines() {
        let p = builtin_profile();
        p.validate().unwrap();
        assert_eq!(p.reads.len(), 2);
        for read in &p.reads {
            assert_eq!(read.qualities, vec![11, 25, 37]);
            assert_eq!(read.cycles.len(), BUILTIN_CYCLES);
            for c in &read.cycles {
                assert!((c.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            }
            assert!(read.cycles[0][2] > read.cycles[99][2]);
            for (c, m) in read.transitions.as_ref().unwrap().iter().enumerate() {
                // rows are distributions, and the chain reproduces the marginals
                let mut next = [0.0f64; 3];
                for (i, row) in m.iter().enumerate() {
                    assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                    assert!(row.iter().all(|&w| w >= 0.0));
                    for (n, w) in next.iter_mut().zip(row) {
                        *n += read.cycles[c][i] * w;
                    }
                }
                for (n, expected) in next.iter().zip(&read.cycles[c + 1]) {
                    assert!((n - expected).abs() < 1e-12);
                }
            }
        }
        // read 2 is slightly worse than read 1
        assert!(p.reads[1].cycles[50][2] < p.reads[0].cycles[50][2]);
    }

    #[test]
    fn builtin_sampler_matches_its_marginals() {
        let p = builtin_profile();
        let sampler = QualitySampler::new(p.read(0), 100).unwrap();
        let mut rng = Rng::for_stream(12, 0);
        let template = vec![b'G'; 100];
        let mut seq = vec![0u8; 100];
        let mut qual = vec![0u8; 100];
        let n = 50_000;
        let mut q37_first = 0u32;
        let mut q37_last = 0u32;
        let mut low_pairs = 0u32;
        let mut low = 0u32;
        for _ in 0..n {
            sampler.apply(&mut rng, &template, &mut seq, &mut qual, 0);
            q37_first += (qual[0] == b'F') as u32;
            q37_last += (qual[99] == b'F') as u32;
            for w in qual.windows(2) {
                if w[0] == b',' {
                    low += 1;
                    low_pairs += (w[1] == b',') as u32;
                }
            }
        }
        let f = q37_first as f64 / n as f64;
        let l = q37_last as f64 / n as f64;
        assert!((f - 0.965).abs() < 0.004, "{f}");
        assert!((l - (0.965 - 0.0005 * 99.0)).abs() < 0.005, "{l}");
        // Q11 calls come in runs: the next base is Q11 far more often than by chance
        assert!(low_pairs as f64 / low as f64 > 0.6);
    }

    #[test]
    fn learning_and_json_round_trip() {
        let fastq = b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nII#I\n@r3\nAC\n+\n#I\n";
        let (read, n) = learn_from_reader(&fastq[..], "test", 0, true).unwrap();
        assert_eq!(n, 3);
        assert_eq!(read.qualities, vec![2, 40]);
        assert_eq!(
            read.cycles,
            vec![
                vec![1.0, 2.0],
                vec![0.0, 3.0],
                vec![1.0, 1.0],
                vec![0.0, 2.0]
            ]
        );
        let t = read.transitions.as_ref().unwrap();
        assert_eq!(t.len(), 3);
        assert_eq!(t[0], vec![vec![0.0, 1.0], vec![0.0, 2.0]]);
        assert_eq!(t[1], vec![vec![0.0, 0.0], vec![1.0, 1.0]]);
        assert_eq!(t[2], vec![vec![0.0, 1.0], vec![0.0, 1.0]]);

        let profile = Profile {
            format: PROFILE_FORMAT.to_owned(),
            version: PROFILE_VERSION,
            description: "a \"quoted\" description".to_owned(),
            reads: vec![read.clone(), read],
        };
        let json = profile.to_json();
        let back: Profile = serde_json::from_str(&json).unwrap();
        back.validate().unwrap();
        assert_eq!(back, profile);

        let (independent, _) = learn_from_reader(&fastq[..], "test", 2, false).unwrap();
        assert!(independent.transitions.is_none());
        assert_eq!(independent.cycles[0], vec![0.0, 2.0]);

        // the built-in profile, with fractional weights, also survives
        let builtin = builtin_profile();
        let back: Profile = serde_json::from_str(&builtin.to_json()).unwrap();
        assert_eq!(back, builtin);
    }

    #[test]
    fn malformed_fastq_is_rejected() {
        for bad in [
            &b"r1\nACGT\n+\nIIII\n"[..],
            &b"@r1\nACGT\n-\nIIII\n"[..],
            &b"@r1\nACGT\n+\nIII\n"[..],
            &b"@r1\nACGT\n+\n"[..],
            &b"@r1\nACGT\n+\nII I\n"[..],
            &b""[..],
        ] {
            assert!(learn_from_reader(bad, "test", 0, true).is_err());
        }
    }

    #[test]
    fn invalid_profiles_are_rejected() {
        let good = Profile {
            format: PROFILE_FORMAT.to_owned(),
            version: PROFILE_VERSION,
            description: String::new(),
            reads: vec![flat(30, 3)],
        };
        good.validate().unwrap();
        let mut p = good.clone();
        p.format = "something else".into();
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.version = 99;
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.reads.clear();
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.reads[0].qualities = vec![94];
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.reads[0].cycles[1] = vec![0.0];
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.reads[0].cycles[1] = vec![1.0, 2.0];
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.reads[0].transitions = Some(vec![vec![vec![1.0]]]);
        assert!(p.validate().is_err());
        let mut p = good.clone();
        p.reads[0].cycles[0] = vec![f64::NAN];
        assert!(p.validate().is_err());
        assert!(serde_json::from_str::<Profile>("{\"format\": 3}").is_err());
    }
}
