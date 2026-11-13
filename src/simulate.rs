//! The simulation itself: sampling fragments, writing reads, and the
//! multi-threaded driver.
//!
//! The fragments are produced in chunks of [`CHUNK_FRAGMENTS`]. Chunk `i`
//! uses its own random stream, derived from the seed and `i`, and the chunks
//! are written in order; the output is therefore the same for any number of
//! threads.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::{Condvar, Mutex, MutexGuard};
use std::time::Instant;

use flate2::Compression;
use flate2::write::GzEncoder;
use serde::Serialize;

use crate::abundance::read_abundances;
use crate::alias::AliasTable;
use crate::bias::{BiasSettings, Biases, accept, expected_unique_fraction};
use crate::error::{Context, Error, Result};
use crate::fasta::Transcriptome;
use crate::fragment::FragmentModel;
use crate::names::{MAX_FRAGMENTS, NameLayout, push_decimal};
use crate::quality::{Profile, QualitySampler, ReadErrors, builtin_profile};
use crate::reads::{AdapterFill, LibraryType, build_read, parse_nucleotides, resolve_adapter};
use crate::rng::{Rng, mix64, threshold_u64};

/// Number of fragments per chunk. Part of the definition of the output: do
/// not change it without changing the version of the program.
pub const CHUNK_FRAGMENTS: u64 = 4096;

/// Longest read that can be simulated.
pub const MAX_READ_LENGTH: usize = 10_000;

/// How the values of the abundance table are interpreted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbundanceType {
    /// Relative molar abundance (TPM or anything proportional to it): the
    /// expected number of fragments is proportional to abundance × effective
    /// length.
    Tpm,
    /// Expected number of fragments (read counts): the expected number of
    /// fragments is proportional to the value itself.
    Counts,
}

impl AbundanceType {
    pub fn name(self) -> &'static str {
        match self {
            AbundanceType::Tpm => "tpm",
            AbundanceType::Counts => "counts",
        }
    }
}

/// Everything about the reads that does not depend on the transcriptome.
#[derive(Clone, Debug)]
pub struct ReadParameters {
    pub num_fragments: u64,
    pub seed: u64,
    pub library_type: LibraryType,
    pub read_length: usize,
    pub fragment_mean: f64,
    pub fragment_sd: f64,
    pub fragment_min: usize,
    pub fragment_max: usize,
    /// Read-through sequence of read 1; may contain `{i7}` and `{i5}`.
    pub adapter1: String,
    /// Read-through sequence of read 2; may contain `{i7}` and `{i5}`.
    pub adapter2: String,
    pub index1: String,
    pub index2: String,
    pub adapter_fill: AdapterFill,
    /// Probability that a base is called `N`.
    pub n_rate: f64,
    pub instrument: String,
    pub run_number: u32,
    pub flowcell: String,
    pub lane: u32,
    pub per_read_truth: bool,
    /// Which candidate fragments the library keeps (see [`crate::bias`]).
    pub biases: BiasSettings,
    /// Number of different molecules in the library; `None` draws a new
    /// molecule for every read pair (no duplicates beyond chance).
    pub library_molecules: Option<u64>,
}

/// Candidates drawn in a row without one being kept, after which the
/// simulation stops: the biases leave almost nothing to keep.
pub const MAX_REJECTIONS: u64 = 100_000_000;

impl Default for ReadParameters {
    fn default() -> Self {
        ReadParameters {
            num_fragments: 1_000_000,
            seed: 1,
            library_type: LibraryType::Reverse,
            read_length: 100,
            fragment_mean: 180.0,
            fragment_sd: 55.0,
            fragment_min: 40,
            fragment_max: 600,
            adapter1: crate::reads::DEFAULT_ADAPTER1.to_owned(),
            adapter2: crate::reads::DEFAULT_ADAPTER2.to_owned(),
            index1: crate::reads::DEFAULT_INDEX1.to_owned(),
            index2: crate::reads::DEFAULT_INDEX2.to_owned(),
            adapter_fill: AdapterFill::PolyG,
            n_rate: 0.0,
            instrument: "A00001".to_owned(),
            run_number: 1,
            flowcell: "H7SIMDSXX".to_owned(),
            lane: 1,
            per_read_truth: false,
            biases: BiasSettings {
                gc_center: None,
                gc_width: 0.15,
                degradation: None,
                hexamer: 0.0,
            },
            library_molecules: None,
        }
    }
}

/// Counters describing what was simulated. All integers, so that adding up
/// the counters of several threads gives the same result in any order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stats {
    pub fragments: u64,
    pub length_sum: u64,
    pub length_square_sum: u128,
    pub length_min: u64,
    pub length_max: u64,
    /// Fragments shorter than the read: both reads run into the adapter.
    pub shorter_than_read: u64,
    /// Read 1 runs past the end of its adapter.
    pub read1_past_adapter: u64,
    /// Read 2 runs past the end of its adapter.
    pub read2_past_adapter: u64,
    /// The drawn length exceeded the transcript; the fragment is the whole
    /// transcript.
    pub whole_transcript: u64,
    /// Read 1 is on the transcript strand.
    pub read1_on_transcript_strand: u64,
    pub substitutions: [u64; 2],
    pub n_calls: [u64; 2],
    /// Candidate fragments drawn (equal to `fragments` without biases).
    pub candidates: u64,
}

impl Default for Stats {
    fn default() -> Self {
        Stats {
            fragments: 0,
            length_sum: 0,
            length_square_sum: 0,
            length_min: u64::MAX,
            length_max: 0,
            shorter_than_read: 0,
            read1_past_adapter: 0,
            read2_past_adapter: 0,
            whole_transcript: 0,
            read1_on_transcript_strand: 0,
            substitutions: [0; 2],
            n_calls: [0; 2],
            candidates: 0,
        }
    }
}

impl Stats {
    pub fn merge(&mut self, other: &Stats) {
        self.fragments += other.fragments;
        self.length_sum += other.length_sum;
        self.length_square_sum += other.length_square_sum;
        self.length_min = self.length_min.min(other.length_min);
        self.length_max = self.length_max.max(other.length_max);
        self.shorter_than_read += other.shorter_than_read;
        self.read1_past_adapter += other.read1_past_adapter;
        self.read2_past_adapter += other.read2_past_adapter;
        self.whole_transcript += other.whole_transcript;
        self.read1_on_transcript_strand += other.read1_on_transcript_strand;
        self.candidates += other.candidates;
        for i in 0..2 {
            self.substitutions[i] += other.substitutions[i];
            self.n_calls[i] += other.n_calls[i];
        }
    }

    pub fn mean_length(&self) -> f64 {
        if self.fragments == 0 {
            0.0
        } else {
            self.length_sum as f64 / self.fragments as f64
        }
    }

    pub fn sd_length(&self) -> f64 {
        if self.fragments == 0 {
            return 0.0;
        }
        let n = self.fragments as f64;
        let mean = self.length_sum as f64 / n;
        (self.length_square_sum as f64 / n - mean * mean)
            .max(0.0)
            .sqrt()
    }
}

/// The uncompressed text of one chunk.
#[derive(Default, Debug)]
pub struct ChunkText {
    pub read1: Vec<u8>,
    pub read2: Vec<u8>,
    /// Lines of the per-read truth table (empty unless requested).
    pub truth: Vec<u8>,
}

impl ChunkText {
    fn clear(&mut self) {
        self.read1.clear();
        self.read2.clear();
        self.truth.clear();
    }
}

/// A fully prepared simulation: transcript weights, fragment model, adapters,
/// quality samplers and read names.
pub struct Simulator<'a> {
    transcriptome: &'a Transcriptome,
    transcripts: AliasTable,
    fragment: FragmentModel,
    library_type: LibraryType,
    read_length: usize,
    adapter1: Vec<u8>,
    adapter2: Vec<u8>,
    fill: AdapterFill,
    quality1: QualitySampler,
    quality2: QualitySampler,
    n_threshold: u64,
    names: NameLayout,
    name_prefix: Vec<u8>,
    suffix1: Vec<u8>,
    suffix2: Vec<u8>,
    seed: u64,
    total: u64,
    per_read_truth: bool,
    /// Effective length of every transcript.
    effective_lengths: Vec<f64>,
    /// Probability that a fragment comes from each transcript.
    probabilities: Vec<f64>,
    biases: Biases,
    molecules: Option<u64>,
    /// Seed of the per-molecule random streams.
    molecule_seed: u64,
}

fn check_name_field(value: &str, what: &str) -> Result<()> {
    if value.is_empty() {
        return Err(Error(format!("the {what} must not be empty")));
    }
    if let Some(c) = value
        .chars()
        .find(|c| !(c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.')))
    {
        return Err(Error(format!(
            "the {what} '{value}' contains '{}'; only letters, digits, '-', '_' and '.' are allowed \
             because it becomes part of the read names",
            c.escape_default()
        )));
    }
    Ok(())
}

impl<'a> Simulator<'a> {
    /// Prepares a simulation.
    ///
    /// `abundances` has one value per transcript; `profile` is the quality
    /// profile to use.
    pub fn new(
        transcriptome: &'a Transcriptome,
        abundances: &[f64],
        abundance_type: AbundanceType,
        profile: &Profile,
        p: &ReadParameters,
    ) -> Result<Self> {
        if transcriptome.is_empty() {
            return Err(Error("the transcriptome is empty".into()));
        }
        if abundances.len() != transcriptome.len() {
            return Err(Error(format!(
                "{} abundances for {} transcripts",
                abundances.len(),
                transcriptome.len()
            )));
        }
        if p.num_fragments == 0 || p.num_fragments > MAX_FRAGMENTS {
            return Err(Error(format!(
                "the number of fragments must be between 1 and {MAX_FRAGMENTS}"
            )));
        }
        if p.read_length == 0 || p.read_length > MAX_READ_LENGTH {
            return Err(Error(format!(
                "the read length must be between 1 and {MAX_READ_LENGTH}"
            )));
        }
        if !(0.0..=1.0).contains(&p.n_rate) {
            return Err(Error(format!(
                "the N rate must be between 0 and 1 (got {})",
                p.n_rate
            )));
        }
        check_name_field(&p.instrument, "instrument name")?;
        check_name_field(&p.flowcell, "flow cell name")?;
        let biases = Biases::new(&p.biases).map_err(Error)?;
        if p.library_molecules == Some(0) {
            return Err(Error("the library must hold at least one molecule".into()));
        }

        let fragment = FragmentModel::new(
            p.fragment_mean,
            p.fragment_sd,
            p.fragment_min,
            p.fragment_max,
        )
        .map_err(Error)?;
        let effective_lengths: Vec<f64> = (0..transcriptome.len())
            .map(|i| fragment.effective_length(transcriptome.length(i)))
            .collect();
        if let Some(bad) = abundances.iter().find(|a| !a.is_finite() || **a < 0.0) {
            return Err(Error(format!(
                "abundance {bad} is not a finite, non-negative number"
            )));
        }
        let weights: Vec<f64> = match abundance_type {
            AbundanceType::Tpm => abundances
                .iter()
                .zip(&effective_lengths)
                .map(|(a, e)| a * e)
                .collect(),
            AbundanceType::Counts => abundances.to_vec(),
        };
        let total_weight: f64 = weights.iter().sum();
        if !total_weight.is_finite() {
            return Err(Error(
                "the abundances are too large to be added up; rescale them".into(),
            ));
        }
        if total_weight <= 0.0 {
            return Err(Error(
                "no transcript has a positive abundance: there is nothing to simulate".into(),
            ));
        }
        let transcripts = AliasTable::new(&weights).map_err(Error)?;
        let probabilities: Vec<f64> = weights.iter().map(|w| w / total_weight).collect();

        let index1 = parse_nucleotides(&p.index1, "index 1").map_err(Error)?;
        let index2 = parse_nucleotides(&p.index2, "index 2").map_err(Error)?;
        let adapter1 =
            resolve_adapter(&p.adapter1, &index1, &index2, "adapter 1").map_err(Error)?;
        let adapter2 =
            resolve_adapter(&p.adapter2, &index1, &index2, "adapter 2").map_err(Error)?;

        let quality1 = QualitySampler::new(profile.read(0), p.read_length).map_err(Error)?;
        let quality2 = QualitySampler::new(profile.read(1), p.read_length).map_err(Error)?;

        let name_prefix = format!(
            "@{}:{}:{}:{}:",
            p.instrument, p.run_number, p.flowcell, p.lane
        )
        .into_bytes();
        let index_field = match (index1.is_empty(), index2.is_empty()) {
            (false, false) => format!("{}+{}", p.index1.to_uppercase(), p.index2.to_uppercase()),
            (false, true) => p.index1.to_uppercase(),
            (true, false) => p.index2.to_uppercase(),
            // without an index, Illumina's software writes the sample number
            (true, true) => "1".to_owned(),
        };
        let suffix1 = format!(" 1:N:0:{index_field}\n").into_bytes();
        let suffix2 = format!(" 2:N:0:{index_field}\n").into_bytes();

        Ok(Simulator {
            transcriptome,
            transcripts,
            fragment,
            library_type: p.library_type,
            read_length: p.read_length,
            adapter1,
            adapter2,
            fill: p.adapter_fill,
            quality1,
            quality2,
            n_threshold: threshold_u64(p.n_rate),
            names: NameLayout::new(p.num_fragments),
            name_prefix,
            suffix1,
            suffix2,
            seed: p.seed,
            total: p.num_fragments,
            per_read_truth: p.per_read_truth,
            effective_lengths,
            probabilities,
            biases,
            molecules: p.library_molecules,
            molecule_seed: mix64(p.seed ^ 0x6D6F_6C65_6375_6C65),
        })
    }

    /// Draws one fragment of the library: transcript, start, length, whether
    /// it is the whole transcript, and whether read 1 is on the transcript
    /// strand. Candidates are drawn as without biases and kept with the
    /// probability the biases give.
    fn draw_fragment(&self, rng: &mut Rng, stats: &mut Stats) -> (usize, usize, usize, bool, bool) {
        let mut rejected = 0u64;
        loop {
            let t = self.transcripts.sample(rng);
            let transcript = self.transcriptome.sequence(t);
            let tlen = transcript.len();
            let (start, flen, whole) = loop {
                let l = self.fragment.sample(rng);
                let p = rng.below(tlen as u64) as usize;
                if l <= tlen {
                    if p + l <= tlen {
                        break (p, l, false);
                    }
                } else if p == 0 {
                    break (0, tlen, true);
                }
            };
            let read1_forward = self.library_type.read1_on_transcript_strand(rng);
            stats.candidates += 1;
            let keep = self.biases.keep_probability(
                &transcript[start..start + flen],
                start,
                tlen,
                read1_forward,
            );
            if accept(rng, keep) {
                return (t, start, flen, whole, read1_forward);
            }
            rejected += 1;
            assert!(
                rejected < MAX_REJECTIONS,
                "{MAX_REJECTIONS} candidate fragments in a row were rejected: the library biases \
                 leave almost nothing to keep (loosen --gc-bias-width or --degradation)"
            );
        }
    }

    pub fn total_fragments(&self) -> u64 {
        self.total
    }

    pub fn num_chunks(&self) -> u64 {
        self.total.div_ceil(CHUNK_FRAGMENTS)
    }

    pub fn fragment_model(&self) -> &FragmentModel {
        &self.fragment
    }

    pub fn effective_lengths(&self) -> &[f64] {
        &self.effective_lengths
    }

    /// Probability that a fragment comes from each transcript.
    pub fn probabilities(&self) -> &[f64] {
        &self.probabilities
    }

    /// The resolved read-through sequences of read 1 and read 2.
    pub fn adapters(&self) -> (&[u8], &[u8]) {
        (&self.adapter1, &self.adapter2)
    }

    /// Whether the reads are longer than the quality profile of read 1 or 2.
    pub fn profile_extended(&self) -> bool {
        self.quality1.extended() || self.quality2.extended()
    }

    /// Produces the text of chunk `chunk`, adds the fragments to the
    /// per-transcript `counts` and updates `stats`.
    pub fn generate_chunk(
        &self,
        chunk: u64,
        out: &mut ChunkText,
        counts: &mut [u64],
        stats: &mut Stats,
    ) {
        out.clear();
        let first = chunk * CHUNK_FRAGMENTS;
        debug_assert!(first < self.total);
        let n = (self.total - first).min(CHUNK_FRAGMENTS);
        let mut rng = Rng::for_stream(self.seed, chunk);
        let mut cursor = self.names.cursor(first);
        let read_len = self.read_length;
        let mut template1 = vec![0u8; read_len];
        let mut template2 = vec![0u8; read_len];
        let mut name: Vec<u8> = Vec::with_capacity(self.name_prefix.len() + 24);
        let record = self.name_prefix.len() + 24 + self.suffix1.len() + 2 * read_len + 4;
        out.read1.reserve(record * n as usize);
        out.read2.reserve(record * n as usize);

        let plain = !self.biases.is_active() && self.molecules.is_none();
        for _ in 0..n {
            let (t, start, flen, whole, read1_forward);
            if plain {
                // 1. transcript, 2. fragment length and start (rejection
                // sampling of a uniform (start, length) pair inside the
                // transcript).
                t = self.transcripts.sample(&mut rng);
                let transcript = self.transcriptome.sequence(t);
                let tlen = transcript.len();
                (start, flen, whole) = loop {
                    let l = self.fragment.sample(&mut rng);
                    let p = rng.below(tlen as u64) as usize;
                    if l <= tlen {
                        if p + l <= tlen {
                            break (p, l, false);
                        }
                    } else if p == 0 {
                        break (0, tlen, true);
                    }
                };
                stats.candidates += 1;
                let fragment = &transcript[start..start + flen];

                // 3. strand, 4. the two error-free reads.
                read1_forward = self.library_type.read1_on_transcript_strand(&mut rng);
                build_read(
                    &mut template1,
                    fragment,
                    read1_forward,
                    &self.adapter1,
                    self.fill,
                    &mut rng,
                );
                build_read(
                    &mut template2,
                    fragment,
                    !read1_forward,
                    &self.adapter2,
                    self.fill,
                    &mut rng,
                );
            } else {
                // A molecule of the library: drawn from its own stream when
                // the library holds a limited number of molecules, so that a
                // duplicate is the same fragment.
                let mut molecule_rng;
                let frng: &mut Rng = match self.molecules {
                    Some(m) => {
                        molecule_rng = Rng::for_stream(self.molecule_seed, rng.below(m));
                        &mut molecule_rng
                    }
                    None => &mut rng,
                };
                (t, start, flen, whole, read1_forward) = self.draw_fragment(frng, stats);
                let fragment = &self.transcriptome.sequence(t)[start..start + flen];
                build_read(
                    &mut template1,
                    fragment,
                    read1_forward,
                    &self.adapter1,
                    self.fill,
                    frng,
                );
                build_read(
                    &mut template2,
                    fragment,
                    !read1_forward,
                    &self.adapter2,
                    self.fill,
                    frng,
                );
            }
            counts[t] += 1;

            // 5. name, 6. qualities and errors.
            let (tile, x, y) = cursor.next_coordinates(&mut rng);
            name.clear();
            name.extend_from_slice(&self.name_prefix);
            push_decimal(&mut name, tile as u64);
            name.push(b':');
            push_decimal(&mut name, x as u64);
            name.push(b':');
            push_decimal(&mut name, y as u64);
            let e1 = write_record(
                &mut out.read1,
                &name,
                &self.suffix1,
                &template1,
                &self.quality1,
                &mut rng,
                self.n_threshold,
            );
            let e2 = write_record(
                &mut out.read2,
                &name,
                &self.suffix2,
                &template2,
                &self.quality2,
                &mut rng,
                self.n_threshold,
            );

            if self.per_read_truth {
                let line = &mut out.truth;
                line.extend_from_slice(&name[1..]);
                line.push(b'\t');
                line.extend_from_slice(self.transcriptome.id(t).as_bytes());
                line.push(b'\t');
                push_decimal(line, start as u64 + 1);
                line.push(b'\t');
                push_decimal(line, flen as u64);
                line.push(b'\t');
                line.push(if read1_forward { b'+' } else { b'-' });
                line.push(b'\n');
            }

            let flen = flen as u64;
            stats.fragments += 1;
            stats.length_sum += flen;
            stats.length_square_sum += (flen as u128) * (flen as u128);
            stats.length_min = stats.length_min.min(flen);
            stats.length_max = stats.length_max.max(flen);
            if (flen as usize) < read_len {
                stats.shorter_than_read += 1;
                if flen as usize + self.adapter1.len() < read_len {
                    stats.read1_past_adapter += 1;
                }
                if flen as usize + self.adapter2.len() < read_len {
                    stats.read2_past_adapter += 1;
                }
            }
            stats.whole_transcript += whole as u64;
            stats.read1_on_transcript_strand += read1_forward as u64;
            stats.substitutions[0] += e1.substitutions as u64;
            stats.substitutions[1] += e2.substitutions as u64;
            stats.n_calls[0] += e1.n_calls as u64;
            stats.n_calls[1] += e2.n_calls as u64;
        }
    }
}

/// Appends one FASTQ record.
#[inline]
fn write_record(
    buf: &mut Vec<u8>,
    name: &[u8],
    suffix: &[u8],
    template: &[u8],
    sampler: &QualitySampler,
    rng: &mut Rng,
    n_threshold: u64,
) -> ReadErrors {
    let r = template.len();
    buf.extend_from_slice(name);
    buf.extend_from_slice(suffix);
    let start = buf.len();
    // sequence, "\n+\n", qualities, "\n"
    buf.resize(start + 2 * r + 4, b'\n');
    let (seq, rest) = buf[start..].split_at_mut(r);
    rest[1] = b'+';
    let qual = &mut rest[3..3 + r];
    sampler.apply(rng, template, seq, qual, n_threshold)
}

// ---------------------------------------------------------------------------
// Running a simulation from files to files
// ---------------------------------------------------------------------------

/// A complete run of `simulate`.
#[derive(Clone, Debug)]
pub struct Config {
    pub transcripts: PathBuf,
    pub abundance: PathBuf,
    pub id_column: Option<String>,
    pub abundance_column: Option<String>,
    pub abundance_type: AbundanceType,
    pub output_prefix: String,
    pub quality_profile: Option<PathBuf>,
    pub reads: ReadParameters,
    pub threads: usize,
    /// Write gzip files (otherwise plain text).
    pub compress: bool,
    /// gzip level, 1 (fastest) to 9 (smallest).
    pub compression_level: u32,
    pub quiet: bool,
}

/// What a run produced.
#[derive(Clone, Debug)]
pub struct Outcome {
    pub stats: Stats,
    pub read1_path: String,
    pub read2_path: String,
    pub truth_path: String,
    pub report_path: String,
    pub per_read_truth_path: Option<String>,
    /// Seconds spent loading the inputs.
    pub load_seconds: f64,
    /// Seconds spent simulating and writing the reads.
    pub simulate_seconds: f64,
}

#[derive(Default)]
struct ChunkOut {
    read1: Vec<u8>,
    read2: Vec<u8>,
    truth: Vec<u8>,
}

/// Compresses `data` as one gzip member into `out` (whose previous content
/// is discarded, but whose allocation is reused).
fn gzip_member(data: &[u8], level: u32, mut out: Vec<u8>) -> Vec<u8> {
    out.clear();
    let mut encoder = GzEncoder::new(out, Compression::new(level));
    encoder
        .write_all(data)
        .expect("writing to memory cannot fail");
    encoder.finish().expect("writing to memory cannot fail")
}

/// Generates chunks; one per thread.
struct Worker<'s, 'a> {
    simulator: &'s Simulator<'a>,
    compression: Option<u32>,
    text: ChunkText,
    counts: Vec<u64>,
    stats: Stats,
}

impl<'s, 'a> Worker<'s, 'a> {
    fn new(simulator: &'s Simulator<'a>, compression: Option<u32>) -> Self {
        Worker {
            simulator,
            compression,
            text: ChunkText::default(),
            counts: vec![0; simulator.transcriptome.len()],
            stats: Stats::default(),
        }
    }

    /// Produces the bytes to write for one chunk. `recycled` is an earlier
    /// result that has been written; its buffers are used again.
    fn produce(&mut self, chunk: u64, mut recycled: ChunkOut) -> ChunkOut {
        self.simulator
            .generate_chunk(chunk, &mut self.text, &mut self.counts, &mut self.stats);
        match self.compression {
            Some(level) => ChunkOut {
                read1: gzip_member(&self.text.read1, level, recycled.read1),
                read2: gzip_member(&self.text.read2, level, recycled.read2),
                truth: if self.text.truth.is_empty() {
                    Vec::new()
                } else {
                    gzip_member(&self.text.truth, level, recycled.truth)
                },
            },
            None => {
                std::mem::swap(&mut self.text.read1, &mut recycled.read1);
                std::mem::swap(&mut self.text.read2, &mut recycled.read2);
                std::mem::swap(&mut self.text.truth, &mut recycled.truth);
                recycled
            }
        }
    }
}

struct Outputs {
    read1: (File, String),
    read2: (File, String),
    truth: Option<(File, String)>,
}

fn create(path: &str) -> Result<File> {
    File::create(path).context(|| format!("cannot create {path}"))
}

fn write_all(target: &mut (File, String), data: &[u8]) -> Result<()> {
    target
        .0
        .write_all(data)
        .context(|| format!("cannot write {}", target.1))
}

impl Outputs {
    fn write(&mut self, chunk: &ChunkOut) -> Result<()> {
        write_all(&mut self.read1, &chunk.read1)?;
        write_all(&mut self.read2, &chunk.read2)?;
        if let Some(t) = &mut self.truth {
            write_all(t, &chunk.truth)?;
        }
        Ok(())
    }
}

struct State {
    next_assign: u64,
    next_write: u64,
    done: BTreeMap<u64, ChunkOut>,
    /// Written chunks whose buffers can be used again.
    spare: Vec<ChunkOut>,
    abort: bool,
}

struct Shared {
    state: Mutex<State>,
    /// Signalled when a worker may take another chunk.
    work: Condvar,
    /// Signalled when a finished chunk was added.
    ready: Condvar,
}

impl Shared {
    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn abort(&self) {
        self.lock().abort = true;
        self.work.notify_all();
        self.ready.notify_all();
    }
}

/// Makes the other threads stop if this thread panics.
struct AbortOnPanic<'a>(&'a Shared);

impl Drop for AbortOnPanic<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.abort();
        }
    }
}

fn run_single(
    simulator: &Simulator<'_>,
    compression: Option<u32>,
    outputs: &mut Outputs,
) -> Result<(Vec<u64>, Stats)> {
    let mut worker = Worker::new(simulator, compression);
    let mut out = ChunkOut::default();
    for chunk in 0..simulator.num_chunks() {
        out = worker.produce(chunk, out);
        outputs.write(&out)?;
    }
    Ok((worker.counts, worker.stats))
}

fn run_parallel(
    simulator: &Simulator<'_>,
    threads: usize,
    compression: Option<u32>,
    outputs: &mut Outputs,
) -> Result<(Vec<u64>, Stats)> {
    let n_chunks = simulator.num_chunks();
    // Chunks that may be in progress or waiting to be written; this bounds
    // the memory used, whatever the number of fragments.
    let window = (threads as u64 * 4).max(8);
    let shared = Shared {
        state: Mutex::new(State {
            next_assign: 0,
            next_write: 0,
            done: BTreeMap::new(),
            spare: Vec::new(),
            abort: false,
        }),
        work: Condvar::new(),
        ready: Condvar::new(),
    };
    let shared = &shared;

    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(move || {
                    let _guard = AbortOnPanic(shared);
                    let mut worker = Worker::new(simulator, compression);
                    loop {
                        let chunk = {
                            let mut st = shared.lock();
                            loop {
                                if st.abort || st.next_assign >= n_chunks {
                                    break None;
                                }
                                if st.next_assign < st.next_write + window {
                                    let c = st.next_assign;
                                    st.next_assign += 1;
                                    break Some((c, st.spare.pop().unwrap_or_default()));
                                }
                                st = shared.work.wait(st).unwrap_or_else(|e| e.into_inner());
                            }
                        };
                        let Some((chunk, spare)) = chunk else { break };
                        let out = worker.produce(chunk, spare);
                        shared.lock().done.insert(chunk, out);
                        shared.ready.notify_all();
                    }
                    (worker.counts, worker.stats)
                })
            })
            .collect();

        let mut result: Result<()> = Ok(());
        for chunk in 0..n_chunks {
            let out = {
                let mut st = shared.lock();
                loop {
                    if st.abort {
                        break None;
                    }
                    if let Some(out) = st.done.remove(&chunk) {
                        break Some(out);
                    }
                    st = shared.ready.wait(st).unwrap_or_else(|e| e.into_inner());
                }
            };
            let Some(out) = out else {
                result = Err(Error("a worker thread stopped unexpectedly".into()));
                break;
            };
            if let Err(e) = outputs.write(&out) {
                result = Err(e);
                break;
            }
            {
                let mut st = shared.lock();
                st.next_write = chunk + 1;
                st.spare.push(out);
            }
            shared.work.notify_all();
        }
        if result.is_err() {
            shared.abort();
        }

        let mut counts = vec![0u64; simulator.transcriptome.len()];
        let mut stats = Stats::default();
        for handle in handles {
            match handle.join() {
                Ok((c, s)) => {
                    for (total, part) in counts.iter_mut().zip(&c) {
                        *total += part;
                    }
                    stats.merge(&s);
                }
                Err(_) => {
                    if result.is_ok() {
                        result = Err(Error("a worker thread stopped unexpectedly".into()));
                    }
                }
            }
        }
        result.map(|()| (counts, stats))
    })
}

#[derive(Serialize)]
struct ReportParameters<'a> {
    transcripts: String,
    abundance: String,
    id_column: &'a str,
    abundance_column: &'a str,
    abundance_type: &'a str,
    output_prefix: &'a str,
    num_fragments: u64,
    seed: u64,
    library_type: &'a str,
    read_length: usize,
    fragment_mean: f64,
    fragment_sd: f64,
    fragment_min: usize,
    fragment_max: usize,
    adapter1: &'a str,
    adapter2: &'a str,
    index1: &'a str,
    index2: &'a str,
    adapter_fill: &'a str,
    quality_profile: String,
    n_rate: f64,
    instrument: &'a str,
    run_number: u32,
    flowcell: &'a str,
    lane: u32,
    compression: &'a str,
    compression_level: Option<u32>,
    per_read_truth: bool,
    gc_bias_center: Option<f64>,
    gc_bias_width: Option<f64>,
    degradation: Option<f64>,
    hexamer_bias: f64,
    library_molecules: Option<u64>,
}

#[derive(Serialize)]
struct ReportInput {
    transcripts: usize,
    transcript_bases: usize,
    abundance_rows: usize,
    transcripts_with_positive_abundance: usize,
}

#[derive(Serialize)]
struct ReportRead {
    bases: u64,
    substitutions: u64,
    substitution_rate: f64,
    n_calls: u64,
}

#[derive(Serialize)]
struct ReportSummary {
    fragments: u64,
    chunks: u64,
    fragments_per_chunk: u64,
    transcripts_with_fragments: usize,
    fragment_length_mean: f64,
    fragment_length_sd: f64,
    fragment_length_min: u64,
    fragment_length_max: u64,
    fragment_length_model_mean: f64,
    fragments_shorter_than_read: u64,
    fraction_shorter_than_read: f64,
    fragments_with_read1_past_adapter: u64,
    fragments_with_read2_past_adapter: u64,
    whole_transcript_fragments: u64,
    read1_on_transcript_strand: u64,
    fraction_read1_on_transcript_strand: f64,
    read1: ReportRead,
    read2: ReportRead,
    /// Candidate fragments drawn; above `fragments` when biases reject some.
    candidate_fragments: u64,
    fraction_of_candidates_kept: f64,
    /// Expected share of read pairs that are the first read of their
    /// molecule (with --library-molecules).
    expected_unique_fraction: Option<f64>,
}

#[derive(Serialize)]
struct ReportOutputs<'a> {
    read1: &'a str,
    read2: &'a str,
    truth: &'a str,
    per_read_truth: Option<&'a str>,
}

#[derive(Serialize)]
struct Report<'a> {
    tool: &'a str,
    version: &'a str,
    command: &'a str,
    parameters: ReportParameters<'a>,
    resolved_adapter1: String,
    resolved_adapter2: String,
    input: ReportInput,
    summary: ReportSummary,
    outputs: ReportOutputs<'a>,
}

fn fraction(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        part as f64 / total as f64
    }
}

/// Runs a simulation: reads the inputs, writes the FASTQ files, the truth
/// table and the report.
pub fn run(config: &Config) -> Result<Outcome> {
    let log = |msg: &str| {
        if !config.quiet {
            eprintln!("[rnaseq-readsim] {msg}");
        }
    };
    if config.threads == 0 {
        return Err(Error("the number of threads must be at least 1".into()));
    }
    if config.compress && !(1..=9).contains(&config.compression_level) {
        return Err(Error(format!(
            "the compression level must be between 1 and 9 (got {})",
            config.compression_level
        )));
    }
    if config.output_prefix.is_empty() {
        return Err(Error("the output prefix must not be empty".into()));
    }

    let started = Instant::now();
    let profile = match &config.quality_profile {
        Some(path) => Profile::from_file(path)?,
        None => builtin_profile(),
    };
    let transcriptome = Transcriptome::from_fasta(&config.transcripts)?;
    log(&format!(
        "read {} transcripts ({} bases) from {}",
        transcriptome.len(),
        transcriptome.total_bases(),
        config.transcripts.display()
    ));
    let abundances = read_abundances(
        &config.abundance,
        &transcriptome,
        config.id_column.as_deref(),
        config.abundance_column.as_deref(),
    )?;
    let positive = abundances.values.iter().filter(|a| **a > 0.0).count();
    log(&format!(
        "read {} abundances from {} (columns '{}' and '{}'); {} transcripts have a positive abundance",
        abundances.rows,
        config.abundance.display(),
        abundances.id_column,
        abundances.abundance_column,
        positive
    ));
    let simulator = Simulator::new(
        &transcriptome,
        &abundances.values,
        config.abundance_type,
        &profile,
        &config.reads,
    )?;
    if simulator.profile_extended() {
        log(&format!(
            "warning: the quality profile has fewer cycles than the read length ({}); \
             its last cycle is repeated",
            config.reads.read_length
        ));
    }
    if config.quality_profile.is_some() && profile.reads.len() == 1 {
        log("note: the quality profile describes one read only; it is used for both reads");
    }
    let load_seconds = started.elapsed().as_secs_f64();

    let prefix = &config.output_prefix;
    let (fastq_ext, table_ext) = if config.compress {
        ("fastq.gz", "tsv.gz")
    } else {
        ("fastq", "tsv")
    };
    let read1_path = format!("{prefix}_1.{fastq_ext}");
    let read2_path = format!("{prefix}_2.{fastq_ext}");
    let truth_path = format!("{prefix}.truth.tsv");
    let report_path = format!("{prefix}.report.json");
    let per_read_truth_path = config
        .reads
        .per_read_truth
        .then(|| format!("{prefix}.reads.{table_ext}"));
    let compression = config.compress.then_some(config.compression_level);

    let mut outputs = Outputs {
        read1: (create(&read1_path)?, read1_path.clone()),
        read2: (create(&read2_path)?, read2_path.clone()),
        truth: match &per_read_truth_path {
            Some(p) => Some((create(p)?, p.clone())),
            None => None,
        },
    };
    if let Some(t) = &mut outputs.truth {
        let header = b"read_name\ttranscript_id\tstart\tfragment_length\tread1_strand\n";
        let data = match compression {
            Some(level) => gzip_member(header, level, Vec::new()),
            None => header.to_vec(),
        };
        write_all(t, &data)?;
    }

    log(&format!(
        "simulating {} fragments ({} library, 2 x {} bases) with {} thread{}",
        simulator.total_fragments(),
        config.reads.library_type.name(),
        config.reads.read_length,
        config.threads,
        if config.threads == 1 { "" } else { "s" }
    ));
    let sim_started = Instant::now();
    let (counts, stats) = if config.threads == 1 {
        run_single(&simulator, compression, &mut outputs)?
    } else {
        run_parallel(&simulator, config.threads, compression, &mut outputs)?
    };
    drop(outputs);
    let simulate_seconds = sim_started.elapsed().as_secs_f64();

    // Truth table.
    {
        let file = create(&truth_path)?;
        let mut w = BufWriter::new(file);
        let total = simulator.total_fragments() as f64;
        let write = |w: &mut BufWriter<File>| -> std::io::Result<()> {
            writeln!(
                w,
                "transcript_id\tlength\teffective_length\tabundance\texpected_fragments\tfragments"
            )?;
            for (i, count) in counts.iter().enumerate() {
                writeln!(
                    w,
                    "{}\t{}\t{:.3}\t{}\t{:.3}\t{}",
                    transcriptome.id(i),
                    transcriptome.length(i),
                    simulator.effective_lengths()[i],
                    abundances.values[i],
                    simulator.probabilities()[i] * total,
                    count
                )?;
            }
            w.flush()
        };
        write(&mut w).context(|| format!("cannot write {truth_path}"))?;
    }

    // Report.
    let read_bases = stats.fragments * config.reads.read_length as u64;
    let (adapter1, adapter2) = simulator.adapters();
    let report = Report {
        tool: env!("CARGO_PKG_NAME"),
        version: env!("CARGO_PKG_VERSION"),
        command: "simulate",
        parameters: ReportParameters {
            transcripts: config.transcripts.display().to_string(),
            abundance: config.abundance.display().to_string(),
            id_column: &abundances.id_column,
            abundance_column: &abundances.abundance_column,
            abundance_type: config.abundance_type.name(),
            output_prefix: prefix,
            num_fragments: config.reads.num_fragments,
            seed: config.reads.seed,
            library_type: config.reads.library_type.name(),
            read_length: config.reads.read_length,
            fragment_mean: config.reads.fragment_mean,
            fragment_sd: config.reads.fragment_sd,
            fragment_min: config.reads.fragment_min,
            fragment_max: config.reads.fragment_max,
            adapter1: &config.reads.adapter1,
            adapter2: &config.reads.adapter2,
            index1: &config.reads.index1,
            index2: &config.reads.index2,
            adapter_fill: config.reads.adapter_fill.name(),
            quality_profile: match &config.quality_profile {
                Some(p) => p.display().to_string(),
                None => "built-in".to_owned(),
            },
            n_rate: config.reads.n_rate,
            instrument: &config.reads.instrument,
            run_number: config.reads.run_number,
            flowcell: &config.reads.flowcell,
            lane: config.reads.lane,
            compression: if config.compress { "gzip" } else { "none" },
            compression_level: compression,
            per_read_truth: config.reads.per_read_truth,
            gc_bias_center: config.reads.biases.gc_center,
            gc_bias_width: config
                .reads
                .biases
                .gc_center
                .map(|_| config.reads.biases.gc_width),
            degradation: config.reads.biases.degradation,
            hexamer_bias: config.reads.biases.hexamer,
            library_molecules: config.reads.library_molecules,
        },
        resolved_adapter1: String::from_utf8_lossy(adapter1).into_owned(),
        resolved_adapter2: String::from_utf8_lossy(adapter2).into_owned(),
        input: ReportInput {
            transcripts: transcriptome.len(),
            transcript_bases: transcriptome.total_bases(),
            abundance_rows: abundances.rows,
            transcripts_with_positive_abundance: positive,
        },
        summary: ReportSummary {
            fragments: stats.fragments,
            chunks: simulator.num_chunks(),
            fragments_per_chunk: CHUNK_FRAGMENTS,
            transcripts_with_fragments: counts.iter().filter(|c| **c > 0).count(),
            fragment_length_mean: stats.mean_length(),
            fragment_length_sd: stats.sd_length(),
            fragment_length_min: stats.length_min,
            fragment_length_max: stats.length_max,
            fragment_length_model_mean: simulator.fragment_model().mean(),
            fragments_shorter_than_read: stats.shorter_than_read,
            fraction_shorter_than_read: fraction(stats.shorter_than_read, stats.fragments),
            fragments_with_read1_past_adapter: stats.read1_past_adapter,
            fragments_with_read2_past_adapter: stats.read2_past_adapter,
            whole_transcript_fragments: stats.whole_transcript,
            read1_on_transcript_strand: stats.read1_on_transcript_strand,
            fraction_read1_on_transcript_strand: fraction(
                stats.read1_on_transcript_strand,
                stats.fragments,
            ),
            read1: ReportRead {
                bases: read_bases,
                substitutions: stats.substitutions[0],
                substitution_rate: fraction(stats.substitutions[0], read_bases),
                n_calls: stats.n_calls[0],
            },
            read2: ReportRead {
                bases: read_bases,
                substitutions: stats.substitutions[1],
                substitution_rate: fraction(stats.substitutions[1], read_bases),
                n_calls: stats.n_calls[1],
            },
            candidate_fragments: stats.candidates,
            fraction_of_candidates_kept: fraction(stats.fragments, stats.candidates),
            expected_unique_fraction: config
                .reads
                .library_molecules
                .map(|m| expected_unique_fraction(stats.fragments, m)),
        },
        outputs: ReportOutputs {
            read1: &read1_path,
            read2: &read2_path,
            truth: &truth_path,
            per_read_truth: per_read_truth_path.as_deref(),
        },
    };
    let mut json = serde_json::to_string_pretty(&report)
        .context(|| "cannot serialise the report".to_owned())?;
    json.push('\n');
    std::fs::write(&report_path, json).context(|| format!("cannot write {report_path}"))?;

    log(&format!(
        "wrote {} read pairs in {:.2} s ({:.0} pairs/s); loading the inputs took {:.2} s",
        stats.fragments,
        simulate_seconds,
        stats.fragments as f64 / simulate_seconds.max(1e-9),
        load_seconds
    ));
    log(&format!(
        "mean fragment length {:.1}; {:.2} % of fragments are shorter than the read (adapter read-through)",
        stats.mean_length(),
        100.0 * fraction(stats.shorter_than_read, stats.fragments)
    ));

    Ok(Outcome {
        stats,
        read1_path,
        read2_path,
        truth_path,
        report_path,
        per_read_truth_path,
        load_seconds,
        simulate_seconds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_transcriptome() -> Transcriptome {
        // deterministic pseudo-random sequences of various lengths
        let mut rng = Rng::for_stream(2024, 0);
        let lengths = [3000usize, 1200, 800, 400, 250, 120, 60, 2000];
        Transcriptome::from_records(lengths.iter().enumerate().map(|(i, &n)| {
            let seq: Vec<u8> = (0..n)
                .map(|_| b"ACGT"[(rng.next_u64() >> 62) as usize])
                .collect();
            (format!("tx{}", i + 1), seq)
        }))
        .unwrap()
    }

    const TOY_ABUNDANCES: [f64; 8] = [10.0, 50.0, 5.0, 100.0, 30.0, 400.0, 3000.0, 0.0];

    fn flat_profile(quality: u8) -> Profile {
        Profile {
            format: crate::quality::PROFILE_FORMAT.to_owned(),
            version: crate::quality::PROFILE_VERSION,
            description: String::new(),
            reads: vec![crate::quality::ReadProfile {
                qualities: vec![quality],
                cycles: vec![vec![1.0]; 4],
                transitions: None,
            }],
        }
    }

    fn simulate_all(sim: &Simulator<'_>, n_tx: usize) -> (ChunkText, Vec<u64>, Stats) {
        let mut all = ChunkText::default();
        let mut text = ChunkText::default();
        let mut counts = vec![0u64; n_tx];
        let mut stats = Stats::default();
        for chunk in 0..sim.num_chunks() {
            sim.generate_chunk(chunk, &mut text, &mut counts, &mut stats);
            all.read1.extend_from_slice(&text.read1);
            all.read2.extend_from_slice(&text.read2);
            all.truth.extend_from_slice(&text.truth);
        }
        (all, counts, stats)
    }

    /// Checks observed counts against expected ones with a chi-square
    /// statistic (generous bound: the seed is fixed, so this cannot flake).
    fn assert_counts_match(counts: &[u64], expected: &[f64]) {
        let mut chi2 = 0.0;
        let mut df = 0usize;
        for (&c, &e) in counts.iter().zip(expected) {
            if e == 0.0 {
                assert_eq!(c, 0, "fragments from a transcript with weight 0");
            } else {
                chi2 += (c as f64 - e).powi(2) / e;
                df += 1;
            }
        }
        // 99.9th percentile of chi-square with 6 degrees of freedom is 22.5
        assert!(df <= 7, "test written for at most 7 categories");
        assert!(chi2 < 22.5, "chi-square {chi2} with {} categories", df);
    }

    #[test]
    fn fragments_are_proportional_to_abundance_times_effective_length() {
        let tx = toy_transcriptome();
        let params = ReadParameters {
            num_fragments: 300_000,
            seed: 7,
            ..Default::default()
        };
        let sim = Simulator::new(
            &tx,
            &TOY_ABUNDANCES,
            AbundanceType::Tpm,
            &builtin_profile(),
            &params,
        )
        .unwrap();
        let (_, counts, stats) = simulate_all(&sim, tx.len());
        assert_eq!(stats.fragments, 300_000);
        assert_eq!(counts.iter().sum::<u64>(), 300_000);
        // expected counts computed independently of the simulator's tables
        let model = FragmentModel::new(180.0, 55.0, 40, 600).unwrap();
        let weights: Vec<f64> = (0..tx.len())
            .map(|i| {
                let len = tx.length(i);
                let eff: f64 = (40..=600usize)
                    .map(|l| {
                        model.probability(l) * if l <= len { (len - l + 1) as f64 } else { 1.0 }
                    })
                    .sum();
                TOY_ABUNDANCES[i] * eff
            })
            .collect();
        let total: f64 = weights.iter().sum();
        let expected: Vec<f64> = weights.iter().map(|w| 300_000.0 * w / total).collect();
        assert_counts_match(&counts, &expected);
        for (p, e) in sim.probabilities().iter().zip(&expected) {
            assert!((p * 300_000.0 - e).abs() < 1e-6 * e.max(1.0));
        }
    }

    #[test]
    fn counts_mode_ignores_length() {
        let tx = toy_transcriptome();
        let params = ReadParameters {
            num_fragments: 200_000,
            seed: 3,
            ..Default::default()
        };
        let counts_in = [100.0, 100.0, 100.0, 100.0, 50.0, 50.0, 0.0, 0.0];
        let sim = Simulator::new(
            &tx,
            &counts_in,
            AbundanceType::Counts,
            &builtin_profile(),
            &params,
        )
        .unwrap();
        let (_, counts, _) = simulate_all(&sim, tx.len());
        let expected: Vec<f64> = counts_in.iter().map(|c| 200_000.0 * c / 500.0).collect();
        assert_counts_match(&counts, &expected);
    }

    #[test]
    fn start_positions_and_lengths_within_a_transcript() {
        // One transcript of 300 bases, fragments of exactly 100: starts must
        // be uniform on 1..=201. One of 60 bases: always the whole transcript.
        let tx = Transcriptome::from_records([
            ("long", "ACGT".repeat(75)),
            ("short", "ACGT".repeat(15)),
        ])
        .unwrap();
        let params = ReadParameters {
            num_fragments: 100_000,
            fragment_mean: 100.0,
            fragment_sd: 0.0,
            per_read_truth: true,
            read_length: 50,
            ..Default::default()
        };
        let sim = Simulator::new(
            &tx,
            &[1.0, 2000.0],
            AbundanceType::Tpm,
            &flat_profile(40),
            &params,
        )
        .unwrap();
        assert_eq!(sim.effective_lengths(), &[201.0, 1.0]);
        let (text, counts, stats) = simulate_all(&sim, 2);
        // expected share of the long transcript: 201 / (201 + 2000)
        let expected_long = 100_000.0 * 201.0 / 2201.0;
        assert!((counts[0] as f64 - expected_long).abs() < 5.0 * expected_long.sqrt());
        assert_eq!(stats.whole_transcript, counts[1]);
        let mut starts = vec![0u32; 202];
        for line in std::str::from_utf8(&text.truth).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            let start: usize = f[2].parse().unwrap();
            let len: usize = f[3].parse().unwrap();
            if f[1] == "long" {
                assert_eq!(len, 100);
                assert!((1..=201).contains(&start));
                starts[start] += 1;
            } else {
                assert_eq!((start, len), (1, 60));
            }
        }
        let per_start = counts[0] as f64 / 201.0;
        for (i, &n) in starts.iter().enumerate().skip(1) {
            assert!(
                (n as f64 - per_start).abs() < 6.0 * per_start.sqrt(),
                "start {i}: {n} fragments, expected about {per_start}"
            );
        }
    }

    /// Splits FASTQ text into (name, sequence, quality) records.
    fn records(text: &[u8]) -> Vec<(String, String, String)> {
        let lines: Vec<&str> = std::str::from_utf8(text).unwrap().lines().collect();
        assert_eq!(lines.len() % 4, 0);
        lines
            .chunks(4)
            .map(|c| {
                assert!(c[0].starts_with('@'));
                assert_eq!(c[2], "+");
                assert_eq!(c[1].len(), c[3].len());
                (c[0].to_owned(), c[1].to_owned(), c[3].to_owned())
            })
            .collect()
    }

    #[test]
    fn reads_match_the_transcripts_for_every_library_type() {
        let tx = toy_transcriptome();
        for library_type in [
            LibraryType::Reverse,
            LibraryType::Forward,
            LibraryType::Unstranded,
        ] {
            let params = ReadParameters {
                num_fragments: 6_000,
                library_type,
                read_length: 75,
                per_read_truth: true,
                fragment_mean: 120.0,
                fragment_sd: 60.0,
                fragment_min: 20,
                ..Default::default()
            };
            // Q93: no substitution in practice (5e-10 per base)
            let sim = Simulator::new(
                &tx,
                &TOY_ABUNDANCES,
                AbundanceType::Tpm,
                &flat_profile(93),
                &params,
            )
            .unwrap();
            let (adapter1, adapter2) = sim.adapters();
            let (text, _, stats) = simulate_all(&sim, tx.len());
            assert_eq!(stats.substitutions, [0, 0]);
            let r1 = records(&text.read1);
            let r2 = records(&text.read2);
            let truth: Vec<Vec<String>> = std::str::from_utf8(&text.truth)
                .unwrap()
                .lines()
                .map(|l| l.split('\t').map(str::to_owned).collect())
                .collect();
            assert_eq!(r1.len(), 6_000);
            assert_eq!(r2.len(), 6_000);
            assert_eq!(truth.len(), 6_000);
            let mut forward = 0u64;
            let mut short = 0u64;
            let mut past_adapter = 0u64;
            for ((a, b), t) in r1.iter().zip(&r2).zip(&truth) {
                // names: same up to the read number; truth has no '@'
                let (name1, tag1) = a.0.split_once(' ').unwrap();
                let (name2, tag2) = b.0.split_once(' ').unwrap();
                assert_eq!(name1, name2);
                assert_eq!(&name1[1..], t[0]);
                assert_eq!(tag1, "1:N:0:CCGCGGTT+CTAGCGCT");
                assert_eq!(tag2, "2:N:0:CCGCGGTT+CTAGCGCT");
                let id = tx.ids().iter().position(|i| *i == t[1]).unwrap();
                let start: usize = t[2].parse::<usize>().unwrap() - 1;
                let len: usize = t[3].parse().unwrap();
                let fragment = &tx.sequence(id)[start..start + len];
                let read1_forward = t[4] == "+";
                match library_type {
                    LibraryType::Reverse => assert!(!read1_forward),
                    LibraryType::Forward => assert!(read1_forward),
                    LibraryType::Unstranded => {}
                }
                forward += read1_forward as u64;
                let rc = crate::reads::reverse_complement(fragment);
                let (insert1, insert2): (&[u8], &[u8]) = if read1_forward {
                    (fragment, &rc)
                } else {
                    (&rc, fragment)
                };
                for (read, insert, adapter) in
                    [(&a.1, insert1, adapter1), (&b.1, insert2, adapter2)]
                {
                    let read = read.as_bytes();
                    assert_eq!(read.len(), 75);
                    let n = len.min(75);
                    assert_eq!(&read[..n], &insert[..n], "insert part of {}", a.0);
                    let k = adapter.len().min(75 - n);
                    assert_eq!(&read[n..n + k], &adapter[..k], "adapter part of {}", a.0);
                    assert!(
                        read[n + k..].iter().all(|&b| b == b'G'),
                        "filler of {}",
                        a.0
                    );
                }
                short += (len < 75) as u64;
                past_adapter += (len + adapter1.len() < 75) as u64;
            }
            assert!(short > 500, "the test needs short fragments, got {short}");
            assert_eq!(stats.shorter_than_read, short);
            assert_eq!(stats.read1_past_adapter, past_adapter);
            assert_eq!(stats.read1_on_transcript_strand, forward);
            if library_type == LibraryType::Unstranded {
                assert!((2_700..3_300).contains(&forward), "{forward}");
            }
        }
    }

    #[test]
    fn chunks_do_not_depend_on_order_or_on_each_other() {
        let tx = toy_transcriptome();
        let params = ReadParameters {
            num_fragments: 3 * CHUNK_FRAGMENTS + 100,
            per_read_truth: true,
            ..Default::default()
        };
        let sim = Simulator::new(
            &tx,
            &TOY_ABUNDANCES,
            AbundanceType::Tpm,
            &builtin_profile(),
            &params,
        )
        .unwrap();
        assert_eq!(sim.num_chunks(), 4);
        let mut counts = vec![0u64; tx.len()];
        let mut stats = Stats::default();
        let mut forward: Vec<ChunkText> = Vec::new();
        for chunk in 0..4 {
            let mut text = ChunkText::default();
            sim.generate_chunk(chunk, &mut text, &mut counts, &mut stats);
            forward.push(text);
        }
        let mut counts2 = vec![0u64; tx.len()];
        let mut stats2 = Stats::default();
        for chunk in (0..4).rev() {
            let mut text = ChunkText::default();
            sim.generate_chunk(chunk, &mut text, &mut counts2, &mut stats2);
            assert_eq!(text.read1, forward[chunk as usize].read1);
            assert_eq!(text.read2, forward[chunk as usize].read2);
            assert_eq!(text.truth, forward[chunk as usize].truth);
        }
        assert_eq!(counts, counts2);
        assert_eq!(stats, stats2);
        // the last chunk is short, and chunks differ from each other
        assert_eq!(records(&forward[3].read1).len(), 100);
        assert_ne!(forward[0].read1, forward[1].read1);
        // a different seed gives different reads
        let other = Simulator::new(
            &tx,
            &TOY_ABUNDANCES,
            AbundanceType::Tpm,
            &builtin_profile(),
            &ReadParameters { seed: 2, ..params },
        )
        .unwrap();
        let mut text = ChunkText::default();
        other.generate_chunk(0, &mut text, &mut counts, &mut stats);
        assert_ne!(text.read1, forward[0].read1);
    }

    #[test]
    fn invalid_parameters_are_reported() {
        let tx = toy_transcriptome();
        let build = |f: &dyn Fn(&mut ReadParameters), abundances: &[f64]| {
            let mut p = ReadParameters::default();
            f(&mut p);
            Simulator::new(&tx, abundances, AbundanceType::Tpm, &builtin_profile(), &p).map(|_| ())
        };
        assert!(build(&|_| {}, &TOY_ABUNDANCES).is_ok());
        assert!(build(&|p| p.num_fragments = 0, &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.read_length = 0, &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.n_rate = 1.5, &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.fragment_min = 0, &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.fragment_max = 10, &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.index1 = "ACGX".into(), &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.adapter1 = "ACGT{i6}".into(), &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.instrument = "A 1".into(), &TOY_ABUNDANCES).is_err());
        assert!(build(&|p| p.flowcell = "".into(), &TOY_ABUNDANCES).is_err());
        assert!(build(&|_| {}, &[0.0; 8]).is_err());
        assert!(build(&|_| {}, &[1.0; 3]).is_err());
        let mut negative = TOY_ABUNDANCES;
        negative[0] = -1.0;
        assert!(build(&|_| {}, &negative).is_err());
        // single index and no index in the read headers
        let header = |i1: &str, i2: &str| {
            let p = ReadParameters {
                num_fragments: 1,
                index1: i1.into(),
                index2: i2.into(),
                ..Default::default()
            };
            let sim = Simulator::new(
                &tx,
                &TOY_ABUNDANCES,
                AbundanceType::Tpm,
                &builtin_profile(),
                &p,
            )
            .unwrap();
            let (text, _, _) = simulate_all(&sim, tx.len());
            records(&text.read2)[0]
                .0
                .split_once(' ')
                .unwrap()
                .1
                .to_owned()
        };
        assert_eq!(header("ATCACG", ""), "2:N:0:ATCACG");
        assert_eq!(header("", ""), "2:N:0:1");
        assert_eq!(header("acgt", "ttaa"), "2:N:0:ACGT+TTAA");
    }

    /// Time of the generation alone (no compression, no files):
    /// `cargo test --release -- --ignored --nocapture throughput`.
    #[test]
    #[ignore = "benchmark"]
    fn throughput() {
        let tx = toy_transcriptome();
        let params = ReadParameters {
            num_fragments: 2_000_000,
            ..Default::default()
        };
        let sim = Simulator::new(
            &tx,
            &TOY_ABUNDANCES,
            AbundanceType::Tpm,
            &builtin_profile(),
            &params,
        )
        .unwrap();
        let mut text = ChunkText::default();
        let mut counts = vec![0u64; tx.len()];
        let mut stats = Stats::default();
        let started = Instant::now();
        for chunk in 0..sim.num_chunks() {
            sim.generate_chunk(chunk, &mut text, &mut counts, &mut stats);
        }
        let seconds = started.elapsed().as_secs_f64();
        eprintln!(
            "{} pairs in {seconds:.2} s: {:.0} ns per pair",
            stats.fragments,
            1e9 * seconds / stats.fragments as f64
        );
    }
}
