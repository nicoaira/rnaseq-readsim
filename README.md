# rnaseq-readsim

A fast, reproducible simulator of Illumina paired-end bulk RNA-seq reads.

Given a transcriptome (FASTA) and an abundance table, `rnaseq-readsim` writes
a pair of FASTQ files that real tools treat like the output of a sequencer:
stranded or unstranded reads, adapter read-through when the fragment is
shorter than the read, per-cycle base qualities with matching substitution
errors, and read names in the instrument's format. It also writes what was
simulated (the "truth") so that the results of a pipeline can be checked.

- One binary, no run-time dependencies; multi-threaded.
- The same inputs, options and seed give **byte-identical** files for any
  number of threads.
- Memory use does not grow with the number of reads.
- Millions of read pairs per second (see [Performance](#performance)).

It is meant for teaching, for testing pipelines and for benchmarks in which a
simple, fully documented model is an advantage. It does not try to reproduce
the biases of real libraries; see [What is not simulated](#what-is-not-simulated).

## Contents

- [Build](#build)
- [Quick start](#quick-start)
- [Output files](#output-files)
- [Options of `simulate`](#options-of-simulate)
- [Quality profiles and the `profile` command](#quality-profiles-and-the-profile-command)
- [The model](#the-model)
- [What is not simulated](#what-is-not-simulated)
- [Adapter sequences](#adapter-sequences)
- [Reproducibility](#reproducibility)
- [Performance](#performance)
- [Validation with real tools](#validation-with-real-tools)
- [Exit status](#exit-status)
- [Development](#development)
- [Licence](#licence)

## Build

A stable Rust toolchain is all that is needed (developed with Rust 1.95;
the code uses the 2024 edition, so Rust 1.85 or later is required).

```bash
git clone <repository URL>
cd rnaseq-readsim
cargo build --release
# the program is target/release/rnaseq-readsim
```

or `cargo install --path .` to put it on your `PATH`. The dependencies are
`clap` (command line), `flate2` with its pure-Rust `zlib-rs` backend (gzip),
and `serde`/`serde_json` (profile and report files). No C compiler is needed.

## Quick start

```bash
# 2 million read pairs, reverse-stranded, 2 x 100 bases, from TPM values
rnaseq-readsim simulate \
    --transcripts gencode.transcripts.fa.gz \
    --abundance abundance.tsv \
    --output-prefix sample1 \
    --num-fragments 2M \
    --seed 1
```

`abundance.tsv` is a tab-separated table with a header:

```text
transcript_id	tpm
ENST00000456328.2	12.5
ENST00000450305.2	0
```

A Salmon `quant.sf`, an RSEM `*.isoforms.results` or a kallisto
`abundance.tsv` can be given directly: the identifier column is the first of
`transcript_id`, `Name`, `target_id` that exists, and the abundance column is
`tpm` in any case (`TPM`, `tpm`). Other columns are chosen with `--id-column`
and `--abundance-column`.

The test data of this repository can be used to try the program:

```bash
rnaseq-readsim simulate -t tests/data/transcripts.fa -a tests/data/abundance.tsv \
    -o toy -n 10000
```

## Output files

| File | Content |
|---|---|
| `<prefix>_1.fastq.gz`, `<prefix>_2.fastq.gz` | Read 1 and read 2, Phred+33, gzip (plain `.fastq` with `--uncompressed`). The gzip files consist of several gzip members, like BGZF files; every common tool reads them. |
| `<prefix>.truth.tsv` | One row per transcript of the FASTA file, in FASTA order: `transcript_id`, `length`, `effective_length`, `abundance` (as given), `expected_fragments` (the model's expectation), `fragments` (simulated). |
| `<prefix>.report.json` | Program version, every parameter including the seed, the resolved adapter sequences, the size of the inputs and summary numbers: mean, standard deviation, minimum and maximum of the fragment lengths, fragments shorter than the read (adapter read-through), reads that ran past the adapter, whole-transcript fragments, strand counts, substitutions and `N` calls per read. |
| `<prefix>.reads.tsv.gz` | Only with `--per-read-truth`: `read_name`, `transcript_id`, `start` (1-based position of the fragment in the transcript), `fragment_length`, `read1_strand` (`+` if read 1 is on the transcript strand). Rows are in the order of the FASTQ files. |

The truth is never written into the read names. A read looks like this:

```text
@A00001:1:H7SIMDSXX:1:1101:21308:1000 1:N:0:CCGCGGTT+CTAGCGCT
CGATATTGGCCACTAACTGCGGTCCCCGAGTGGGAGAATTATCAAATCGCCTGTCCAATGCTCGGATCACCCGTGCACCCACCGGGCACCAGTCAAACAC
+
FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```

## Options of `simulate`

`rnaseq-readsim simulate --help` prints the same information.

### Input

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--transcripts <FASTA>` | required | Transcript sequences, plain or gzip (recognised by content, not by name). The identifier is the header up to the first white space or `\|`, so `>ENST00000832824.1\|ENSG00000290825.2\|…` and `>ENST00000832824.1` both give `ENST00000832824.1`. Lower case is accepted; `U` is read as `T`; other IUPAC codes become `N`. Duplicate identifiers and empty sequences are errors. |
| `-a`, `--abundance <TSV>` | required | Tab-separated table with a header line, plain or gzip. Empty lines and lines starting with `#` are skipped. Identifiers are cut like FASTA headers. Transcripts that are not in the table get abundance 0. Identifiers that are not in the FASTA file are an error (the first five are listed), as are duplicated identifiers and negative or non-numeric values. |
| `--id-column <NAME>` | first of `transcript_id`, `Name`, `target_id` | Column with the transcript identifiers. |
| `--abundance-column <NAME>` | `tpm` (any case) | Column with the abundances. |
| `--abundance-type <tpm\|counts>` | `tpm` | `tpm`: the values are relative molar abundances (TPM, or anything proportional to the number of molecules); the expected number of fragments is proportional to abundance × effective length. `counts`: the values are expected numbers of fragments (for example Salmon's `NumReads`); the expected number of fragments is proportional to the value. |

### Output

| Option | Default | Meaning |
|---|---|---|
| `-o`, `--output-prefix <PREFIX>` | required | Prefix of all output files. May contain directories, which must exist. Existing files are overwritten. |
| `--uncompressed` | off | Write plain FASTQ (and a plain per-read table). |
| `--compression-level <1-9>` | 4 | gzip level; 1 is the fastest, 9 gives the smallest files. |
| `--per-read-truth` | off | Also write `<prefix>.reads.tsv.gz`. |

### Library

| Option | Default | Meaning |
|---|---|---|
| `-n`, `--num-fragments <N>` | 1000000 | Number of fragments, that is, read pairs. `k`, `M` and `G` suffixes are accepted (`500k`, `2M`, `1.5M`). At most 50 billion. |
| `-l`, `--library-type <reverse\|forward\|unstranded>` | `reverse` | `reverse`: read 1 is the reverse complement of the transcript and read 2 is on the transcript strand (dUTP protocols, "fr-firststrand"; Salmon calls it ISR, featureCounts `-s 2`). `forward`: read 1 is on the transcript strand (ISF, `-s 1`). `unstranded`: each fragment is `forward` or `reverse` with probability 1/2 (IU). |
| `-r`, `--read-length <BASES>` | 100 | Length of both reads (1 to 100,000). |
| `--fragment-mean <BASES>` | 180 | Mean of the normal fragment-length distribution before truncation. |
| `--fragment-sd <BASES>` | 55 | Its standard deviation; 0 gives a single length. |
| `--fragment-min <BASES>` | 40 | Shortest fragment. |
| `--fragment-max <BASES>` | 600 | Longest fragment. |

### Adapters

| Option | Default | Meaning |
|---|---|---|
| `--adapter1 <SEQUENCE>` | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC{i7}ATCTCGTATGCCGTCTTCTGCTTG` | What read 1 contains after the last base of a fragment shorter than the read. `{i7}` and `{i5}` are replaced by `--index1` and `--index2`. An empty string means no adapter (the read continues with the filler). |
| `--adapter2 <SEQUENCE>` | `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT{i5}GTGTAGATCTCGGTGGTCGCCGTATCATT` | The same for read 2. |
| `--index1 <SEQUENCE>` | `CCGCGGTT` | Index 1 (i7). Also written in the read headers. May be empty. |
| `--index2 <SEQUENCE>` | `CTAGCGCT` | Index 2 (i5) in the orientation in which read 2 reads through it. Written in the read headers after a `+`. May be empty (single-index library). |
| `--adapter-fill <polyG\|polyA\|random>` | `polyG` | Bases of a read past the end of the adapter. `polyG` is what two-colour instruments (NextSeq, NovaSeq) report when there is no signal; `polyA` or `random` are closer to four-colour instruments. |

See [Adapter sequences](#adapter-sequences) for the origin of the defaults and
for single-index and Nextera libraries.

### Qualities and errors

| Option | Default | Meaning |
|---|---|---|
| `-q`, `--quality-profile <JSON>` | built-in | Profile written by `rnaseq-readsim profile`. |
| `--n-rate <PROBABILITY>` | 0 | Probability that a base is called `N` with quality 2 (`#`), independently for every base. |

### Read names

| Option | Default | Meaning |
|---|---|---|
| `--instrument <NAME>` | `A00001` | Instrument field. |
| `--run-number <NUMBER>` | 1 | Run field. |
| `--flowcell <NAME>` | `H7SIMDSXX` | Flow cell field. |
| `--lane <NUMBER>` | 1 | Lane field. |

Names have the form `@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>`
followed by ` 1:N:0:<index1>+<index2>` for read 1 and ` 2:N:0:…` for read 2
(only `<index1>` for a single-index library). Tile numbers (1101 … 1178,
1201 … 1678, 2101 … 2678: two surfaces, six swaths of 78 tiles) increase
through the file; within a tile the y coordinate increases in discrete rows
and x increases within a row, as in a FASTQ file of a patterned flow cell.
The fragments are spread evenly over the 936 tiles and every name is unique.

### Run

| Option | Default | Meaning |
|---|---|---|
| `-s`, `--seed <SEED>` | 1 | Seed of the random number generator (0 to 2^64 − 1). |
| `-p`, `--threads <THREADS>` | all cores | Number of threads. The output does not depend on it. With 1, everything runs on a single thread; otherwise there are that many worker threads and the main thread only writes. |
| `--quiet` | off | No progress messages on standard error. |

## Quality profiles and the `profile` command

A quality profile gives, for read 1 and read 2 and for every cycle, the
distribution of the Phred quality value, and optionally the distribution
conditional on the quality of the previous cycle (a first-order Markov chain:
a low-quality base is more likely to be followed by another one).

**Built-in profile.** Used when `--quality-profile` is not given. It is
synthetic, not learned from a run, and imitates a recent patterned-flow-cell
instrument with binned qualities: only Q11 (`,`), Q25 (`:`) and Q37 (`F`)
occur. In read 1 the fraction of Q37 calls is 96.5 % at the first cycle and
falls by 0.05 percentage points per cycle (91.55 % at cycle 100); in read 2 it
starts at 95.5 % and falls by 0.07 points per cycle (88.57 % at cycle 100);
never below 50 %. The other calls are 60 % Q25 and 40 % Q11. A Q11 call
repeats itself with an extra probability of 0.7 and a Q25 call with 0.5, which
produces short runs of low-quality bases while keeping the per-cycle fractions
exactly as stated. The resulting substitution rate is about 0.2 % in read 1
and 0.3 % in read 2 for 100 cycles. The profile is defined for 300 cycles.

**Learning a profile from real reads.**

```bash
rnaseq-readsim profile -1 real_R1.fastq.gz -2 real_R2.fastq.gz -o profile.json
rnaseq-readsim simulate ... --quality-profile profile.json
```

| Option | Default | Meaning |
|---|---|---|
| `-1`, `--read1 <FASTQ>` | required | FASTQ file of read 1 (or of single-end reads), plain or gzip, Phred+33. |
| `-2`, `--read2 <FASTQ>` | none | FASTQ file of read 2. Without it the profile describes one read, which `simulate` uses for both reads. |
| `-o`, `--output <JSON>` | required | Profile file to write. |
| `--max-reads <N>` | 1000000 | Number of reads used from the start of each file; 0 uses all of them. |
| `--independent` | off | Record only the per-cycle distributions, without the dependence on the previous cycle. |
| `--quiet` | off | No progress messages. |

The profile is a JSON file with one cycle per line: for each read the list of
quality values that occur (`qualities`), the count of each value at each cycle
(`cycles`) and, unless `--independent` was given, the count of each value at a
cycle for each value at the previous cycle (`transitions`). A profile of
binned qualities takes a few kilobytes, one with about 40 distinct values and
60 cycles about 0.4 MB. Counts or frequencies are both accepted, so a profile
can also be written by hand.

When reads are shorter than the profile, its first cycles are used; when they
are longer, the last cycle is repeated and a warning is printed. If a quality
value was never seen as the previous value at some cycle, the per-cycle
distribution is used instead.

## The model

Notation: transcript *t* has length *len(t)* and abundance *a(t)* from the
table. All random choices of different fragments are independent.

**Fragment-length distribution.** *g(L)* is a normal distribution restricted
to the integers from `--fragment-min` to `--fragment-max`:
*g(L) ∝ exp(−(L − mean)² / (2 sd²))*, renormalised to sum to 1.

**Placements.** A fragment of length *L* can start at
*n(t, L) = len(t) − L + 1* positions of transcript *t* if *L ≤ len(t)*. If
*L > len(t)*, the fragment is the whole transcript: *n(t, L) = 1*.

**Effective length.** *e(t) = Σ_L g(L) · n(t, L)*. For a transcript longer
than `--fragment-max` this is *len(t) − mean of g + 1*. It is the
`effective_length` column of the truth table.

**One fragment** is drawn as follows.

1. *Transcript.* With `--abundance-type tpm`: transcript *t* with probability
   *a(t) · e(t) / Σ_u a(u) · e(u)*. With `counts`: with probability
   *a(t) / Σ_u a(u)*. The number of fragments per transcript is therefore
   multinomial; its expectation is the `expected_fragments` column.
2. *Length and start.* Given *t*, length *L* with probability
   *g(L) · n(t, L) / e(t)*; then, if *L ≤ len(t)*, a start position uniform
   over the *n(t, L)* possible ones; otherwise the fragment is the whole
   transcript. (Implemented by rejection: draw *L* from *g* and a position
   *p* uniform in the transcript; accept if the fragment fits, or if
   *L > len(t)* and *p* is the first base.)
3. *Strand.* By `--library-type`. The read on the transcript strand starts at
   the 5′ end of the fragment; the other read is the reverse complement of
   the fragment and so starts at its 3′ end. The two reads always point
   towards each other ("inward").
4. *Read sequences.* Each read is the first `--read-length` bases of its
   strand of the fragment. If the fragment is shorter than the read, the read
   continues with the adapter of that read, and after the end of the adapter
   with the filler (`--adapter-fill`).
5. *Qualities.* For every cycle a quality is drawn from the profile (given
   the quality of the previous cycle if the profile records that).
6. *Errors.* Each base, including adapter and filler bases, is replaced with
   probability *10^(−Q/10)* by one of the three other bases, each with
   probability 1/3. With `--n-rate` > 0 a base is instead called `N` with
   quality 2 with that probability.

In `tpm` mode, steps 1 and 2 together mean that every combination of a
molecule, a length and a possible start has probability proportional to
*a(t) · g(L)*: transcripts are fragmented uniformly and fragments are
size-selected with *g*. Three consequences are worth knowing:

- The lengths of the simulated fragments follow *g(L) · Σ_t a(t) · n(t, L)*,
  not *g* itself: long fragments are slightly rarer because they fit at fewer
  positions. For human transcripts and the default parameters the mean is
  about 179 instead of 180.8 (the mean of the truncated normal). The report
  gives both numbers.
- Transcripts shorter than the drawn length yield their whole sequence, which
  may be shorter than `--fragment-min`. Such fragments are counted in the
  report (`whole_transcript_fragments`); with real human abundances they were
  about 0.001 % of the fragments (10 of 2 million).
- A transcript of abundance 0 never yields a fragment.

**Random numbers.** xoshiro256++ generators seeded through SplitMix64, both
implemented in the program so that results do not depend on the version of a
library. Discrete distributions are sampled with alias tables.

## What is not simulated

- **Sequence-specific and positional biases**: no random-hexamer priming
  bias at the start of reads, no GC bias, no 5′/3′ coverage bias, no
  degradation. Coverage is uniform along each transcript apart from the edge
  effect of the fragment length.
- **PCR duplicates**: none beyond the duplicates that arise by chance when
  many fragments are drawn from a short or highly expressed transcript.
- **Indels**: errors are substitutions only. There are no homopolymer errors.
- **Poly(A) tails**: transcripts are used as they are in the FASTA file. Add
  tails to the sequences if you need them.
- **Anything that is not in the FASTA file or has abundance 0**: no rRNA, no
  mitochondrial excess, no intronic or intergenic reads, no genomic DNA, no
  contamination, no adapter dimers, no PhiX, unless you provide sequences and
  abundances for them.
- **Sequence variation**: no SNPs, no allele-specific expression, no RNA
  editing.
- **Flow-cell structure**: one lane; qualities and errors do not depend on
  the tile or the position; no optical duplicates; no index hopping. The
  tile and coordinates in the read names are cosmetic.
- **Quality details**: the quality of a base depends on the cycle and, at
  most, on the previous quality. Qualities do not depend on the sequence, are
  not lower in the adapter or the poly-G filler, and the two reads of a pair
  are independent. A first-order chain gives fewer entirely bad reads than
  real data (see [Validation](#validation-with-real-tools)). Substitutions do
  not depend on the base or its neighbours.
- **`N` in the reference**: an ambiguous base of a transcript is written as
  `N` in the read, with a quality drawn like any other.
- **Single-end reads** and reads of different lengths are not implemented.
- **Biological variability between samples**: each run is one multinomial
  sample of the given abundances. To simulate replicates or conditions, give
  each sample its own abundance table and its own seed.

## Adapter sequences

The defaults are the Illumina TruSeq adapters with unique dual indexes, as
used by most stranded RNA-seq kits:

| Read | Runs into | Origin |
|---|---|---|
| 1 | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` + i7 index + `ATCTCGTATGCCGTCTTCTGCTTG` | The index (i7) adapter `(A)GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG`, including the `A` added by A-tailing. |
| 2 | `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` + i5 index + `GTGTAGATCTCGGTGGTCGCCGTATCATT` | The reverse complement of the i5 adapter `AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT`. |

Source: Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22
(pages "TruSeq Single Indexes" and "IDT for Illumina–TruSeq DNA and RNA UD
Indexes",
<https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/Illumina_DNA/AdapterSequencesIntro.htm>),
which also gives `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and
`AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` as the sequences to trim from read 1 and
read 2. The Cutadapt user guide lists the same two sequences and explains the
leading `A` (<https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq>).
Both were consulted in November 2026. Note that the full read-1 sequence has
one more base, `C`, between the 33-base trimming sequence and the index. The
oligonucleotide sequences are © Illumina, Inc.

The default indexes are those of the index pair UDI0001 of the same document:
i7 `CCGCGGTT`, i5 `AGCGCTAG`. Read 2 reads through the i5 index on the
opposite strand, so `--index2` is its reverse complement, `CTAGCGCT`; this is
also how instruments that read index 2 on the reverse strand (the "reverse
complement" workflow) write it in the FASTQ header. On forward-strand
instruments the header would show `AGCGCTAG` instead; the header is cosmetic.

With the defaults and a read length of 100, the adapters (66 and 70 bases)
are longer than what remains of the read after the shortest fragment
(40 bases), so the filler is only reached by whole-transcript fragments
shorter than 34 bases.

Other libraries:

```bash
# TruSeq single index (read 2 runs into the universal adapter, which has no index)
--index1 ATCACG --index2 '' \
--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# Nextera / Illumina DNA Prep style (same Illumina document)
--adapter1 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC{i7}ATCTCGTATGCCGTCTTCTGCTTG' \
--adapter2 'CTGTCTCTTATACACATCTGACGCTGCCGACGA{i5}GTGTAGATCTCGGTGGTCGCCGTATCATT'

# no adapters at all: short fragments are followed by the filler only
--adapter1 '' --adapter2 ''
```

## Reproducibility

The fragments are produced in chunks of 4,096. Chunk *i* has its own random
stream, derived from the seed and *i*, and the chunks are written in order.
Each chunk is compressed as one gzip member, so the compressed bytes do not
depend on the threads either. Consequently:

- the same program version, inputs, options and seed give byte-identical
  FASTQ files, truth tables and reports for any `--threads`;
- with a different `--compression-level` or with `--uncompressed`, the reads
  are the same and only the encoding of the files differs;
- the report does not contain the date, the run time or the number of
  threads; it does contain the file names as they were given.

Results are not guaranteed to be identical between program versions, and
have only been compared on one platform (x86-64 Linux).

## Performance

Measured with `scripts/benchmark.sh` on a desktop computer with 8 cores
(16 hardware threads) and an NVMe solid-state disk, on GENCODE v48 human
transcripts (412,044 sequences, 412 million bases) with abundances from a
real sample; 2 × 100 bases, default options, 4 million read pairs per run.
Each number is the best of three runs, in read pairs per second, and covers
simulating and writing the reads. Other programs were running on the computer
during the measurements, and single runs varied by up to about 30 %.

@@PERFORMANCE_TABLE@@

Reading the inputs (a 617 MB uncompressed FASTA file and a 412,044-row
table) takes a further 1 to 1.5 s per run. For comparison, the simulator of
RSEM was timed on the same computer at 1 million read pairs in 11 s on one
core, about 90,000 pairs per second, uncompressed.

Peak memory (maximum resident set size) with these inputs was 0.64 GB with
1 thread and 0.68 GB with 8 threads, and was the same for 200,000, 2 million
and 20 million read pairs. About 0.41 GB of it is the transcript sequences.

The generation of the reads takes about 1.1 to 1.6 µs per pair and thread
(measured without input and output on a small transcriptome,
`cargo test --release -- --ignored --nocapture throughput`); with gzip output
most of the time is spent compressing.

## Validation with real tools

A sample of 2 million read pairs was simulated from the human transcripts and
abundances described under [Performance](#performance) (2 × 100 bases,
default fragment lengths, built-in quality profile, seed 2026), once for each
library type, and analysed with released versions of common tools. The
simulator's report gave a mean fragment length of 178.8 bases and 7.14 % of
fragments shorter than the read.

### Salmon 2.7.0 (`salmon quant -l A`, untrimmed reads)

The index contained the same transcripts; transcripts with identical
sequences are collapsed by Salmon, so their truth was added up
(`scripts/compare_salmon.py`). Correlations are for log2(x + 1) of simulated
fragments against Salmon's `NumReads`, for transcripts or genes with at least
10 simulated fragments.

| Library type | Detected | Compatible fragments | Mapping rate | Transcripts: Pearson / Spearman | Genes: Pearson / Spearman |
|---|---|---|---|---|---|
| `reverse` | ISR | 99.989 % | 97.10 % | 0.914 / 0.958 (n = 19,364) | 0.995 / 0.998 (n = 11,233) |
| `forward` | ISF | 99.990 % | 97.11 % | 0.914 / 0.958 (n = 19,364) | 0.995 / 0.998 (n = 11,233) |
| `unstranded` | IU | 100 % | 97.11 % | 0.918 / 0.960 (n = 19,426) | 0.995 / 0.998 (n = 11,222) |

For input abundance (rescaled to TPM) against Salmon's TPM, the `reverse`
sample gave Pearson 0.947 and Spearman 0.958 at transcript level, and 0.987
and 0.984 at gene level. The median absolute relative difference between
simulated fragments and `NumReads` was 5.8 % for transcripts and 2.7 % for
genes. The lower agreement at transcript level reflects the difficulty of
assigning 2 million read pairs to the isoforms of a gene, not a property of
the reads: at gene level the agreement is close.

The 2.9 % of unmapped pairs are the fragments so short that most of the read
is adapter. After Trim Galore (below), Salmon mapped 99.99 % of the pairs of
the `reverse` sample (still ISR, 99.996 % compatible); the median relative
difference at gene level fell below 0.01 %, and the gene-level correlations
were 0.995 (Pearson) and 0.998 (Spearman).

### Trim Galore 2.3.0 (`trim_galore --paired`)

Trim Galore detected the adapter type by itself: "Using Illumina adapter for
trimming (count: 43081). Second best hit was smallRNA (count: 7)". Its
poly-G detection stayed off (8 of 1,000,000 reads with poly-G tails), as
expected since no read reached the filler.

It reported 36.6 % (read 1) and 38.5 % (read 2) of "reads with adapters".
This number cannot be compared directly with the 7.14 % of fragments shorter
than the read, because with the default minimum overlap of 1 base every read
that happens to end in `A` counts. The comparison is meaningful for longer
matches, using the fragment lengths from `--per-read-truth`:

| Bases removed | Trim Galore, read 1 | Trim Galore, read 2 | Simulated fragments that short |
|---|---|---|---|
| at least 5 | 6.17 % | 6.18 % | 6.13 % |
| at least 10 | 5.00 % | 4.99 % | 5.03 % |
| at least 20 | 3.21 % | 3.21 % | 3.24 % |

For each removed length from 10 to 60 bases, the number of reads differed
from the number of fragments of the corresponding length by at most 57 (of
about 400 to 4,000 reads per length).

### FastQC 0.12.1

Both files were processed without error or warning. Every module gave PASS
for read 1 and for read 2:

| Module | Result | Would a good real stranded library look like this? |
|---|---|---|
| Basic Statistics | PASS | Yes. |
| Per base sequence quality | PASS | Yes: median Q37 at every position, mean falling from 36.4 to 35.5 (read 1) and from 36.2 to 35.0 (read 2). Real runs usually show a small dip in the first cycles as well. |
| Per tile sequence quality | PASS | Yes, trivially: tiles do not differ. |
| Per sequence quality scores | PASS | Yes, but the distribution is narrower than in real data: almost all reads have a mean quality of 34 to 36, and there is no tail of bad reads. |
| Per base sequence content | PASS | **No.** Real RNA-seq libraries made with random hexamers show a biased composition in the first 10 to 12 bases, and FastQC warns or fails. This bias is not simulated. |
| Per sequence GC content | PASS | Yes; real libraries pass or warn. |
| Per base N content | PASS | Yes (no `N` calls by default). |
| Sequence Length Distribution | PASS | Yes for untrimmed reads. |
| Sequence Duplication Levels | PASS | At this depth, yes: 78 % of the reads would remain after deduplication. At 20 million read pairs of the same sample FastQC gave WARN (52 % remaining), as real RNA-seq libraries of that depth commonly do. In real data part of the duplication comes from PCR; here all of it comes from drawing many fragments from highly expressed transcripts. |
| Overrepresented sequences | PASS | Yes for a clean poly(A) library; real libraries sometimes show rRNA or adapter-dimer sequences. |
| Adapter Content | PASS | Yes: the Illumina Universal Adapter curve rises steadily towards the 3′ end and reaches 4.4 % (FastQC warns from 5 %), as expected with 5.3 % of fragments shorter than 89 bases. |

### Quality profile learned from real reads

`profile` was run on 200,000 real read pairs of 63 bases from an older
four-colour instrument (38 distinct quality values from 2 to 41), and 200,000
pairs were simulated with the result. The mean quality per cycle of the
simulated reads differed from the real ones by at most 0.03. The first-order
model reproduced only part of the tail of bad reads: 0.7 % of simulated
read 1 had a mean quality below 25, against 1.1 % of the real ones (1.1 %
against 1.6 % for read 2).

### Determinism

The 2-million-pair `reverse` sample was produced with 1, 3 and 8 threads:
the MD5 sums of both FASTQ files and of the truth table were identical, and
so were the reports. The automated tests check the same on a small example
for gzip and plain output.

### Not checked

- No genome aligner (STAR, HISAT2) was run on the reads.
- No tool other than those above was tried.
- Only human transcripts and read lengths of 63 and 100 bases were used in
  these runs.
- Platforms other than x86-64 Linux.

## Exit status

| Status | Meaning |
|---|---|
| 0 | Success. |
| 1 | The run failed: a file could not be read or written, or an input or a parameter value is invalid. A one-line message starting with `error:` says why. |
| 2 | The command line could not be parsed (unknown option, missing required option, value of the wrong kind). |

## Development

```bash
cargo test                       # unit tests and end-to-end tests
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

The unit tests cover the random number generators (against published
reference values), alias sampling, FASTA and table parsing, the fragment
model, strand logic, adapter read-through, the relation between quality and
error rate, read names, and the simulator (fragments per transcript against
the model, start positions, every read compared with its transcript for the
three library types, independence of the chunks). The end-to-end tests in
`tests/cli.rs` run the program on a tiny made-up transcriptome
(`tests/data`, see its README) and check the files, the error messages and
that the output is identical for different numbers of threads.

`scripts/compare_salmon.py` compares a Salmon `quant.sf` with a truth table
and `scripts/benchmark.sh` measures the speed; both were used for the numbers
above.

The source is organised as a small library (`src/lib.rs`) and a command-line
front end (`src/main.rs`): `rng`, `alias`, `fragment`, `fasta`, `abundance`,
`reads` (strand, adapters), `quality` (profiles, errors), `names`, and
`simulate` (the simulator and its multi-threaded driver).

## Licence

MIT; see [LICENSE](LICENSE). The adapter sequences are Illumina's and are
reproduced for the purpose of simulating and trimming reads of Illumina
libraries.

## References

- Illumina, Inc. *Illumina Adapter Sequences*, document 1000000002694 v22.
- M. Martin. Cutadapt removes adapter sequences from high-throughput
  sequencing reads. *EMBnet.journal* 17 (2011): 10–12. doi:10.14806/ej.17.1.200
- D. Blackman, S. Vigna. Scrambled linear pseudorandom number generators.
  *ACM Trans. Math. Softw.* 47 (2021): 36. doi:10.1145/3460772
- G. L. Steele, D. Lea, C. H. Flood. Fast splittable pseudorandom number
  generators. *OOPSLA* 2014. doi:10.1145/2660193.2660195
- M. D. Vose. A linear algorithm for generating random numbers with a given
  distribution. *IEEE Trans. Softw. Eng.* 17 (1991): 972–975.
  doi:10.1109/32.92917
- B. Li, C. N. Dewey. RSEM: accurate transcript quantification from RNA-Seq
  data with or without a reference genome. *BMC Bioinformatics* 12 (2011):
  323. doi:10.1186/1471-2105-12-323 (effective length and the generative
  model of fragments)
