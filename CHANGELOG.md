# Changelog

## 0.2.0

- `simulate`: optional library biases, all off by default (the output of
  earlier options and seeds is unchanged): `--gc-bias-center` and
  `--gc-bias-width` (fragment GC preference), `--degradation` (3' bias of
  degraded RNA after poly(A) selection), `--hexamer-bias` (random-hexamer
  priming bias at the start of the reads) and `--library-molecules` (limited
  library complexity: PCR duplicates of the same fragment). Biases are applied
  by rejecting candidate fragments, so they change the number of fragments per
  transcript as in a real library. The report lists the new parameters,
  `candidate_fragments`, `fraction_of_candidates_kept` and
  `expected_unique_fraction`.

## 0.1.0

First version.

- `simulate`: paired-end reads from a transcript FASTA and an abundance table
  (TPM-like values or counts); reverse-stranded, forward-stranded or
  unstranded libraries; truncated-normal fragment lengths; adapter
  read-through with configurable adapters, indexes and filler; per-cycle
  quality profiles with substitution errors and optional `N` calls;
  instrument-style read names; gzip or plain output; per-transcript truth
  table, JSON run report and optional per-read truth table.
- `profile`: learns a quality profile (per-cycle distributions, optionally
  with first-order dependence on the previous cycle) from FASTQ files.
- Output is byte-identical for a given seed, whatever the number of threads.
