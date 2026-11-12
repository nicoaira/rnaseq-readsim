# Changelog

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
