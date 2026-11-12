# Test data

A tiny, made-up transcriptome for the integration tests: eight transcripts of
random sequence (70 to 2,400 bases, 6,970 bases in total) with headers in the
style of GENCODE, and an abundance table. None of the sequences is from a real
genome; the identifiers (`TOY…`) and gene names are invented.

- `transcripts.fa`: the sequences.
- `abundance.tsv`: columns `transcript_id` and `tpm`; `TOY0008.1` has
  abundance 0.
