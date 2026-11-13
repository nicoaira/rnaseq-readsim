#!/usr/bin/env bash
# Measure the speed of `rnaseq-readsim simulate`.
#
# Usage: scripts/benchmark.sh TRANSCRIPTS.fa ABUNDANCE.tsv OUTPUT_DIR [PAIRS] [REPEATS]
#
# Runs 1, 4 and 8 threads (or the list in READSIM_THREADS), uncompressed and
# gzip (level 4 and level 1; or the list in READSIM_MODES: uncompressed,
# gzip-4, gzip-1), and prints for each configuration the best of REPEATS runs,
# all runs, the input loading time and the peak memory, in read pairs per
# second and kB. The time is the one the program reports for simulating and
# writing the reads; reading the inputs is reported separately. Extra options
# (for the abundance table, or library biases) can be given in the
# environment variable READSIM_OPTIONS.
set -euo pipefail

transcripts=$1
abundance=$2
outdir=$3
pairs=${4:-2000000}
repeats=${5:-3}
bin=${READSIM_BIN:-target/release/rnaseq-readsim}

mkdir -p "$outdir"
printf 'threads\toutput\tpairs\tbest_pairs_per_s\tall_pairs_per_s\tload_s\tmax_rss_kb\n'
for threads in ${READSIM_THREADS:-1 4 8}; do
  for mode in "uncompressed:--uncompressed" "gzip-4:--compression-level 4" "gzip-1:--compression-level 1"; do
    label=${mode%%:*}
    options=${mode#*:}
    case " ${READSIM_MODES:-uncompressed gzip-4 gzip-1} " in *" $label "*) ;; *) continue ;; esac
    rates=()
    load=""
    rss=0
    for _ in $(seq "$repeats"); do
      # shellcheck disable=SC2086
      out=$(/usr/bin/time -f 'max_rss_kb %M' "$bin" simulate -t "$transcripts" -a "$abundance" -o "$outdir/bench" \
        -n "$pairs" -p "$threads" $options ${READSIM_OPTIONS:-} 2>&1)
      line=$(grep 'read pairs in' <<<"$out")
      r=$(sed -nE 's/^max_rss_kb ([0-9]+)$/\1/p' <<<"$out")
      [ "$r" -gt "$rss" ] && rss=$r
      rates+=("$(sed -E 's/.*\(([0-9]+) pairs\/s\).*/\1/' <<<"$line")")
      load=$(sed -E 's/.*took ([0-9.]+) s.*/\1/' <<<"$line")
      rm -f "$outdir"/bench_[12].fastq "$outdir"/bench_[12].fastq.gz \
        "$outdir"/bench.truth.tsv "$outdir"/bench.report.json
    done
    best=$(printf '%s\n' "${rates[@]}" | sort -n | tail -1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$threads" "$label" "$pairs" "$best" "${rates[*]}" "$load" "$rss"
  done
done
