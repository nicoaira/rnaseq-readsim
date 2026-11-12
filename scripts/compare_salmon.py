#!/usr/bin/env python3
"""Compare a Salmon quantification of simulated reads with the simulator's truth table.

Usage:
    compare_salmon.py TRUTH.tsv QUANT.sf [--duplicates duplicate_clusters.tsv]
                      [--tx2gene TX2GENE.tsv] [--min-fragments 10]

TRUTH.tsv is <prefix>.truth.tsv written by `rnaseq-readsim simulate`; QUANT.sf
is Salmon's output. If the Salmon index was built without --keepDuplicates,
transcripts with identical sequences were collapsed; pass the index's
duplicate_clusters.tsv so that the truth of the dropped copies is added to the
transcript that was kept. TX2GENE.tsv has two tab-separated columns
(transcript, gene) and no header.

Correlations are computed on log2(x + 1) for the features (transcripts or
genes) with at least --min-fragments simulated fragments, and Spearman's rank
correlation for the same features. Only the Python standard library is used.
"""
import argparse
import csv
import math
import sys


def read_table(path, key, columns):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row[key].split("|")[0]: [float(row[c]) for c in columns] for row in reader}


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            result[order[k]] = mean_rank
        i = j + 1
    return result


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def report(label, truth, estimate, min_fragments, unit):
    """truth and estimate: dict feature -> (count, tpm)."""
    features = [f for f, (count, _) in truth.items() if count >= min_fragments]
    if len(features) < 3:
        print(f"{label}: fewer than 3 features with >= {min_fragments} fragments")
        return
    for column, name in ((0, "fragments vs NumReads"), (1, "abundance vs TPM")):
        x = [truth[f][column] for f in features]
        y = [estimate.get(f, (0.0, 0.0))[column] for f in features]
        lx = [math.log2(v + 1) for v in x]
        ly = [math.log2(v + 1) for v in y]
        print(
            f"{label}\t{name}\t{unit} with >= {min_fragments:g} simulated fragments: {len(features)}"
            f"\tPearson(log2(x+1)) = {pearson(lx, ly):.4f}\tSpearman = {spearman(x, y):.4f}"
        )
    x = [truth[f][0] for f in features]
    y = [estimate.get(f, (0.0, 0.0))[0] for f in features]
    rel = sorted(abs(b - a) / a for a, b in zip(x, y))
    print(
        f"{label}\tfragments vs NumReads\tmedian absolute relative difference = "
        f"{rel[len(rel) // 2]:.4f}"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("truth")
    parser.add_argument("quant")
    parser.add_argument("--duplicates")
    parser.add_argument("--tx2gene")
    parser.add_argument("--min-fragments", type=float, default=10)
    args = parser.parse_args()

    truth_rows = read_table(args.truth, "transcript_id", ["fragments", "abundance"])
    quant = read_table(args.quant, "Name", ["NumReads", "TPM"])

    # the simulator's abundances, rescaled to sum to one million like TPM
    total_abundance = sum(v[1] for v in truth_rows.values())
    truth = {t: [c, a * 1e6 / total_abundance] for t, (c, a) in truth_rows.items()}

    if args.duplicates:
        moved = 0
        with open(args.duplicates, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                kept, dropped = row["RetainedRef"], row["DuplicateRef"]
                if dropped in truth and kept in truth:
                    truth[kept][0] += truth[dropped][0]
                    truth[kept][1] += truth[dropped][1]
                    del truth[dropped]
                    moved += 1
        print(f"# truth of {moved} duplicate transcripts added to the retained copies")

    missing = [t for t in truth if t not in quant and truth[t][0] > 0]
    if missing:
        print(
            f"# {len(missing)} transcripts with simulated fragments are not in the Salmon output",
            file=sys.stderr,
        )
    simulated = sum(v[0] for v in truth.values())
    estimated = sum(v[0] for v in quant.values())
    print(f"# simulated fragments: {simulated:.0f}; Salmon NumReads total: {estimated:.0f}")

    report("transcript", truth, quant, args.min_fragments, "transcripts")

    if args.tx2gene:
        gene_of = {}
        with open(args.tx2gene) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    gene_of[fields[0]] = fields[1]

        def by_gene(table):
            genes = {}
            for transcript, (count, tpm) in table.items():
                entry = genes.setdefault(gene_of.get(transcript, transcript), [0.0, 0.0])
                entry[0] += count
                entry[1] += tpm
            return genes

        report("gene", by_gene(truth), by_gene(quant), args.min_fragments, "genes")


if __name__ == "__main__":
    main()
