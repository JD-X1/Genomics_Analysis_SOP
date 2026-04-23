#!/usr/bin/env python3
"""
Compute assembly statistics for metagenome bins produced by MaxBin2, CONCOCT,
and MetaBAT2. Outputs a CSV with per-bin metrics.

Usage:
    python3 bin_stats.py <binning_output_dir> [output_csv]

Outputs to <binning_output_dir>/bin_stats_summary.csv if no output path given.
"""

import csv
import glob
import os
import sys
from statistics import median


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


# ---------------------------------------------------------------------------
# Per-bin statistics
# ---------------------------------------------------------------------------

def gc_content(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq) - seq.count("N")
    return round(100.0 * gc / total, 4) if total > 0 else 0.0


def nx_lx(lengths, fraction=0.5):
    """Return (Nx length, Lx count) given sorted-descending contig lengths."""
    target = sum(lengths) * fraction
    running = 0
    for i, l in enumerate(lengths, 1):
        running += l
        if running >= target:
            return l, i
    return 0, len(lengths)


def bin_stats(fasta_path):
    seqs = [(hdr, seq) for hdr, seq in parse_fasta(fasta_path) if seq]
    if not seqs:
        return None

    lengths = sorted([len(s) for _, s in seqs], reverse=True)
    total_len = sum(lengths)
    gc = round(sum(gc_content(s) * len(s) for _, s in seqs) / total_len, 4)

    n50, l50 = nx_lx(lengths, 0.5)
    n90, l90 = nx_lx(lengths, 0.9)

    return {
        "num_contigs":       len(lengths),
        "total_length":      total_len,
        "largest_contig":    lengths[0],
        "mean_length":       round(total_len / len(lengths), 2),
        "median_length":     int(median(lengths)),
        "n50":               n50,
        "l50":               l50,
        "n90":               n90,
        "l90":               l90,
        "gc_pct":            gc,
    }


# ---------------------------------------------------------------------------
# Bin discovery
# ---------------------------------------------------------------------------

BINNER_PATTERNS = {
    "maxbin2":  ("maxbin2",            "*.fasta"),
    "concoct":  ("concoct/fasta_bins", "*.fa"),
    "metabat2": ("metabat2",           "bin.*.fa"),
    "quickbin": ("quickbin",           "bin*.fa"),
}


def discover_bins(binning_dir):
    """Yield (sample, binner, bin_name, fasta_path) for all bins."""
    for sample_dir in sorted(glob.glob(os.path.join(binning_dir, "*"))):
        if not os.path.isdir(sample_dir):
            continue
        sample = os.path.basename(sample_dir)
        for binner, (subdir, pattern) in BINNER_PATTERNS.items():
            bin_dir = os.path.join(sample_dir, subdir)
            if not os.path.isdir(bin_dir):
                continue
            for fasta in sorted(glob.glob(os.path.join(bin_dir, pattern))):
                bin_name = os.path.splitext(os.path.basename(fasta))[0]
                yield sample, binner, bin_name, fasta


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

FIELDS = [
    "sample", "binner", "bin_name",
    "num_contigs", "total_length", "largest_contig",
    "mean_length", "median_length",
    "n50", "l50", "n90", "l90",
    "gc_pct",
]


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    binning_dir = sys.argv[1]
    out_csv = sys.argv[2] if len(sys.argv) > 2 else os.path.join(binning_dir, "bin_stats_summary.csv")

    rows = []
    for sample, binner, bin_name, fasta in discover_bins(binning_dir):
        stats = bin_stats(fasta)
        if stats is None:
            print(f"[WARN] empty FASTA, skipping: {fasta}", file=sys.stderr)
            continue
        rows.append({"sample": sample, "binner": binner, "bin_name": bin_name, **stats})

    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} bins → {out_csv}")

    # Brief summary to stdout
    from collections import defaultdict
    counts = defaultdict(lambda: defaultdict(int))
    for r in rows:
        counts[r["sample"]][r["binner"]] += 1
    print(f"\n{'Sample':<40} {'MaxBin2':>8} {'CONCOCT':>8} {'MetaBAT2':>9} {'QuickBin':>9}")
    print("-" * 78)
    for sample in sorted(counts):
        b = counts[sample]
        print(f"{sample:<40} {b.get('maxbin2',0):>8} {b.get('concoct',0):>8} {b.get('metabat2',0):>9} {b.get('quickbin',0):>9}")


if __name__ == "__main__":
    main()
