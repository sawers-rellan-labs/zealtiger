#!/usr/bin/env python3
"""Binomial-thin the MolBreeding GBTS GATK table to target mean depths, for the
coverage sweep (call quality + compute vs coverage).

For each marker, keep each read independently with p = min(1, target/actual):
    new_ref ~ Binomial(REF_COUNT, p)
    new_alt ~ Binomial(ALT_COUNT, p)
This reproduces the real sampling noise you'd get at the target depth (a true het
at 5x can show 5/0), not just a rescaling. Seeded for reproducibility.
Uses stdlib random.binomialvariate (Python 3.12+); no numpy dependency.

Input : data/molbreeding_45k/gatk_table_SNP_wsfilt_v5.tsv  (~110x, wideseq-filtered)
Output: data/molbreeding_45k/gatk_table_d<X>_v5.tsv  for X in the depth grid
        (same GATK schema; feeds fit_rtiger_molbreeding.R as set "d<X>")
        + prints realized mean depth and #distinct (alt,total) pairs per depth
          (the latter is RTIGER's per-iteration cost driver).

Usage: python3 molbreeding_downsample_gatk_table.py
"""
import random
from pathlib import Path

SRC = Path("data/molbreeding_45k/gatk_table_SNP_wsfilt_v5.tsv")
OUTDIR = Path("data/molbreeding_45k")
DEPTHS = [3, 5, 10, 15, 20, 30]
SEED = 42

binom = random.binomialvariate


def load(src):
    meta, ref, alt = [], [], []
    with src.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            # SAMPLE CONTIG POSITION REF_NUC ALT_NUC REF_COUNT ALT_COUNT TOTAL_COUNT
            meta.append((c[0], c[1], c[2], c[3], c[4]))
            ref.append(int(c[5]))
            alt.append(int(c[6]))
    return header, meta, ref, alt


def main():
    random.seed(SEED)
    header, meta, ref, alt = load(SRC)
    tot = [r + a for r, a in zip(ref, alt)]
    n = len(meta)
    uncapped_pairs = len({(a, t) for a, t in zip(alt, tot)})
    print(f"source: {SRC.name}  ({n:,} rows)")
    print(f"  uncapped: mean depth {sum(tot)/n:.1f}x | distinct (alt,total) pairs {uncapped_pairs:,}")
    print(f"{'depth':>6} {'mean_x':>7} {'distinct_pairs':>15} {'covered%':>9}")

    for X in DEPTHS:
        out = OUTDIR / f"gatk_table_d{X}_v5.tsv"
        pairs = set()
        sum_t = 0
        covered = 0
        with out.open("w") as w:
            w.write("\t".join(header) + "\n")
            for (s, contig, pos, rb, ab), r, a, t in zip(meta, ref, alt, tot):
                p = 1.0 if t <= X else X / t
                r2 = binom(r, p) if r else 0
                a2 = binom(a, p) if a else 0
                t2 = r2 + a2
                sum_t += t2
                if t2 > 0:
                    covered += 1
                pairs.add((a2, t2))
                w.write(f"{s}\t{contig}\t{pos}\t{rb}\t{ab}\t{r2}\t{a2}\t{t2}\n")
        print(f"{X:>5}x {sum_t/n:>7.1f} {len(pairs):>15,} {100*covered/n:>8.1f}%")


if __name__ == "__main__":
    main()