#!/usr/bin/env python3
"""Compare the MolBreeding GBTS coverage sweep: {MLE, MoM} x {3,5,10,15,20,30x}
against the 110x MLE gold, on the wideseq marker grid.

For each call set: mean donor dosage, #donor blocks/sample, and per-marker state
concordance vs the 110x MLE calls. Also MoM-vs-MLE agreement at each depth, and
fit wall-time (from fit_progress.log) + distinct (k,n) pairs (from the downsample).

Usage: python3 compare_coverage_sweep.py
"""
import csv
from bisect import bisect_right
from pathlib import Path

GBTS = Path("data/molbreeding_gbts")
GOLD = GBTS / "SNP_wsfilt" / "calls_common_schema.csv"          # 110x MLE
GRID = Path("data/molbreeding_45k/wideseq_keep_v5.tsv")
DEPTHS = [3, 5, 10, 15, 20, 30]
PAIRS = {3: 66, 5: 113, 10: 222, 15: 353, 20: 504, 30: 817, 110: 7555}


def load_calls(path):
    """(name, chr) -> (sorted starts, ends, states)."""
    segs = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            k = (r["name"], int(r["chr"]))
            segs.setdefault(k, []).append(
                (int(r["start_bp"]), int(r["end_bp"]), int(r["state"])))
    idx = {}
    for k, v in segs.items():
        v.sort()
        idx[k] = ([s[0] for s in v], [s[1] for s in v], [s[2] for s in v])
    return idx


def state_at(idx, name, chrom, pos):
    k = (name, chrom)
    if k not in idx:
        return None
    starts, ends, states = idx[k]
    i = bisect_right(starts, pos) - 1
    if i >= 0 and starts[i] <= pos <= ends[i]:
        return states[i]
    return None


def load_grid():
    g = []
    with open(GRID) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            g.append((int(c[0].replace("chr", "")), int(c[1])))
    return g


def dosage_and_segs(path):
    num = den = 0.0
    nseg = 0
    donor_blocks = 0
    prev = None  # (name, chr) last donor-state to count blocks
    rows = sorted(csv.DictReader(open(path)),
                  key=lambda r: (r["name"], int(r["chr"]), int(r["start_bp"])))
    last_key = None
    last_donor = False
    for r in rows:
        l = int(r["end_bp"]) - int(r["start_bp"]) + 1
        st = int(r["state"])
        num += l * st / 2
        den += l
        nseg += 1
        key = (r["name"], r["chr"])
        donor = st > 0
        if key != last_key:
            last_donor = False
        if donor and not last_donor:
            donor_blocks += 1
        last_donor = donor
        last_key = key
    return num / den, nseg, donor_blocks


def concordance(idx_a, idx_b, grid, names):
    agree = total = 0
    for nm in names:
        for chrom, pos in grid:
            sa = state_at(idx_a, nm, chrom, pos)
            sb = state_at(idx_b, nm, chrom, pos)
            if sa is None or sb is None:
                continue
            total += 1
            if sa == sb:
                agree += 1
    return agree / total if total else float("nan")


def fit_time(path):
    p = Path(path)
    if not p.exists():
        return None
    last = None
    for line in p.read_text().splitlines():
        if "elapsed=" in line:
            last = float(line.split("elapsed=")[1].split()[0])
    return last


def main():
    grid = load_grid()
    gold = load_calls(GOLD)
    names = sorted({k[0] for k in gold})
    gd, gseg, gblk = dosage_and_segs(GOLD)
    gt = fit_time(GBTS / "SNP_wsfilt" / "fit" / "fit_progress.log")
    gt_str = f"{gt:.0f}s" if gt is not None else "NA"
    print(f"GOLD = 110x MLE: dosage {gd:.3f}, {gseg} segs, {gblk/len(names):.1f} blocks/samp, "
          f"{gt_str}, {PAIRS[110]} pairs, {len(names)} samples\n")

    hdr = f"{'depth':>5} {'est':>4} {'pairs':>6} {'time_s':>7} {'dosage':>7} {'segs':>5} {'blk/s':>6} {'conc_vs_110xMLE':>16}"
    print(hdr)
    print("-" * len(hdr))
    mom110 = GBTS / "SNP_wsfilt_mom" / "calls_common_schema.csv"
    if mom110.exists():
        d, s, b = dosage_and_segs(mom110)
        c = concordance(load_calls(mom110), gold, grid, names)
        t = fit_time(GBTS / "SNP_wsfilt_mom" / "fit" / "fit_progress.log")
        t_str = f"{t:.1f}" if t is not None else "NA"
        print(f"{'110':>5} {'MoM':>4} {PAIRS[110]:>6} {t_str:>7} {d:>7.3f} {s:>5} {b/len(names):>6.1f} {c:>16.4f}")

    rows_out = []
    for est, sub in (("MLE", ""), ("MoM", "_mom")):
        for X in DEPTHS:
            path = GBTS / f"d{X}{sub}" / "calls_common_schema.csv"
            if not path.exists():
                continue
            d, s, b = dosage_and_segs(path)
            c = concordance(load_calls(path), gold, grid, names)
            t = fit_time(GBTS / f"d{X}{sub}" / "fit" / "fit_progress.log")
            ts = f"{t:.1f}" if t is not None else "NA"
            print(f"{X:>5} {est:>4} {PAIRS[X]:>6} {ts:>7} {d:>7.3f} {s:>5} {b/len(names):>6.1f} {c:>16.4f}")
            rows_out.append((est, X, c))

    # MoM vs MLE agreement at each depth
    print("\nMoM-vs-MLE per-marker agreement at each depth:")
    for X in DEPTHS:
        a = GBTS / f"d{X}" / "calls_common_schema.csv"
        m = GBTS / f"d{X}_mom" / "calls_common_schema.csv"
        if a.exists() and m.exists():
            c = concordance(load_calls(a), load_calls(m), grid, names)
            print(f"  d{X}x: {c:.4f}")
    am = GBTS / "SNP_wsfilt" / "calls_common_schema.csv"
    mm = GBTS / "SNP_wsfilt_mom" / "calls_common_schema.csv"
    if mm.exists():
        print(f"  110x: {concordance(load_calls(am), load_calls(mm), grid, names):.4f}")


if __name__ == "__main__":
    main()
