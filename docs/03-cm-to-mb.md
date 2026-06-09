# Stage 3 — Convert donor segments cM → Mb

**Code:** `R/03_segments_to_mb.R`, driven by `run_all.R --sim`.
**Output:** `results/nil_segments.csv` / `.rds` (one row per donor segment with
`start_cm, end_cm, start_bp, end_bp, mb`).

## What the code does

Per chromosome, build a monotone **cM → bp** interpolator from the **full clean
anchor set** (Stage 1, not the thinned grid):

- `make_cm_to_bp` collapses cM ties to mean bp (so the relation is a function),
  then `approxfun(cm, bp, rule = 2)` — linear interpolation, clamped at the ends.
- Segment Mb = `(bp(end_cm) − bp(start_cm)) / 1e6`.
- `assert_within_chr_length` guarantees no segment exceeds its chromosome's
  physical span.

Linear interpolation is the defensible baseline. A monotone spline
(`splinefun(method = "monoH.FC")`) is an optional smoother variant; it matters
mainly in recombination-cold pericentromeres, where small cM intervals map to
large Mb spans — the dominant driver of the largest physical introgressions.

## Why the full anchor set, not the grid

The 2,296-marker thinned grid is sparser than the 19,486 clean anchors. Using all
anchors makes the interpolation stable, especially across pericentromeres. The
grid is reserved for the optional marker-resolution / `convert2geno` sensitivity
check.

## Coordinate system

Everything is **B73 NAM v5 bp**, the same assembly as the MaizeGDB v5 columns and
as the real BZea genotype data (Schnable 2023 variants). So simulated and
observed segments are directly comparable with no liftover.
