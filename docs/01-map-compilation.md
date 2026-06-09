# Stage 1 — Compile the B73 NAM v5 genetic map

**Code:** `R/01_compile_map.R`, driven by `run_all.R` (and `nil_introgression.qmd` §1).
**Input:** `data/maizegdb_consensus_map/consensus_chr01.txt … chr10.txt`
**Output:** `results/maize_map_v5_clean.rds`, `results/maize_map_v5_grid2500.rds`,
`results/marey_maps.pdf`

Stage 1 is the QC gate: every Mb number downstream depends on the cM↔Mb map.

## What the code does

1. **Parse** (`parse_consensus_v5`) — read each consensus file (skip the free-text
   title line), `str_squish` the header, and locate the **B73 NAM v5** columns
   *by name*. This is mandatory because the v5 block sits at different column
   indices per file (chr01: cols 16–19; chr02: cols 12–15) and the chromosome
   labels are lowercase (`chr1`). Rows lacking v5 start/end are dropped; physical
   anchor `bp = (start + end) / 2`.
2. **Enforce monotonicity** (`enforce_monotone`) — within each chromosome, sort by
   bp and keep only markers whose cM is non-decreasing (running-max filter). This
   removes mismapped/paralogous loci that would invert the Marey map.
3. **Thin** (`thin_markers`) — allocate a ~2,500-marker quota per chromosome in
   proportion to cM span, picking markers evenly spaced in cM (a cap, not a
   guarantee).
4. **Marey QC plot** — cM vs Mb per chromosome; eyeballed before proceeding.

## Results

| | value |
|---|---|
| Files parsed | 10/10 |
| Raw v5 anchors | 20,051 |
| Clean monotone anchors | 19,486 (dropped 565, **2.8%**; per-chr 0.5–6.8%) |
| Thinned grid | 2,296 markers |
| **Total map length** | **1,783 cM** |

The Marey plots are textbook: steep recombination at telomeres, long
recombination-cold pericentromeres, strictly monotone, no inversions. Mb spans
match known v5 chromosome sizes (chr1 ≈ 308 Mb → chr10 ≈ 152 Mb).

## Decision / caveat

The total **1,783 cM exceeds the expected single-population range (~1,400–1,600 cM)**.
This is the signature of a MaizeGDB **composite consensus map** (it pools
recombination across many populations, inflating cM). The map is internally
consistent and monotone — not corrupted — but **donor segment sizes in Mb scale
inversely with map length**: on this map they run smaller than they would on a
~1,450 cM single-population map. Recorded in `run_metadata.json` (`total_cM`).
Decision (with Fausto): **proceed as-is**, report it as a map property.

## Note for the simulation

The simulation does **not** use the marker grid as a map. It uses
`chr_cm_lengths` = max cM per chromosome (the `L` vector for `simcross`). The
full clean anchor set is used for cM→Mb interpolation (Stage 3); the 2,500 grid
is only for an optional `convert2geno` sensitivity check.
