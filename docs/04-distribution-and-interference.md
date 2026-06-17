# Stage 4 — Size distribution + interference sensitivity

**Repo:** [sawers-rellan-labs/zealtiger](https://github.com/sawers-rellan-labs/zealtiger)

**Code:** `R/04_summarize.R`, driven by `run_all.R --sim`; also
`nil_introgression_1400.qmd` §4.
**Output:** `results/nil_summary.*`, `results/figures/introgression_size.{png,pdf}`,
`results/figures/interference_comparison.{png,pdf}`,
`results/interference_comparison.csv`, `results/run_metadata.json`

## Per-NIL summary

`summarise_nils` collapses to one row per line (total donor Mb, segment count),
filled to the full 1,400 so zero-introgression lines are represented.

## Introgression-size distribution (primary, m = 10)

| quantity | value |
|---|---|
| Mean total donor genome / NIL | **296 Mb** (~14% of the ~2.1 Gb genome ≈ the union fraction) |
| Segments / line | ~10 |
| Median donor segment | **10.7 Mb** |
| 90th percentile | 94.5 Mb |
| Max segment | 244 Mb |

Heavily right-skewed: an exponential-like bulk plus a long pericentromeric tail
(small cM intervals in cold regions → huge Mb spans), exactly as predicted.

## Crossover interference (spec §8.2) — decision: report both

`simcross` uses the Stahl model; its default is **m = 10, p = 0** (strong
positive interference, biologically defensible for maize). We ran the default
**and** a no-interference pass (m = 0) to bound the tail.

| model | median (Mb) | p90 | p99 | max |
|---|---|---|---|---|
| default (m=10) | 10.7 | 94.5 | 174 | 244 |
| no interference (m=0) | 9.5 | 94.7 | 184 | 296 |

Interference preserves the *mean* CO count but clusters crossovers, so it
slightly lengthens the median and **trims the far tail** (max 244 vs 296 Mb).
Effect is modest. Both interference settings are recorded in `run_metadata.json`.

## Reproducibility

`run_metadata.json` records the seed, marker count, total cM, genome Mb,
`simcross` version, both interference runs' donor dosage and mean total Mb, and
the crossover statistics ([05](05-crossovers.md)).
