# Documentation

📦 **Repository:** <https://github.com/sawers-rellan-labs/zealtiger>

Step-by-step record of the pipeline, the decisions made, and the gotchas hit.
Read in order; each file maps to one stage of the work.

**Rendered notebook:** [`nil_introgression_1400.html`](nil_introgression_1400.html) —
self-contained HTML of Stages 1–4 (map QC → 1400-NIL simulation → cM→Mb →
size-distribution figures → interference comparison → crossovers) plus the
RTIGER r=3 benchmark results. Regenerate with
`quarto render nil_introgression_1400.qmd` (the source `.qmd` is at the repo root).

**Rendered notebook:** [`molbreeding_vs_skim.html`](molbreeding_vs_skim.html) — concordance
of the MolBreeding target sequencing calls (source #5, **real ~110× allele counts**) vs the SNP50K skim
RTIGER calls (source #2) for the NILs on both platforms (donor dosage + B73-bulk control,
confusion matrix, segment scale, introgression-size distribution, donor composition).
**Frozen** direct comparison at a borrowed r=8. Regenerate with
`quarto render molbreeding_vs_skim.qmd`.

**Simulation-calibrated suite** (per-platform truth-derived rigidity; methodology
[`15-rigidity-and-scoring.md`](15-rigidity-and-scoring.md)):

- [`skim_sim_calibration.html`](skim_sim_calibration.html) — SNP50K skim (r=5):
  coverage-model calibration (layer 2), recovery vs truth (layer 3), three-`f_i`
  bias decomposition. `quarto render skim_sim_calibration.qmd`.
- [`molb_sim_calibration.html`](molb_sim_calibration.html) — MolBreeding target sequencing (r=2):
  same three layers. `quarto render molb_sim_calibration.qmd`.
- [`molbreeding_vs_skim_calibrated.html`](molbreeding_vs_skim_calibrated.html) —
  calibrated target sequencing-vs-skim comparison on the 14 shared NILs (calls at skim r=5 /
  target sequencing r=2; null = n=1500 simulation). `quarto render molbreeding_vs_skim_calibrated.qmd`.

| Doc | Step |
|---|---|
| [00-original-spec.md](00-original-spec.md) | The handover specification (input; provenance) |
| [00-skeleton-notebook.qmd.txt](00-skeleton-notebook.qmd.txt) | The handover skeleton notebook (input) |
| [01-map-compilation.md](01-map-compilation.md) | Stage 1 — compile the B73 v5 genetic map |
| [02-pedigree-and-simulation.md](02-pedigree-and-simulation.md) | Stage 2 — BC2S3 pedigree + 1400-NIL simulation |
| [03-cm-to-mb.md](03-cm-to-mb.md) | Stage 3 — convert donor segments cM → Mb |
| [04-distribution-and-interference.md](04-distribution-and-interference.md) | Stage 4 — size distribution + interference sensitivity |
| [05-crossovers.md](05-crossovers.md) | Crossovers per genome (for RTIGER autotune) |
| [06-rtiger-benchmark-dataset.md](06-rtiger-benchmark-dataset.md) | Generate the RTIGER benchmark with known truth |
| [07-rtiger-fitting.md](07-rtiger-fitting.md) | Fit RTIGER (Julia backend) on the benchmark |
| [08-scoring.md](08-scoring.md) | Score RTIGER calls vs ground truth |
| [09-rigidity-selection.md](09-rigidity-selection.md) | Choosing the rigidity parameter from first principles |
| [10-rtiger-fitting-performance.md](10-rtiger-fitting-performance.md) | Why fitting is slow (quadratic in rigidity), autotune cost, scaling |
| [11-benchmarking-other-tools.md](11-benchmarking-other-tools.md) | Benchmark *your own* ancestry caller: scoring schema + input-adapter gotchas |
| [12-het-excess-diagnosis.md](12-het-excess-diagnosis.md) | Diagnosing a heterozygote excess in real ancestry calls (selection vs heterosis vs paralogs vs caller) |
| [13-single-locus-validation.md](13-single-locus-validation.md) | Validate the simulation vs the BC2S3 single-locus expectation (f_i mean/variance, Hotelling, three-f_i bias decomposition) |
| [14-interference-and-rigidity.md](14-interference-and-rigidity.md) | Where m=10 comes from, why the map can't fit it, and how m sets the truth-optimal rigidity r |
| [15-rigidity-and-scoring.md](15-rigidity-and-scoring.md) | Deriving per-platform rigidity r from a truth sweep (vs RTIGER's `optimize_R`); FP/FN definitions side by side; the skim/target sequencing sweep tables |

## Pipeline at a glance

```
data/maizegdb_consensus_map/   (input)
        │  Stage 1  R/01  ─────────────────────────────►  results/maize_map_v5_clean.rds
        │                                                  results/marey_maps.pdf
        ▼  Stage 2  R/02  (simcross, 1400 NILs, no selection)
   donor segments (cM, union of mat/pat)
        │  Stage 3  R/03  (monotone cM→bp per chromosome)
        ▼
   donor segments (Mb)
        │  Stage 4  R/04  ─────────────────────────────►  results/nil_segments.csv
                                                           results/figures/*.png
                                                           results/run_metadata.json

   Benchmark track:
   R/05 + make_rtiger_benchmark.R  ──►  results/rtiger_benchmark/{counts,truth_*,expDesign}
   fit_rtiger.R (RTIGER + Julia)   ──►  results/rtiger_benchmark/rtiger_out/
   R/06 + score_rtiger.R           ──►  results/rtiger_benchmark/scoring/
```

## Key facts to carry forward

- **Map is a composite consensus map**: total **1,783 cM** (above the ~1,400–1,600
  single-population range). Mb segment sizes scale inversely with this — see
  [01](01-map-compilation.md).
- **Validation anchor**: genome-wide donor *dosage* converges to **0.125**
  (allele frequency), distinct from the *union* introgression fraction ≈0.14.
  See [02](02-pedigree-and-simulation.md).
- **SNP50K coverage model** (Fausto's fit): λ=0.43×, missingness
  `π + (1−π)e^(−kλ)` with π=0.161, k=1.042. See [06](06-rtiger-benchmark-dataset.md).
- **RTIGER rigidity is in *informative* markers**; recommended **r ≈ 8–10**.
  See [09](09-rigidity-selection.md).
- **Reproducibility**: global seed `20260609`, recorded in `run_metadata.json`.
