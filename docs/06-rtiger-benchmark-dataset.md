# RTIGER benchmark dataset (ground truth known)

**Code:** `R/05_make_rtiger_input.R`, driven by `make_rtiger_benchmark.R`.
**Output:** `results/rtiger_benchmark/` — `counts/sim_NNNN.tsv`, `expDesign.csv`,
`seqlengths.csv`, `truth_segments.csv`, `truth_markers.rds`, `sample_stats.csv`,
`params.json`, `README.md`.

Generate **100 simulated BC2S3 NILs genotyped at ~50K B73 v5 markers** under the
real **SNP50K** regime, recording the true introgression structure so RTIGER can
be scored (precision/recall/FDR, CO recovery, boundary error).

## Markers and orientation

- ~50K markers spread across the v5 chromosomes (∝ physical length, ~42.6 kb
  spacing). The real markers are teosinte-vs-B73 polymorphic sites (Schnable
  2023, 27.6M SNPs, B73 NAM v5), so all are donor/recurrent-informative.
- **Donor = teosinte = ALT; recurrent = B73 = REF.** Each marker's diploid
  genotype comes from the simulated mosaic: dosage 0/1/2 = homozygous-REF /
  het / homozygous-ALT.

## Coverage + missingness model (Fausto's SNP50K fit)

Missingness is **coverage-driven**, not independent noise. From Fausto's
exp-floor fit: `missing(λ) = π + (1−π)·e^(−kλ)` with, for SNP50K,
**λ = 0.43×, π = 0.161, k = 1.042 → ~0.70 missing**.

`draw_allele_counts` matches **both** moments exactly by construction:
`present_prob = (1−π)(1−e^(−kλ))`; present-site depth = `1 + Poisson(λ/present_prob − 1)`
(≈1.4 reads at λ=0.43, so a het site usually shows a single read and *looks
homozygous* — the regime that stresses the HMM). ALT count ~ Binomial(depth,
p_alt) with p_alt = {0, 0.5, 1} for dosage {0,1,2} and sequencing error 0.005.
Per-sample λ ~ Gamma (mean pinned to 0.43), **floored at 0.15** (the raw lower
tail produced λ≈0.01 / 99%-missing samples that real genotyping would exclude and
that break RTIGER's `2×rigidity` requirement — see [07](07-rtiger-fitting.md)).

Realized over 100 NILs: mean coverage **0.436**, mean missingness **0.712**
(per-sample 0.31–0.88), donor dosage 0.124.

## File format (matches RTIGER + the real `allelic_counts50K.tsv`)

6 columns, tab-separated, **no header**:
`Chr  Position  RefBase  RefCount  AltBase  AltCount`. Missing/uncovered sites
are kept as rows with `0  …  0` counts.

## Ground truth

- `truth_segments.csv` — runs of constant dosage state per sample
  (`name, chr, start_bp, end_bp, state(0/1/2), n_markers`).
- `truth_markers.rds` — per-marker true dosage for marker-level scoring.
- `sample_stats.csv` — per-sample λ, realized coverage/missingness, detectable
  COs (~21.5/line at this resolution vs 33.3 in the continuous mosaic — the gap
  is resolution + missingness loss), donor fraction.

To sweep difficulty, edit `lambda_mean` / `pi_floor` / `k_decay` and rerun.
