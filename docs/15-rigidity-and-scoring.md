# Rigidity selection and segment scoring — methodology notes

**Date:** 2026-06-17
**Context:** Redo of `molbreeding_vs_skim`. How the per-platform RTIGER rigidity
`r` is **derived** (not autotuned), and the two FP/FN definitions in play — ours
(counted vs the simulated truth) and RTIGER's `optimize_R` (estimated under an
iid model). Companion to `13-single-locus-validation.md` (single-locus / three-`f_i`)
and `14-interference-and-rigidity.md` (m, the map, the truth-optimal r).

---

## 1. Why we derive r from a truth sweep — not from `optimize_R`

The rigidity `r` is the minimum run length (in **informative markers**) RTIGER
will keep as a segment. We set it **per platform** by minimizing the number of
**counted** erroneous segments (FP+FN) against the n=1500 BC2S3 simulation truth,
re-degraded through each platform's calibrated coverage model
(`docs/13`, `docs/14`). The rigidity is a **resolution** knob, not a fit
parameter; the sweep picks the resolution that recovers the simulated truth most
faithfully.

We do **not** use RTIGER's built-in `optimize_R`. As documented in
`docs/09-rigidity-selection.md` and `docs/14-interference-and-rigidity.md`,
`optimize_R` over-rigidifies: it estimates erroneous-segment *rates* from an
**iid per-marker LLR** model (§2), which understates how strongly correlated real
genotyping errors are, so the objective keeps climbing and recommends a far
larger r than the truth sweep does (on the skim sim it recommends `r=256` while
the truth-optimal is `r≈2–5` — see memory `optimize-R-overrigidifies-nil-sim`).

Per-platform is necessary because coverage and informative-marker density differ.
The truth-optimal r scales with physical resolution:

    true_r(platform) ≈ L_floor(m) × informative_marker_density(platform)

The skim (~6.7 informative markers/Mb at 0.4×, ~70% missing) and GBTS
(~4.3 informative markers/Mb on the wideseq-filtered panel, but ~110× and
near-complete) sit in different regimes, so a shared r (the borrowed `r=8` of the
frozen `molbreeding_vs_skim.qmd`) does **not** give a matched physical floor.
See `docs/14` §4.

---

## 2. FP and FN, side by side

Both methods talk about "false" and "missed" segments, but one **counts** them
against a known truth and the other **estimates their probability** under a model.

### Ours — `R/06_score_rtiger.R::score_segments` (counted, vs truth)

- An **event** = a maximal run of states ∈ `target` (e.g. donor-present = states
  {1,2}, het = {1}, hom-donor = {2}).
- A truth event is a **TP** if a called event **reciprocally overlaps** it
  ≥ `min_overlap` (default 0.5) in **both** directions — i.e.
  `ov_len/len_truth ≥ 0.5` **and** `ov_len/len_called ≥ 0.5`
  (`score_rtiger.R:152-159`).
- **FN** = an unmatched **truth** event (a real segment the caller **missed**):
  `FN = n_truth_events − TP`.
- **FP** = an unmatched **called** event (a **false** segment the caller invented):
  `FP = n_called_events − n_matched_called`.
- Precision = TP/(TP+FP), recall = TP/(TP+FN), F1 = harmonic mean.

This requires the simulated truth. The sweep error metric is `FP+FN` on the
donor-present class (states 1+2).

### RTIGER `optimize_R` — `~/Desktop/RTIGER/R/Optimization_Rvalue.R` (estimated, no truth)

- Objective = expected number of erroneous **segments**,
  `SE_total = SE+ + SE−`.
- `exactFPR` = P(a uniform length-`u` region is **falsely broken** into a spurious
  segment) — computed by summing a per-marker log-likelihood ratio
  `log P(marker | other state) − log P(marker | flanking state)` and asking how
  often the running sum crosses the rigidity / transition threshold.
- `exactFNR` = P(a true central segment of a given length is **missed** /
  absorbed).
- The per-marker LLR increments are convolved as **iid** — the assumption that is
  violated in practice (errors cluster: mismapping, paralogs, local coverage
  dropouts).
- Here `r` and segment length are in **markers**; `total_length` and the CO-rate
  prior are in bp / Mb.
- **No truth is used** — everything is a probability under the emission model.

### One-line distinction

Ours **counts** errors against the simulated truth; `optimize_R` **estimates**
their probability under an iid LLR model — and because real errors are correlated,
the iid estimate inflates the apparent cost of a low r, so `optimize_R`
over-rigidifies.

---

## 3. Rigidity sweeps (chosen r marked)

Columns: `r` | resolution median (Mb, window-span median) | resolution mean
(= r · average informative-marker spacing AM, Mb) | F1 for the **introgression**
class (donor-present, states 1+2) | F1 het (state 1) | F1 hom-donor (state 2).
F1s are introgression / het / hom-donor — **not** REF.

Resolution conventions (`docs/14`; scripts `agent/informative_resolution.R`,
`agent/gbts_resolution.R`): **mean = r·AM** (AM = average informative-marker
spacing) and **median = the window-span median** — *not* `r·median(gap)` (which
underestimates) and *not* the harmonic mean (degenerate).

### Skim sweep (SNP50K, ~0.4×; informative density 6.7/Mb)

Source: `results/sim_calibration/skim_rigidity_sweep.csv`.

| r | res median | res mean (r/6.7) | F1 intro | F1 het | F1 hom-donor |
|---|---|---|---|---|---|
| 2 | 0.07 | 0.29 | 0.945 | 0.900 | 0.968 |
| 3 | 0.20 | 0.44 | 0.944 | 0.896 | 0.970 |
| **5 (chosen)** | 0.45 | 0.73 | 0.946 | 0.889 | 0.966 |
| 8 | 0.83 | 1.17 | 0.944 | 0.870 | 0.957 |
| 10 | 1.09 | 1.47 | 0.936 | 0.861 | 0.950 |

Counted FP+FN over r = 2,3,5,8,10 → **140 / 142 / 137 / 141 / 160**: a flat
plateau from r=2 to r=8 with the **argmin at r=5**. CO bias drifts
−1.55 → −2.08 as r grows (more rigidity → more crossovers absorbed). We take
**r=5**: it is the error argmin and keeps the donor-present F1 at its peak (0.946)
while holding het F1 reasonable; pushing lower buys nothing on FP+FN and pushing
higher degrades het and CO recovery.

### GBTS sweep (MolBreeding target-capture, ~110× → capped 20×; wsfilt density 4.3/Mb)

Source: `results/sim_calibration/gbts_rigidity_sweep.csv`.

| r | res median | res mean (r/4.3) | F1 intro | F1 het | F1 hom-donor |
|---|---|---|---|---|---|
| **2 (chosen)** | 0.094 | 0.46 | 0.985 | 0.964 | 0.990 |
| 3 | 0.248 | 0.69 | 0.981 | 0.961 | 0.987 |
| 5 | 0.577 | 1.15 | 0.976 | 0.941 | 0.980 |
| 8 | 1.082 | 1.84 | 0.953 | 0.903 | 0.967 |
| 10 | 1.413 | 2.30 | 0.938 | 0.889 | 0.960 |

Counted FP+FN over r = 2,3,5,8,10 → **39 / 48 / 62 / 118 / 153**: **monotone**
increasing, **argmin at r=2**. CO bias drifts −0.48 → −2.06. We take **r=2**:
the error minimum and the F1 peak across all three classes. GBTS tolerates the
finer resolution because ~110× coverage resolves het and crossovers the skim
misses — so the binding constraint is the interference floor, not coverage, and
r=2 sits right at it.

### Why the two platforms land on different r

The skim's coarse informative density (6.7/Mb but ~70% missing, low confidence per
marker) makes its error curve **flat** down to r=2 — there is no resolution to be
gained below ~0.5 Mb, so r=5 is chosen at the plateau argmin. GBTS's near-complete
high-confidence calls make the error curve **monotone**, rewarding the finest r=2.
Same simulated biology, two coverage regimes, two derived resolutions — exactly
the per-platform tuning `docs/14` argues for.
