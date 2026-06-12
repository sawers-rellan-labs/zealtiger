# RTIGER coverage: estimator/cap benchmarking & run configuration

**Status.** The run config — **MLE, r=8, cap to 20× (`fit_rtiger_molbreeding.R` defaults)** —
is validated against **simulated ground truth** at both operating coverages, skim 0.4× and
seqcapture 20× ([§3b](#3b-validation-against-simulated-truth-at-the-two-operating-coverages)).
Capping at ≤20× makes the MLE-vs-MoM estimator choice moot, so there is no depth-based routing;
a routing alternative was explored against a non-truth reference and dropped ([§4](#4-run-configuration-defaults)).

## 1. Problem

RTIGER (faustovrz fork) fits a 3-state BetaBinomial HMM by EM. Its per-iteration cost
is driven by the **number of distinct `(alt_count, total_count)` pairs**, not by marker
count — both the E-step (`getlogpsi`, memoized BetaBinomial `logpdf` over distinct pairs)
and the M-step (`emissionUpdateState`, a numerical `optimize()` of the dispersion τ that
sums over distinct pairs, per state, per iteration) assume few distinct pairs. That holds
for low-coverage data and collapses for high coverage.

Distinct pairs scale with coverage, so **compute scales with coverage**:

| input | distinct (k,n) pairs | RTIGER fit (16 samples, 9,157 markers, r=8) |
|---|---|---|
| 0/1/2 pseudo-counts (fixed depth) | 3 | 13 s |
| GBTS real counts, downsampled 3× | 66 | 4 s |
| GBTS real counts, downsampled 30× | 817 | 18 s |
| GBTS real counts, **110×** | 7,555 | **240 s** |

This is why the old array-style 0/1/2 pseudo-count fit ran in seconds while real ~110×
GBTS counts take minutes: same markers, same samples, ~57× the distinct pairs.

## 2. Three approaches

- **MLE (stock).** Numerical `optimize()` of τ. Robust at small `n` (uses the exact
  discrete pmf), but its cost is the bottleneck at high coverage.
- **MoM (`mom_mstep.jl`).** Closed-form method-of-moments for τ (Williams/Pearson
  intra-class-correlation estimator from the gamma-weighted variance). Same BetaBinomial
  model, drop-in override of `emissionUpdateState`; removes the `optimize()` loop.
- **Cap (binomial thinning).** Down-thin counts to a target depth before fitting:
  `new_ref ~ Binomial(ref, p)`, `p = min(1, cap/total)`. Collapses distinct pairs →
  speeds *both* steps, and routes high-coverage data into the regime RTIGER already
  handles fast and well.

### Coverage as a Phred confidence (why a cap of ~10–20× is defensible)

The hard call in 3-state ancestry is het-vs-homozygote. A true het sequenced at depth `n`
mis-reads as homozygous (all reads one allele) with probability `2·0.5ⁿ`, i.e. a Phred
confidence `Q ≈ 3.01·(n−1)` — **each read adds ~3 Phred points**:

| depth | P(het reads as hom) | Q |
|---:|---|---:|
| 3× | 0.25 | 6 |
| 5× | 0.06 | 12 |
| 10× | 1/512 | 27 |
| 20× | 1.9×10⁻⁶ | 57 |
| 30× | 1.9×10⁻⁹ | 87 |

10× (Q≈27) is already ~variant-grade; 20× (Q≈57) makes each marker individually decisive
(matters at breakpoints / short introgressions where HMM smoothing can't help). Beyond
~20–30× the per-marker call adds nothing for Mb-scale ancestry.

## 3. Empirical sweep (MolBreeding GBTS, provisional)

Binomial-thinned the wideseq-filtered GBTS table (110×, 9,157 markers, 16 samples) to a
depth grid, fit r=8 with MLE and MoM, scored **per-marker state concordance vs the 110× MLE
fit**. ⚠️ The reference is the 110× MLE *fit*, **not ground truth** — see the caveat in §4.

| depth | est | pairs | fit time | dosage | conc. vs 110× MLE |
|---:|---|---:|---:|---:|---:|
| 110× | MLE | 7555 | 240 s | 0.069 | — (reference) |
| 110× | MoM | 7555 | 2.0 s | 0.092 | **0.924** |
| 30× | MLE | 817 | 18.2 s | 0.069 | 0.9987 |
| 30× | MoM | 817 | 0.9 s | 0.071 | 0.9912 |
| 20× | MLE | 504 | 14.5 s | 0.068 | 0.9984 |
| 20× | MoM | 504 | 0.4 s | 0.079 | 0.9930 |
| 10× | MLE | 222 | 6.6 s | 0.068 | 0.9948 |
| 10× | MoM | 222 | 0.2 s | 0.068 | 0.9911 |
| 5× | MLE | 113 | 6.3 s | 0.066 | 0.9921 |
| 3× | MLE | 66 | 4.2 s | 0.065 | 0.9890 |
| 3× | MoM | 66 | 0.2 s | 0.063 | 0.9872 |

MoM-vs-MLE agreement at each depth: 3× 0.997, 10× 0.994, 30× 0.991, **110× 0.924**.

**Observations (against this reference):**
1. **Capping is near-lossless.** Even 3× MLE is 98.9% concordant with 110× MLE; 20–30× MLE
   is 99.85%+ at 13–16× the speed. For 3-state ancestry the 110× depth is largely wasted.
2. **MoM matches MLE on capped data** (≤30×: ~99% vs reference, sub-second).
3. **MoM diverges from MLE on raw 110×** (92.4%, dosage 0.092 vs 0.069 — over-calls donor).
   Likely mechanism: real counts aren't perfectly BetaBinomial (overdispersion/outliers);
   that misspecification is invisible under low-coverage binomial noise but visible at 110×,
   where MoM's variance-matching over-weights outlier markers (worst at the clamped boundary
   states μ=0.01/0.99) while MLE is more robust. **But see §4 — we cannot currently tell
   whether this divergence is MoM erring or MoM correcting MLE.**

## 3b. Validation against simulated TRUTH at the two operating coverages

The §3 sweep scored vs the 110× MLE *fit*. The decision to trust RTIGER on the two
datasets actually used — SNP50K skim (~0.4×) and seqcapture capped to ~20× — is justified
against **simulated ground truth** (BC2S3 mosaics, simcross), one consistent metric, MLE r=8.
Skim: `results/rtiger_benchmark` (λ=0.43, realized 0.42×, ~70% missing, 50K markers, scored in
`scoring_auto_r8`). Seqcapture: `results/rtiger_benchmark_20x` (λ=20, realized 20.0×, 0.3% missing,
20K markers; `make_rtiger_benchmark.R` with `LAMBDA_MEAN=20 PI_FLOOR=0.002`).

| RTIGER r=8 MLE vs truth | skim 0.42× (70% miss) | seqcapture 20× (0.3% miss) |
|---|---|---|
| per-marker accuracy | 0.9986 | 0.99990 |
| recall: recurrent / het / donor | 0.9997 / 0.974 / 0.996 | 1.000 / 0.999 / 0.9996 |
| segment F1 — homozygous donor | 0.939 | 0.978 |
| segment F1 — heterozygous | 0.830 | 0.953 |
| segment F1 — donor-present | 0.917 | 0.975 |
| breakpoint median error | ~85–127 kb | ~0 kb |
| crossover bias / rmse | −2.71 / 3.73 | −0.84 / 1.41 |
| fit time | (low-coverage, fast) | 26 s (100 samples) |

**Both regimes are trustworthy; 20× is near-perfect.** The het call (the hard discrimination)
improves from F1 0.83→0.95 and breakpoints from ~100 kb to marker-exact. Caveat: the grids differ
(skim 50K markers, 20× 20K), so this is not a pure coverage isolation — but 20× resolving
breakpoints to ~0 kb *despite* coarser markers confirms coverage is the binding constraint.
This (plus the §3 cap sweep and §2 Q-values) is the basis for the run decision in §4.

## 4. Run configuration (defaults)

`fit_rtiger_molbreeding.R` defaults to the validated config — **`--estimator=mle --cap=20`**:
MLE E-step, binomial-thin to ~20×. The cap no-ops when a set's mean depth is already below 20×
(e.g. the 0.4× skim is untouched); `--cap=0` disables it. `--estimator=mom` switches to the
closed-form MoM (faster, but trust only on capped data — it diverges on raw high coverage, §3);
`--r` sets rigidity (default 8). There is no depth-based auto-routing: capping at ≤20× makes the
MLE-vs-MoM choice moot, so the estimator is a simple default, not a function of depth.

This config is justified by §3b (truth-validated at both operating coverages), the §3 cap sweep
(≤20× is ~99.85% concordant with the deepest data), and the §2 Q-values (20× ⇒ Q≈57 per het marker).

### A note on a depth-routing alternative (not used)

An earlier exploration considered routing the estimator by mean depth (MLE ≤9×, MoM 10–30×,
cap+MLE >30×). Those thresholds were eyeballed from §3 (scored vs the 110× MLE *fit*, not truth),
and capping at ≤20× removes the need for them, so the routing was dropped. §5 sketches how it
*would* be calibrated objectively if ever revisited.

### The limitation that motivates §5

§3 scores against the **110× MLE fit, which is itself an estimate, not truth**. Concordance
to a reference fit is circular: it cannot distinguish "MoM is wrong" from "MoM is right and
the 110× MLE is off." The single dataset is also one operating point (BC2S3 NILs, donor
dosage ~0.07, one marker density, one taxon set). The routing defaults should not be trusted
beyond this point until benchmarked against ground truth.

## 5. Toward methodic calibration (planned)

Use the project's **simulation ground truth** (Source 1, simcross BC2S3 mosaics) instead of
a reference fit, with one consistent objective — the same philosophy as RTIGER's rigidity
autotune (objectively minimize erroneous segments), but grounded in known truth:

1. **Truth:** simulate BC2S3 mosaics → known per-marker ancestry (`truth_segments.csv`,
   `make_joint_benchmark.R`).
2. **Coverage grid:** generate allele counts at 1,2,3,5,8,10,15,20,30,50,110× per true mosaic.
3. **Realistic noise — the linchpin:** the simulated counts MUST reproduce the real count
   structure (calibrated overdispersion, the real depth distribution, contamination/outlier
   markers). A clean binomial simulation would make MoM look fine at 110× and yield the wrong
   boundaries — the §3 MoM divergence is *caused* by real overdispersion that a naive sim omits.
   Calibrate to the GBTS count distribution (and reuse the fitted λ/π/k coverage models).
4. **Fit** {MLE, MoM} × {no cap, cap 10/20/30} at each depth.
5. **Score vs truth** with one uniform metric (segment F1 / per-marker accuracy / breakpoint
   distance) + record cost (fit time, #pairs).
6. **Defaults emerge:** at each depth, best (estimator, cap) by the metric under a cost budget;
   the regime thresholds are where the winner changes — not hand-set.
7. **Robustness:** repeat across donor dosage (BC1/BC2/BC3), marker density, taxon divergence,
   sample count — at least spot-check, and document the calibration regime.
8. **Bonus:** test whether RTIGER's internal autotune-style criterion (no truth needed)
   tracks the truth-based error; if so, route at runtime by that criterion rather than a
   hardcoded depth threshold.

Until §5 is done, treat §4 as a hypothesis. Simulation truth will, in particular, settle
whether MoM's 110× divergence is error or correction — the one thing §3 cannot.

---

*Scripts:* `molbreeding_downsample_gatk_table.py` (binomial thinning), `mom_mstep.jl`
(MoM M-step override), `fit_rtiger_molbreeding.R` (`--estimator`/`--cap` autoselect),
`compare_coverage_sweep.py` (concordance). Raw sweep numbers: `agent/coverage_sweep_findings.md`.
