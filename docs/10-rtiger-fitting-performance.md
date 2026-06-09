# RTIGER fitting performance — why it's slow and how it scales

Notes on the runtime/memory cost of `fit_rtiger.R`, why `autotune=TRUE` is much
slower than a fixed rigidity, and what to expect at full (1,439-sample) scale.
Companion to [07-rtiger-fitting.md](07-rtiger-fitting.md) (mechanics) and
[09-rigidity-selection.md](09-rigidity-selection.md) (choosing `r`).

## What is being computed

RTIGER fits a 3-state HMM (AA / AB / BB) by **Baum-Welch (EM)**: repeatedly
sweep **forward-backward** over every marker of every sample, re-estimate
parameters, repeat until the likelihood converges (up to `max.iter`, default
50). Then **Viterbi**-decode each sample and export per-sample segments. Cost per
EM iteration is roughly:

```
(samples × markers) × states² × 2   (forward + backward)
```

The fit is **single-threaded** (Julia backend via JuliaCall, embedded in the R
process — so it shows as ~100% of *one* core under the R process, with no
separate `julia` process).

## Measured breakdown (100 samples, r=10) — EM is cheap; export dominates

Instrumenting the fit (per-EM-iteration `progress_log`, see
[07](07-rtiger-fitting.md)) and timing the phases gave the **real** picture,
which differs from the a-priori EM cost model in the next section:

- **EM (parameter fit): ~10 s.** Converges in **5 iterations** (`delta`
  30.8 → 0.003 < `eps=0.01`, ~2 s/iter). Not the bottleneck.
- **Viterbi decode:** folded into the Julia call; fast.
- **Per-sample export (only when `save.results=TRUE`): the dominant phase.**
  `export2IGV` + `plotGenotype` write **3 bigWig tracks + 1 genotype PDF per
  sample** (300 bigWigs + 100 PDFs for 100 samples), all **R-side** — ~3 min.
  These are *visualization artifacts, not the analysis*.

The earlier "~30 min" run was the **unoptimized `main` branch**; the
`faustovrz/RTIGER@optimize-julia-core` fork cut the EM/E-step cost, but export
still dominates when `save.results=TRUE`.

**Fast benchmark path — skip the export.** Run `fit_rtiger.R … 0` (third arg
`save_results=FALSE`) and score from the saved `.rds` (the scorer reads segments
from `@Viterbi`). The fit then collapses to EM + Viterbi, with **identical**
scores:

| | `save.results=TRUE` | `save.results=FALSE` |
|---|---|---|
| wall time (100 samples, r=10) | ~3 min (≈30 min on unoptimized `main`) | **29 s** |
| outputs | 400 BEDs + 300 bigWigs + 100 PDFs | 1 `.rds` + 1 progress log |
| segment scores (homo / donor / het, P·R) | 0.999/0.870, 0.979/0.854, 0.874/0.764 | **identical** |

So for benchmarking, the `progress_log` EM ETA now covers essentially the whole
runtime, because the dominant export phase is gone.

## Where the EM-internal time goes (ranked)

The model below governs the **EM loop's** cost — it matters for large `r` and for
`autotune` (which fits at many `r`), *not* for the r=10 wall-clock, which is
export-dominated (see the measured breakdown above).

1. **Rigidity expands the state space — the dominant cost.** "Minimum segment
   length = `r` markers" is enforced by expanding each genotype state into a
   chain of ~`r` duration sub-states, so effective states ≈ **3·r**.
   Forward-backward is **quadratic in state count**, so cost scales as **(3r)²**:

   | r | ≈ states | work/marker (∝ states²) | relative to r=10 |
   |---|---|---|---|
   | 5 | 15 | 225 | 0.25× |
   | 10 | 30 | 900 | 1× |
   | 20 | 60 | 3,600 | 4× |
   | 128 | 384 | 147,456 | ~160× |

2. **EM iterates the full sweep many times** (~20–50) until convergence — a
   multiplier on top of factor 1.

3. **Sequence length is the full ~50K markers; low coverage does not help.** The
   ~70% missing positions are kept as uninformative (NA-emission) slots *in the
   chain*, so the HMM still processes ~50K positions per sample.

4. **Samples scale linearly** — 100 samples = 100× one sample, summed in the
   E-step. Real, but the *smallest* multiplier; it is dwarfed by the quadratic
   state expansion (1) × EM iterations (2). **Multiple samples is not the main
   reason a fit is slow.**

This is the *per-iteration* EM cost; in practice EM converges in ~5 iterations
(~10 s at r=10 on the optimized fork — see measured breakdown). The (3r)² scaling
is what makes large-`r` and `autotune` fits expensive, not the r=10 EM loop.

## Why `autotune=TRUE` is much slower

Autotune does not run at a single `r`. It fits the model **repeatedly across a
rigidity ladder** from the initial `rigidity` up toward `max_rigidity`, scoring
each, to pick the optimum. Because high-`r` fits are quadratically expensive
(r=128 ≈ 160× an r=10 fit, per the table), the ladder is dominated by its
largest-`r` rungs. A 100-sample autotune run was killed after ~22 min without
finishing.

## Fits run in this project

| Fit | samples | rigidity setting | outcome |
|---|---|---|---|
| Smoke test | 3 | `autotune=TRUE`, initial r=20 | finished (autotune-chosen r not retained) |
| First 100-sample | 100 | `autotune=TRUE`, initial r=20 | failed: a chromosome's *informative* markers < 2r (sparse low-λ samples) |
| Autotune 100-sample | 100 | `autotune=TRUE`, initial r=10, `max_rigidity=128` | **killed** (too slow; never settled on a final r) |
| Fixed, `save.results=TRUE` | 100 | `autotune=FALSE`, r=10 | ~3 min (dominated by bigWig/PDF export) |
| **Fixed, `save.results=FALSE`** | 100 | `autotune=FALSE`, r=10 | **29 s**, identical scores — the benchmark configuration |

## Scaling to the real data (1,439 samples)

- **Runtime**: ≈14× this 100-sample fit (linear in samples) at the same `r`.
- **Memory**: ≈14× → ~25–30 GB for the joint object (~25 MB/sample at 50K
  markers, plus ~1 GB base). This is the regime Fausto's fork targets.

**Recommendations**
- For benchmarking, run with **`save.results=FALSE`** and score from the `.rds`
  (`@Viterbi`) — skips the dominant per-sample export (~3 min → ~29 s, identical
  scores). Only set `save.results=TRUE` when you actually want the bigWig
  tracks / genotype plots.
- Prefer a **fixed `r`** over autotune; pick `r` from first principles
  ([09](09-rigidity-selection.md)) — autotune's cost rarely justifies it here.
- **Batch samples** to bound peak memory (the HMM parameters are shared, but the
  per-sample matrices accumulate).
- If sweeping `r ∈ {5, 10, 20}`: r=5 is ~4× faster than r=10, r=20 ~4× slower
  (EM cost ∝ (3r)²) — and with `save.results=FALSE` each is seconds-to-minutes.
