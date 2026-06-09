# Scoring RTIGER calls against ground truth

**Code:** `R/06_score_rtiger.R`, driven by `score_rtiger.R`.
**Output:** `results/rtiger_benchmark/scoring/` — `marker_scores.csv`,
`segment_scores.csv`, `crossover_scores.csv`, `crossover_recovery.png`.

## Inputs

```bash
Rscript score_rtiger.R --rtiger results/rtiger_benchmark/rtiger_out   # dir of CompleteBlock-*.bed
Rscript score_rtiger.R --rtiger <results.rds>                          # a saved RTIGER object
Rscript score_rtiger.R --demo                                         # self-test, no RTIGER needed
```

`read_rtiger_segments` reads either the per-sample `CompleteBlock-state-*.bed`
files or a fitted RTIGER object (segments extracted from its `@Viterbi` slot),
normalizing to a tidy table (`name, chr, start_bp, end_bp, state` in 0/1/2). State
mapping as in [07](07-rtiger-fitting.md).

## Metrics

1. **Marker-level** (`score_markers`) — overlap-join each truth marker to the
   called state; confusion matrix (true × called over 0/1/2), per-state
   recall/precision, overall accuracy, uncovered fraction. Resolution-independent.
2. **Segment-level** (`score_segments`) — precision / recall / F1 / FDR by
   **reciprocal overlap ≥ 0.5**, for three event classes: **homozygous-donor**
   (state 2), **donor-present** (states 1+2), **heterozygous** (state 1); plus
   median start/end boundary error (kb).
3. **Crossovers** (`score_cos`) — true vs called COs per line (bias, RMSE) and a
   `crossover_recovery.png` scatter.

## Results — RTIGER on the SNP50K benchmark (r = 3)

Best-rigidity fit (`r=3`, the empirical optimum from the sweep —
[09](09-rigidity-selection.md)); 100 NILs at λ=0.43×, ~71% missing;
`save.results=FALSE`, scored from the `.rds`.

**Marker-level:** overall accuracy **0.999**, 0.36% uncovered.

| state | recall | precision |
|---|---|---|
| 0 homozygous recurrent | 1.000 | 0.999 |
| 1 heterozygous | 0.974 | 0.990 |
| 2 homozygous donor | 0.998 | 0.997 |

**Segment-level** (reciprocal overlap ≥ 0.5):

| class | n_truth | precision | recall | F1 | FDR | boundary err |
|---|---|---|---|---|---|---|
| homozygous donor | 884 | 0.987 | 0.920 | **0.952** | 0.013 | ~85 kb |
| donor present (het+homo) | 1400 | 0.981 | 0.869 | **0.921** | 0.019 | ~85 kb |
| heterozygous | 516 | 0.962 | 0.775 | **0.858** | 0.039 | ~127 kb |

**Crossovers:** true 21.8 vs called 19.4 per line (bias −2.4, RMSE 3.3).

Takeaways: marker-level genotyping is near-perfect (**99.9%**); donor-segment
**precision ≥0.98, FDR ≤2%** — almost no false introgressions; **recall ~0.87–0.92**,
the missing ~8–13% being sub-Mb segments below the coverage floor
([09](09-rigidity-selection.md)), not a tuning shortfall. Boundaries land within
~2 marker spacings (~85 kb; ~127 kb for het). The full `r∈{2,3,5,10,20}` sweep
that selects `r=3` is in [09](09-rigidity-selection.md).

## Self-test (validates the scorer with no Julia)

`--demo` scores two synthetic call sets against truth:

| call set | homo-donor recall | donor-present recall | CO bias | start err |
|---|---|---|---|---|
| perfect (= truth) | 1.000 | 1.000 | 0.0 | 0 kb |
| degraded mock | 0.889 | 0.878 | −2.5 | 22.7 kb |

Perfect input scores a perfect 1.0 everywhere (no logic bug); the degraded mock
(small segments merged, boundaries jittered) drops recall and gains boundary
error exactly as expected. `degrade_truth` builds the mock.

## Interpreting against the coverage floor

At SNP50K coverage there is a hard detection floor (~0.5–1.5 Mb); sub-floor
introgressions are systematically missed regardless of RTIGER settings. Read the
segment recall numbers together with [09](09-rigidity-selection.md).
