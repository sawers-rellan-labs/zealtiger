# Choosing the RTIGER rigidity parameter

## Units: rigidity is in *informative* markers

RTIGER's `rigidity` (min markers per segment) counts **informative** markers
(covered, non-`0/0`), not the 50K input positions. Evidence: the fit aborted when
a chromosome had 16 *informative* markers vs `2×rigidity=40`, despite ~3,585 total
positions ([07](07-rtiger-fitting.md)). So `r` must be reasoned about in
informative-marker units.

## Derivation from the generated data

Informative density at SNP50K (λ=0.43, ~71% missing): **~6.6 informative
markers/Mb** (per-sample 2.8–16.5). So `r` informative markers ⇒ minimum
detectable segment ≈ `r / 6.6` Mb. Following the paper's rule (set `r` as high as
possible while the missed-segment fraction stays ~1–2%):

| r (inf. markers) | min detectable | homo-donor missed | het missed |
|---|---|---|---|
| 5 | ~0.76 Mb | 6.2% | 14.3% |
| 8 | ~1.2 Mb | 10.2% | 19.6% |
| **10** | **~1.5 Mb** | **12.4%** | **24.8%** |
| 15 | ~2.3 Mb | 18.4% | 33.1% |
| 20 | ~3.0 Mb | 23.2% | 40.1% |

**Recommendation: `r ≈ 8–10` informative markers.**
- Captures the bulk of introgressions (median ~8 Mb).
- Strong FP control: at error 0.005 a spurious run of ≥8 same-error markers is
  essentially impossible, and ~8–10 markers is the minimum to distinguish het
  from homozygous at ~1.4-read depth.
- Satisfies `2r ≤` the sparsest chromosome's informative count (2×10 ≪ 441).

Lower `r=5` favors recall of small donor blocks (more FP); higher `r=15` favors
clean FDR/het calls.

## The coverage floor (a result, not a tuning failure)

The paper's "1–2% missed" target is **unreachable at SNP50K coverage**. Their
dense WGS had ~4,150 informative markers/Mb (`r=250` → 60 kb at ~1% FN); SNP50K
has ~630× lower density (6.6/Mb), so **even `r=5` misses ~6–9% of
introgressions** — everything below ~0.5–0.8 Mb. There is a hard detection floor
around **0.5–1.5 Mb** at 0.43×/70% missing; sub-floor introgressions are lost for
any `r`. The benchmark quantifies exactly this.

Also: a fixed `r` gives *variable physical resolution* across samples (3.6 Mb at
λ=0.15, 0.6 Mb at λ=1.4). For uniform, conservative FP control, calibrate `r` to
the lowest-coverage sample you intend to keep.

## Recommended procedure

Don't trust a single autotuned number — **sweep** it. Fit at `r ∈ {5, 10, 20}`
(fixed, `autotune = FALSE`, fast) and score each against the ground truth
([08](08-scoring.md)); the right `r` is the one maximizing F1 on the benchmark
(expected ~8–12, consistent with the derivation above).
