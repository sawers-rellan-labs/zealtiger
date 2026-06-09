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

**A-priori estimate (superseded by the empirical sweep below): `r ≈ 8–10`.**
Reasoning at the time: ~8–10 markers seemed the minimum to separate het from
homozygous at ~1.4-read depth while keeping false positives negligible, and
`2r ≤` the sparsest chromosome's informative count (2×10 ≪ 441). The measured
sweep shows this was **too conservative** — see next section.

## Empirical optimum — the rigidity sweep (this is the answer)

Fitting the 100-NIL benchmark at `r ∈ {2,3,5,10,20}` (fixed, `save.results=FALSE`,
~30 s each) and scoring against the ground truth ([08](08-scoring.md)) gives the
**F1 by class**:

| r | homo-donor F1 | donor-present F1 | het F1 | CO bias (true 21.8) |
|---|---|---|---|---|
| **2**  | **0.953** | **0.922** | **0.858** | −2.38 |
| **3**  | 0.952 | 0.921 | 0.858 | −2.44 |
| 5  | 0.949 | 0.920 | 0.851 | −2.58 |
| 10 | 0.930 | 0.912 | 0.815 | −2.91 |
| 20 | 0.872 | 0.873 | 0.717 | −3.97 |

**F1 rises monotonically as `r` drops, then plateaus at `r ≈ 2–3`** (r=2 vs r=3
differ by ≤0.001 in every class). The knee is sharp — most of the gain is already
in by r=5 — and there is **no measurable benefit below r=3**.

Why the a-priori "8–10" was too high: it assumed FP control *required* ~8–10
markers. Empirically **precision does not break down at low r** — at r=2 it is
still 0.99 (homo-donor) / 0.98 (donor-present) / **0.96 (het)**. So the recall
gained by lowering `r` is nearly free, and the optimum is recall-driven and low.
The opposite failure shows up at *high* r: r=20 forces long segments and
*manufactures* spurious het (precision 0.73). Crossover recovery is likewise best
at low r (bias −2.4 at r=2–3 vs −4.0 at r=20).

**Recommendation: `r = 3`** — on the plateau, round, and robust (anything in 2–5
is within ~0.01 F1). `2r = 6 ≪ 441` sparsest-chr informative markers, so it is
always admissible.

## The coverage floor (a result, not a tuning failure)

The paper's "1–2% missed" target is **unreachable at SNP50K coverage**. Their
dense WGS had ~4,150 informative markers/Mb (`r=250` → 60 kb at ~1% FN); SNP50K
has ~630× lower density (6.6/Mb), so **even `r=5` misses ~6–9% of
introgressions** — everything below ~0.5–0.8 Mb. There is a hard detection floor
around **0.5–1.5 Mb** at 0.43×/70% missing; sub-floor introgressions are lost for
any `r`. The benchmark quantifies exactly this: at the optimum `r=3`, ~8% of
homozygous-donor and ~13% of donor-present segments are still missed — and the
F1 plateau below r=3 confirms these are the floor (lowering `r` further does not
recover them), not a tuning shortfall.

Also: a fixed `r` gives *variable physical resolution* across samples (3.6 Mb at
λ=0.15, 0.6 Mb at λ=1.4). For uniform, conservative FP control, calibrate `r` to
the lowest-coverage sample you intend to keep.

## Recommended procedure

Don't trust a single autotuned number — **sweep** it. Fit at fixed `r` over a
range (`autotune = FALSE`, `save.results = FALSE` → ~30 s each) and score against
the ground truth ([08](08-scoring.md)); pick the `r` that maximizes F1. Here that
sweep landed on **`r = 3`** (plateau 2–3), well below the a-priori guess — a
reminder to measure rather than assume. For a *new* coverage/marker regime, rerun
the sweep: the optimal `r` tracks informative-marker density, so it will differ.
