# Parameter search: what an analytic fragment-size law buys, and why `r` is not an MLE — methodology notes

**Date:** 2026-06-26
**Context:** Whether a closed-form fragment-size distribution could replace the
Monte-Carlo rigidity calibration with a maximum-likelihood search, why the
rigidity `r` cannot be estimated by MLE, and whether ABC MLE (Dean et al. 2014)
applies. Companion to `09-rigidity-selection.md`, `14-interference-and-rigidity.md`,
`15-rigidity-and-scoring.md` (how `r` is derived from a truth sweep), and
`04-distribution-and-interference.md` (the fragment-size law). Memory:
`rigidity-not-mle`, `bc2s3-fragment-size-expectation`, `optimize-R-overrigidifies-nil-sim`.

---

## 1. The original idea

An analytic formula for the fragment-size distribution `f(x)` could, in principle,
accelerate a parameter search: instead of running a Monte-Carlo simulation per
candidate parameter value and scoring the match, you write the **likelihood of the
observed segment sizes** `∏ f(xᵢ; θ)` and maximize it over `θ` directly — one
optimization, no per-candidate simulation.

The question was whether this lets us replace the current rigidity calibration
(grid + golden-section sweep of `r`, minimizing KS distance between caller block
sizes and the simulated truth; `agent/ks_sweep_*.R/.py`) with an MLE.

## 2. A single-generation Gamma fit buys almost nothing

What we have is `Gamma(k=1.15, λ=4.72/Morgan)` fit to the n=1500 BC2S3 simulation
(`results/nil_segments.csv`; see `04-distribution-and-interference.md` and memory
`bc2s3-fragment-size-expectation`). That is `f(x)` at **one** generation, **not**
the family `f(x; g)`.

It only (a) smooths the Null in plots and (b) could turn the `r`-objective from a
two-sample KS (against the noisy n=1500 sample) into a one-sample fit against the
fitted curve. Neither is worth much: n=1500 is already large enough that the
Monte-Carlo noise in the reference is small, and the dominant cost of the `r`-sweep
is the **RTIGER refits at each `r`** (the Julia HMM fits in `evalD`), not the
reference simulation. So a single-generation analytic curve gives essentially **no
speedup**.

The Monte-Carlo search the formula *would* accelerate is the **generation /
effective-meioses inference** — see §4.

## 3. Why `r` cannot be a maximum-likelihood estimate

This is the core point. `r` is **not a parameter of the generative model**, so
there is nothing for an MLE to estimate.

1. **`r` is a property of the caller, not of the data-generating process.** The
   fragment-size law is `f(x; g)` — it contains the number of meioses `g`, no `r`.
   MLE estimates generative parameters from data; there is no likelihood `L(r)`
   over fragment sizes to maximize. There is no "true `r*` that generated the
   reads."

2. **Inside RTIGER's own read-count likelihood, `r` enters as a hard
   regularization constraint, not a free parameter.** Less rigidity = more freedom
   to switch state = weakly **higher** data likelihood. So maximizing the caller's
   likelihood over `r` is monotone toward `r → minimum` — over-segmentation,
   spurious tiny calls. A regularization strength cannot be chosen by the
   likelihood it regularizes; the same reason you cannot pick a lasso `λ` by
   minimizing training error, or the number of clusters by maximizing within-model
   likelihood. (This is also why RTIGER's built-in `optimize_R` over-rigidifies in
   the opposite direction under its iid-LLR model — `09`, `14`, memory
   `optimize-R-overrigidifies-nil-sim`.)

3. **Therefore `r` is set by matching the caller's output size-distribution to an
   external ground truth** (the simulation), via KS distance — a **minimum-distance
   / goodness-of-fit estimator**, not MLE. Even scoring the called sizes under the
   analytic `f` is not a proper `r`-likelihood: called segments are not iid draws
   from `f`; they are a deterministic function of `r` applied to fixed reads.

   *(Aside: a fully-Bayesian marginal-likelihood / empirical-Bayes "evidence" over
   `r` would be a legitimate data-only, Occam-penalized estimator — but RTIGER does
   not expose it, and it is type-II MLE, not the plain MLE imagined.)*

## 4. Where MLE *would* apply — and why it is not worth building now

The legitimate generative parameter is the **generation / effective number of
meioses `g_eff`** (and a recombination-suppression factor) — relevant to the open
observation that all callers run **longer** than the model (mean 34–39 vs 24 cM),
read as recombination suppression / lower effective `g` rather than donor-biased
selection.

To MLE `g_eff` on observed block sizes you need the **family** `Gamma(k(g), λ(g))`:
`λ` scales roughly linearly with the number of meioses, `k>1` from interference +
the haplotype-union. We have one anchor (the fit at the sim's known `g`); a couple
more simcross runs at other `g` would pin down the `g`-dependence. Then a single
MLE (with a confidence interval) replaces a grid of simulations.

**But for the only call needed now — discrete BC2S2 vs BC2S3 (2 simulations;
`bc2s2_vs_bc2s3.qmd`) — Monte Carlo is fine** and the analytic family is not worth
building. Reserve the `f(x; g)` route for if/when a *continuous* `g_eff` /
recombination-suppression estimate is wanted.

## 5. ABC MLE (Dean et al. 2014) — not applicable

`data/dean2014.pdf` — Dean, Singh, Jasra, Peters, *Parameter Estimation for Hidden
Markov Models with Intractable Likelihoods* (Scand. J. Statist. 41:970–987). It is
the **frequentist** ABC variant ("ABC MLE": maximize the ABC likelihood
`P_θ(Y₁∈B^ε_{ŷ₁},…,Yₙ∈B^ε_{ŷₙ})`, estimated by SMC), for HMMs where you can
**simulate but cannot evaluate** the observation density `g_θ(y|x)` (e.g. α-stable
stochastic volatility). It does **not** apply here, for three reasons:

1. **`r` is not a generative parameter.** ABC MLE estimates `θ*` — the parameter
   that generated the data. ABC solves likelihood *intractability*; it does nothing
   about the *non-existence* of a generative likelihood for a regularization knob
   (§3). The `r→min` monotonicity is untouched.
2. **Our emission density is not intractable.** The paper's premise is an
   un-evaluable emission density. Our callers' emissions are Binomial/BetaBinomial
   — fully evaluable. We do not have the problem ABC is for.
3. **For the one place a generative MLE *does* fit (`g_eff`, §4), ABC MLE is still
   Monte Carlo** (SMC) — it would replace simcross-grid Monte Carlo with ABC/SMC
   Monte Carlo (no compute win), and the paper itself proves noisy-ABC's Fisher
   information is **strictly less** than the true MLE's (efficiency loss). ABC is
   the consolation prize for when an analytic likelihood is unavailable, not an
   upgrade.

**Connection worth noting:** the existing KS-vs-simulated-truth `r` calibration is
*already* a minimum-distance / indirect-inference estimator (match a summary
statistic — the block-size distribution — within a tolerance; the
Gourieroux 1993 / Heggland & Frigessi 2004 flavor the paper cites). Dean et al. is
adjacent: it would supply asymptotic consistency/normality theory for that style of
estimator, but gives **no speedup, no way to MLE `r`, and a known efficiency loss**.

## 6. Conclusion

- Keep `r` calibration as **KS-vs-simulated-truth** (minimum-distance). `r` is a
  caller regularization hyperparameter and cannot be MLE'd.
- A single-generation Gamma fit is a convenience for the Null, **not** a route to a
  faster `r` search.
- The only Monte-Carlo search an analytic fragment-size law could replace is the
  **generative `g_eff`** estimate, which needs `Gamma(k(g), λ(g))`; for the
  discrete BC2S2-vs-BC2S3 call (2 sims) it is not worth building.
- ABC MLE (Dean et al. 2014) is the right tool for intractable-emission generative
  MLE — a problem we do not have. Not applicable.
