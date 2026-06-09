# Fitting RTIGER on the benchmark

**Code:** `fit_rtiger.R`.
**Output:** `results/rtiger_benchmark/rtiger_out/` — per-sample
`CompleteBlock-state-<name>.bed` (chr, start, end, `AA`/`AB`/`BB`), state/count
tracks, plots, and `rtiger_result.rds` (the fitted S4 object).

## Package: the fork, not CRAN

Install **`faustovrz/RTIGER`** (Fausto's time/memory-optimized fork), v2.1.0:

```r
remotes::install_github("faustovrz/RTIGER")   # remove any CRAN RTIGER first
```

RTIGER's HMM/Viterbi runs on a **Julia backend** via `JuliaCall`
(`setupJulia()` + `sourceJulia()`). `generateObject()` loads data without Julia
(handy for a format check); the fit needs Julia.

## Julia architecture gotcha (Apple Silicon)

JuliaCall loads `libjulia` **into the R process**, so Julia must match R's
architecture. R here is **arm64**; the `/Applications/Julia-1.11.app` install is
**x86_64** and failed with `incompatible architecture … have 'x86_64', need
'arm64'`. Fix: use the **arm64 Julia** managed by `juliaup`
(`~/.julia/juliaup/julia-1.12.6+…aarch64…`) and pass its `bin` dir as
`JULIA_HOME`. `fit_rtiger.R` defaults to that path (override via the `JULIA_HOME`
env var). First run precompiles the Julia packages (a few minutes).

## State mapping

RTIGER labels: `AA`/`pat` = homozygous P1 = REF = recurrent → state **0**;
`AB`/`het` → **1**; `BB`/`mat` = homozygous P2 = ALT = donor → **2**.

## The `2 × rigidity` failure (and fix)

The first 100-sample fit aborted: *"Some of your observations is smaller than 2
times rigidity."* Cause: **RTIGER's rigidity is counted in *informative*
markers** (covered, non-`0/0`), and the unrealistic λ≈0.01 samples (from the raw
Gamma lower tail) had as few as 16 informative markers on a chromosome — below
`2 × rigidity`. Fix:
1. Floor per-sample λ at 0.15 in the benchmark ([06](06-rtiger-benchmark-dataset.md)),
   removing samples real genotyping would exclude — the sparsest now has ~441
   informative markers on its thinnest chromosome.
2. In `fit_rtiger.R`, lower the initial `rigidity` to 10 and cap
   `max_rigidity = 128` (so `2 × max_rigidity = 256 <` the sparsest chromosome).

## Run

```bash
JULIA_HOME=<arm64 julia bin> Rscript fit_rtiger.R 0   # 0 = all samples; N = first N
```

`autotune = TRUE` is opaque and slow on 100 samples (an initial joint fit, then a
rigidity ladder, then a final fit; ~25 MB/sample, single-threaded, no progress
bar). For benchmarking, prefer a **fixed rigidity** (`autotune = FALSE`,
`rigidity = 10`) and/or a sweep over r — see [09](09-rigidity-selection.md).
A 3-sample smoke fit validates the toolchain in seconds.
