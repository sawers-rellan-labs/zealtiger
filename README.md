# BC2S3 NIL introgression simulation & RTIGER benchmark

Forward-simulate maize **BC2S3 near-isogenic lines** (teosinte donor × B73
recurrent) with [`simcross`](https://github.com/kbroman/simcross), estimate the
donor introgression-size distribution in Mb, and generate ground-truth datasets
to **benchmark RTIGER** introgression/crossover detection (for the BZea project).

## Layout

```
.
├── R/                  # function libraries (sourced by the scripts/notebook)
│   ├── 01_compile_map.R        # map parsing + monotone Marey map + thinning
│   ├── 02_simulate.R           # pedigree, mosaic extraction, crossovers
│   ├── 03_segments_to_mb.R     # cM -> Mb interpolation
│   ├── 04_summarize.R          # distribution, CO summary, assertions
│   ├── 05_make_rtiger_input.R  # simulate genotypes -> RTIGER allele counts
│   └── 06_score_rtiger.R       # score RTIGER calls vs ground truth
├── run_all.R                   # Stages 1-4: map -> simulate -> Mb -> figures
├── make_rtiger_benchmark.R     # generate the 100-NIL RTIGER benchmark
├── fit_rtiger.R                # fit RTIGER (Julia backend) on the benchmark
├── score_rtiger.R              # score a fit against ground truth
├── nil_introgression.qmd       # narrative notebook (Stages 1-4)
├── docs/                       # step-by-step documentation (start at docs/README.md)
├── data/        (untracked)    # input: MaizeGDB consensus map
└── results/     (untracked)    # all generated output
```

`data/` and `results/` are git-ignored. See **[docs/README.md](docs/README.md)**
for the pipeline, decisions, and how to reproduce each step.

## Quick start

```bash
Rscript run_all.R --sim          # Stages 1-4 (map QC, 1400-NIL sim, figures)
Rscript make_rtiger_benchmark.R  # 100-NIL RTIGER benchmark (ground truth known)
Rscript fit_rtiger.R 0           # fit RTIGER (needs an arm64 Julia; see docs/07)
Rscript score_rtiger.R --rtiger results/rtiger_benchmark/rtiger_out
```

Environment: R 4.5.2 (arm64), `simcross` 0.8, `RTIGER` 2.1.0
(`faustovrz/RTIGER` fork), Julia 1.12 (arm64). Seed `20260609`.
