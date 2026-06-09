# BC2S3 NIL introgression simulation & RTIGER benchmark

Forward-simulate maize **BC2S3 near-isogenic lines** (teosinte donor × B73
recurrent) with [`simcross`](https://github.com/kbroman/simcross), estimate the
donor introgression-size distribution in Mb, and generate ground-truth datasets
to **benchmark RTIGER** introgression/crossover detection (for the BZea project).

📖 **Docs site & rendered notebook:** <https://sawers-rellan-labs.github.io/zealtiger/>

## Role: one of four ancestry-call sources for the BC2S3 population

This repo is the **theoretical** source — simulated crosses on the consensus
genetic map — giving the expected introgression-size distribution and per-locus
genotype expectations (donor allele frequency **0.125**, heterozygosity **~3%** at
BC2S3; see `docs/02`), plus a **ground-truth benchmark** for evaluating ancestry
callers at the real SNP50K coverage regime. The population will be genotyped four
independent ways and cross-validated:

1. **Theoretical** — this repo (simulation + consensus map).
2. **SNP50K skim (~0.4×) via RTIGER** (optimized fork).
3. **SNP50K skim (~0.4×) via a Gaussian-mixture + HMM caller** (wideseq method).
4. **BRBseq (>10×)** — empirical near-truth.

Sources 2–4 are scored against the simulation truth here (`docs/08`, `docs/11`);
heterozygote-excess and other discrepancies are diagnosed in `docs/12`.

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
├── nil_introgression.qmd       # narrative notebook source (Stages 1-4 + RTIGER benchmark)
├── docs/                       # documentation — browsable at the Pages site
│   ├── README.md                     # docs index
│   ├── 01..05-*.md                   # simulation: map -> sim -> cM->Mb -> size distribution -> crossovers
│   ├── 06..11-*.md                   # RTIGER benchmark: dataset, fit, scoring, rigidity, performance, other tools
│   ├── 12-het-excess-diagnosis.md    # interpreting a heterozygote excess in real ancestry calls
│   └── nil_introgression.html        # RENDERED notebook — the introgression-size distribution figures live here
├── data/        (untracked)    # input: MaizeGDB consensus map
├── results/     (untracked)    # generated output: figures/, *.rds/*.csv, rtiger_benchmark/
└── agent/       (untracked)    # internal roadmap / handoff notes
```

`data/`, `results/`, and `agent/` are git-ignored. See
**[docs/README.md](docs/README.md)** for the pipeline and how to reproduce each
step.

**Where's the introgression-size distribution?** In the rendered notebook
**[`docs/nil_introgression.html`](https://sawers-rellan-labs.github.io/zealtiger/nil_introgression.html)**
(Stage 4 — per-line total donor Mb, per-segment Mb histogram, and the segment-size
ECDF), summarized in **`docs/04-distribution-and-interference.md`**. The figure
files are regenerated to `results/figures/introgression_size.{png,pdf}`.

## Quick start

```bash
Rscript run_all.R --sim          # Stages 1-4 (map QC, 1400-NIL sim, figures)
Rscript make_rtiger_benchmark.R  # 100-NIL RTIGER benchmark (ground truth known)
Rscript fit_rtiger.R 0           # fit RTIGER: all samples, r=3, save.results=FALSE
                                 #   (needs an arm64 Julia; see docs/07, docs/10)
Rscript score_rtiger.R --rtiger results/rtiger_benchmark/rtiger_out_r3/rtiger_result.rds
```

Empirical optimum **r = 3** and the fast `save.results=FALSE` path are explained
in `docs/09`–`docs/10`. Environment: R 4.5.2 (arm64), `simcross` 0.8, `RTIGER`
2.1.0 (`faustovrz/RTIGER@optimize-julia-core`), Julia 1.12 (arm64). Seed `20260609`.
