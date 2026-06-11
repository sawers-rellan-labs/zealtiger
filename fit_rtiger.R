#!/usr/bin/env Rscript
# Fit RTIGER on the simulated benchmark and save results.
#
# Usage:
#   Rscript fit_rtiger.R [n_samples]   # n_samples omitted/0 = all
#
# Writes per-sample results under results/rtiger_benchmark/rtiger_out/ and the
# fitted object to rtiger_out/rtiger_result.rds (score with score_rtiger.R).

suppressMessages({library(RTIGER); library(readr)})

JULIA_HOME <- Sys.getenv("JULIA_HOME",
  "/Applications/Julia-1.11.app/Contents/Resources/julia/bin")
setupJulia(JULIA_HOME = JULIA_HOME)
sourceJulia()

# args: [n_samples] [rigidity] [save_results] [threads]
#   n=0 -> all samples; rigidity default 3 (minimum-CO-bias point of the sweep,
#   docs/09 — autotune objectively selects r=8, see agent/autotune_rtiger_benchmark.R);
#   save_results default 0/FALSE; threads default 1.
# save_results=FALSE skips the dominant per-sample export (bigWig tracks + genotype
# PDFs ~ minutes for 100 samples) — those are visualization artifacts, not the
# analysis. The fitted object is always saved to rtiger_result.rds, and the
# scorer reads segments straight from it, so benchmark runs need save_results=0.
# threads: RTIGER's Julia core parallelizes per-sample EM; threads=4 cuts walltime
# ~35% on this 100-NIL benchmark. Results are deterministic in the thread count.
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 0L
r <- if (length(args) >= 2) as.integer(args[2]) else 3L
save_it <- if (length(args) >= 3) as.logical(as.integer(args[3])) else FALSE
threads <- if (length(args) >= 4) as.integer(args[4]) else 1L

ed <- read_csv("results/rtiger_benchmark/expDesign.csv", show_col_types = FALSE)
if (n > 0L) ed <- ed[seq_len(n), ]
sl <- read_csv("results/rtiger_benchmark/seqlengths.csv", show_col_types = FALSE)
seqv <- setNames(sl$len, sl$chr_label)

outdir <- sprintf("results/rtiger_benchmark/rtiger_out_r%d", r)
if (n > 0L) outdir <- paste0(outdir, "_n", n)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Fixed rigidity (r informative markers; ~r/6.6 Mb resolution at SNP50K coverage,
# see docs/09). autotune is slow/opaque on 100 samples; sweep r instead.
# progress_log: per-EM-iteration ETA log (faustovrz/RTIGER@optimize-julia-core);
# tail it for a live ETA.
prog <- file.path(outdir, "fit_progress.log")
if (file.exists(prog)) file.remove(prog)
cat(sprintf("Fitting RTIGER on %d samples at rigidity=%d, threads=%d (save.results=%s) -> %s\n",
            nrow(ed), r, threads, save_it, outdir))
cat(sprintf("Progress log (tail for ETA): %s\n", prog))
res <- RTIGER(
  expDesign    = ed,
  outputdir    = outdir,
  seqlengths   = seqv,
  rigidity     = r,
  autotune     = FALSE,
  progress_log = prog,
  save.results = save_it,
  threads      = threads,
  verbose      = TRUE
)
saveRDS(res, file.path(outdir, "rtiger_result.rds"))
cat("DONE. samples fitted:", length(res@Viterbi),
    "| result saved to", file.path(outdir, "rtiger_result.rds"), "\n")
