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

# args: [n_samples] [rigidity]   (n=0 -> all; rigidity default 10, fixed)
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 0L
r <- if (length(args) >= 2) as.integer(args[2]) else 10L

ed <- read_csv("results/rtiger_benchmark/expDesign.csv", show_col_types = FALSE)
if (n > 0L) ed <- ed[seq_len(n), ]
sl <- read_csv("results/rtiger_benchmark/seqlengths.csv", show_col_types = FALSE)
seqv <- setNames(sl$len, sl$chr_label)

outdir <- sprintf("results/rtiger_benchmark/rtiger_out_r%d", r)
if (n > 0L) outdir <- paste0(outdir, "_n", n)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Fixed rigidity (r informative markers; ~r/6.6 Mb resolution at SNP50K coverage,
# see docs/09). autotune is slow/opaque on 100 samples; sweep r instead.
cat(sprintf("Fitting RTIGER on %d samples at rigidity=%d -> %s\n", nrow(ed), r, outdir))
res <- RTIGER(
  expDesign    = ed,
  outputdir    = outdir,
  seqlengths   = seqv,
  rigidity     = r,
  autotune     = FALSE,
  save.results = TRUE,
  verbose      = TRUE
)
saveRDS(res, file.path(outdir, "rtiger_result.rds"))
cat("DONE. samples fitted:", length(res@Viterbi),
    "| result saved to", file.path(outdir, "rtiger_result.rds"), "\n")
