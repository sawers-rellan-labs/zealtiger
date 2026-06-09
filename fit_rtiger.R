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

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 0L

ed <- read_csv("results/rtiger_benchmark/expDesign.csv", show_col_types = FALSE)
if (n > 0L) ed <- ed[seq_len(n), ]
sl <- read_csv("results/rtiger_benchmark/seqlengths.csv", show_col_types = FALSE)
seqv <- setNames(sl$len, sl$chr_label)

outdir <- "results/rtiger_benchmark/rtiger_out"
if (n > 0L) outdir <- paste0(outdir, "_n", n)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Fitting RTIGER on %d samples -> %s\n", nrow(ed), outdir))
res <- RTIGER(
  expDesign        = ed,
  outputdir        = outdir,
  seqlengths       = seqv,
  rigidity         = 10,      # initial guess; autotune refines (low for sparse data)
  autotune         = TRUE,
  max_rigidity     = 128,     # cap so 2*R stays below the sparsest chr's markers
  average_coverage = 0.43,    # SNP50K mean (known)
  save.results     = TRUE,
  verbose          = TRUE
)
saveRDS(res, file.path(outdir, "rtiger_result.rds"))
cat("DONE. samples fitted:", length(res@Viterbi),
    "| result saved to", file.path(outdir, "rtiger_result.rds"), "\n")
