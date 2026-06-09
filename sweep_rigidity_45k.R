#!/usr/bin/env Rscript
# Sweep RTIGER rigidity on the 45K-array benchmark and pick the F1-optimal r for
# the MolBreeding run (mirrors the SNP50K rigidity selection, docs/09).
#
# Run make_45k_array_benchmark.R first. Output: results/rtiger_45k_benchmark/rigidity_sweep.csv
#
# Usage: Rscript sweep_rigidity_45k.R [r1 r2 ...]   # default 2 3 4 5 6 8 10

suppressMessages({library(tidyverse); library(RTIGER)})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

# JULIA_HOME: env var, else juliaup launcher's Sys.BINDIR (version-independent)
JULIA_HOME <- Sys.getenv("JULIA_HOME", "")
if (!nzchar(JULIA_HOME) || !file.exists(file.path(JULIA_HOME, "julia"))) {
  launcher <- path.expand("~/.juliaup/bin/julia")
  if (file.exists(launcher))
    JULIA_HOME <- system(paste(shQuote(launcher), "--startup-file=no -e",
                               shQuote("print(Sys.BINDIR)")), intern = TRUE)
}
setupJulia(JULIA_HOME = JULIA_HOME); sourceJulia()

args <- commandArgs(trailingOnly = TRUE)
rs   <- if (length(args) >= 1) as.integer(args) else c(2L,3L,4L,5L,6L,8L,10L)

dir <- "results/rtiger_45k_benchmark"
ed  <- read_csv(file.path(dir, "expDesign.csv"), show_col_types = FALSE)
sl  <- read_csv(file.path(dir, "seqlengths.csv"), show_col_types = FALSE)
seqv       <- setNames(sl$len, sl$chr_label)
truth_seg  <- read_csv(file.path(dir, "truth_segments.csv"), show_col_types = FALSE)
truth_mark <- readRDS(file.path(dir, "truth_markers.rds"))

rows <- list()
for (r in rs) {
  od <- file.path(dir, sprintf("fit_r%d", r)); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  fit <- tryCatch(
    RTIGER(ed[, c("files","name")], outputdir = od, seqlengths = seqv,
           rigidity = r, autotune = FALSE, save.results = FALSE, verbose = FALSE),
    error = function(e) { message(sprintf("r=%d failed: %s", r, conditionMessage(e))); NULL })
  if (is.null(fit)) next
  called <- segments_from_rtiger_object(fit)
  sm <- score_markers(assign_called_state(called, truth_mark))
  f1 <- function(t) score_segments(called, truth_seg, target = t)$summary$F1
  row <- tibble(r = r, marker_acc = sm$accuracy,
                F1_donor_present = f1(c(1L,2L)), F1_hom_donor = f1(2L), F1_het = f1(1L))
  rows[[as.character(r)]] <- row
  cat(sprintf("r=%2d  marker_acc=%.4f  F1_donor_present=%.3f  F1_hom_donor=%.3f  F1_het=%.3f\n",
              r, row$marker_acc, row$F1_donor_present, row$F1_hom_donor, row$F1_het))
}

tab <- bind_rows(rows)
readr::write_csv(tab, file.path(dir, "rigidity_sweep.csv"))
best <- tab$r[which.max(tab$F1_donor_present)]
cat(sprintf("\nBest rigidity (max donor-present F1): r = %d  →  use for the MolBreeding run.\n", best))
cat(sprintf("Sweep table → %s\n", file.path(dir, "rigidity_sweep.csv")))
