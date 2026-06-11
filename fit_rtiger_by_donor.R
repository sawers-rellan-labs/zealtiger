#!/usr/bin/env Rscript
# Fit RTIGER per donor accession on the real SNP50K data.
#
# Rationale: a single joint EM over all 1,439 samples pools the HMM emission /
# transition parameters across donors that differ in sequence divergence from
# B73 (luxurians vs diploperennis vs parviglumis ...), which mis-specifies the
# model. Splitting the joint fit BY DONOR ACCESSION removes cross-donor pooling
# while keeping joint strength within each donor. It also caps the largest single
# fit at ~100 samples (Zx.0580) instead of 1,439 -> runs comfortably on 24 GB.
#
# Input is produced by make_rtiger_50K_input.sh:
#   data/rtiger_50K/expDesign_all.csv   (files,name,donor)
#   data/rtiger_50K/seqlengths.csv      (chr_label,len)
#   data/rtiger_50K/counts/<donor>/<sample>.tsv
#
# Output, per accession:
#   data/rtiger_50K/fits/<donor>/rtiger_result.rds   (fitted object; score from it)
#   data/rtiger_50K/fits/<donor>/fit_progress.log    (per-EM-iteration ETA)
# plus data/rtiger_50K/fits/run_summary.csv          (donor,n,seconds,status).
#
# Usage:
#   Rscript fit_rtiger_by_donor.R [donor] [rigidity] [save_results] [min_n] [threads]
#     donor        : run only this accession (e.g. Zx.0580); "" / omitted = all
#     rigidity     : default 8 (RTIGER autotune pick on the SNP50K simulation, and
#                    the value both empirical platforms autotune to; see
#                    agent/autotune_rtiger_benchmark.R and nil_introgression.qmd)
#     save_results : default 0/FALSE (skip per-sample bigWig/PDF export)
#     min_n        : skip donor groups smaller than this (default 1 = fit all)
#     threads      : Julia EM threads (default 1; 4 cuts walltime ~35% per group)
#
# To capture peak RSS for the resource table, run the whole thing under:
#   /usr/bin/time -l Rscript fit_rtiger_by_donor.R          # macOS -> "maximum resident set size"
# or fit the largest group alone first:
#   /usr/bin/time -l Rscript fit_rtiger_by_donor.R Zx.0580

suppressMessages({library(RTIGER); library(readr); library(dplyr)})

# Resolve JULIA_HOME: env var if it points at a real julia, else ask the juliaup
# launcher for its Sys.BINDIR (version-independent — survives Julia upgrades).
JULIA_HOME <- Sys.getenv("JULIA_HOME", "")
if (!nzchar(JULIA_HOME) || !file.exists(file.path(JULIA_HOME, "julia"))) {
  launcher <- path.expand("~/.juliaup/bin/julia")
  if (file.exists(launcher)) {
    JULIA_HOME <- system(paste(shQuote(launcher), "--startup-file=no -e",
                               shQuote("print(Sys.BINDIR)")), intern = TRUE)
  }
}
setupJulia(JULIA_HOME = JULIA_HOME)
sourceJulia()

args        <- commandArgs(trailingOnly = TRUE)
donor_only  <- if (length(args) >= 1 && nzchar(args[1])) args[1] else NA_character_
r           <- if (length(args) >= 2) as.integer(args[2]) else 8L
save_it     <- if (length(args) >= 3) as.logical(as.integer(args[3])) else FALSE
min_n       <- if (length(args) >= 4) as.integer(args[4]) else 1L
threads     <- if (length(args) >= 5) as.integer(args[5]) else 1L

root    <- "data/rtiger_50K"
ed_all  <- read_csv(file.path(root, "expDesign_all.csv"), show_col_types = FALSE)
sl      <- read_csv(file.path(root, "seqlengths.csv"),    show_col_types = FALSE)
seqv    <- setNames(sl$len, sl$chr_label)

# Coverage QC: RTIGER aborts a joint fit if ANY sample has a chromosome with
# < 2*rigidity covered markers. Drop failed/empty libraries (see qc_coverage_50K.sh).
qc_path <- "agent/coverage_qc.csv"
if (file.exists(qc_path)) {
  qc <- read_csv(qc_path, show_col_types = FALSE)
  qc_pass <- qc$n_chr_cov == 10 & qc$min_chr_cov >= 2L * r
  dropped <- qc$sample[!qc_pass]
  n0 <- nrow(ed_all)
  ed_all <- ed_all[!(ed_all$name %in% dropped), ]
  cat(sprintf("Coverage QC (r=%d, need >=%d covered/chr on all 10 chr): dropped %d/%d failed libraries.\n",
              r, 2L * r, n0 - nrow(ed_all), n0))
} else {
  cat("WARNING: agent/coverage_qc.csv missing — fitting WITHOUT QC; a failed library may abort its donor group.\n")
}

# Group sizes; biggest first so the memory-stress fit (and any OOM) shows early.
sizes   <- ed_all |> count(donor, name = "n") |> arrange(desc(n))
if (!is.na(donor_only)) sizes <- sizes |> filter(donor == donor_only)
sizes   <- sizes |> filter(n >= min_n)
if (nrow(sizes) == 0) stop("No donor groups to fit (check args / min_n).")

fits_dir <- file.path(root, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Per-donor RTIGER: %d accession(s), rigidity=%d, threads=%d, save.results=%s, min_n=%d\n",
            nrow(sizes), r, threads, save_it, min_n))
cat(sprintf("Largest group: %s (n=%d) — the memory/runtime stress fit.\n",
            sizes$donor[1], sizes$n[1]))

summary_rows <- vector("list", nrow(sizes))

for (i in seq_len(nrow(sizes))) {
  d   <- sizes$donor[i]
  n   <- sizes$n[i]
  ed  <- ed_all |> filter(donor == d) |> select(files, name) |> as.data.frame()

  outdir <- file.path(fits_dir, d)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prog <- file.path(outdir, "fit_progress.log")
  if (file.exists(prog)) file.remove(prog)

  cat(sprintf("\n[%d/%d] %s  n=%d  -> %s\n", i, nrow(sizes), d, n, outdir))
  t0 <- Sys.time()
  status <- "ok"
  res <- tryCatch(
    RTIGER(
      expDesign    = ed,
      outputdir    = outdir,
      seqlengths   = seqv,
      rigidity     = r,
      autotune     = FALSE,
      progress_log = prog,
      save.results = save_it,
      threads      = threads,
      verbose      = TRUE
    ),
    error = function(e) { status <<- paste0("ERROR: ", conditionMessage(e)); NULL }
  )
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (!is.null(res)) {
    saveRDS(res, file.path(outdir, "rtiger_result.rds"))
    cat(sprintf("    done in %.1fs (%d samples fitted)\n", secs, length(res@Viterbi)))
  } else {
    cat(sprintf("    FAILED in %.1fs: %s\n", secs, status))
  }
  summary_rows[[i]] <- data.frame(donor = d, n = n, seconds = round(secs, 1),
                                  status = status, stringsAsFactors = FALSE)
  # Write the summary incrementally so progress survives an interrupt.
  write_csv(dplyr::bind_rows(summary_rows[seq_len(i)]),
            file.path(fits_dir, "run_summary.csv"))
}

cat(sprintf("\nALL DONE. %d accessions; summary -> %s\n",
            nrow(sizes), file.path(fits_dir, "run_summary.csv")))
