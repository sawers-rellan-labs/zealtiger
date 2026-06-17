#!/usr/bin/env Rscript
# Fit RTIGER per TEOSINTE TAXON on the real SNP50K data (Zx/Zv/Zd/Zl/Zh).
#
# Rationale (refines fit_rtiger_by_donor.R): a single joint EM over all ~1,400 samples
# pools the HMM emission parameters across taxa of differing divergence from B73
# (parviglumis vs diploperennis vs ...), mis-specifying the model; fitting per DONOR
# avoids that but wastes joint strength. Grouping by TAXON is the middle ground —
# samples within a taxon share divergence, so pooling is valid and gains power, while
# the largest fit (Zx ~575) still stays well under a 1,400-sample joint.
#
# Input from make_rtiger_50K_input.sh: data/rtiger_50K/{expDesign_all.csv,seqlengths.csv,
# counts/<donor>/<sample>.tsv}. Taxon = prefix of donor_id before ".".
#
# Output per taxon: data/rtiger_50K/fits_taxa/<taxon>/rtiger_result.rds (+ fit_progress.log)
# plus data/rtiger_50K/fits_taxa/run_summary.csv (taxon,n,seconds,status).
#
# Usage:
#   Rscript fit_rtiger_by_taxa.R [taxon] [rigidity] [save_results] [min_n] [threads]
#     taxon        : run only this taxon (e.g. Zx); "" / omitted = all
#     rigidity     : default 8 (RTIGER autotune pick; see nil_introgression_1400.qmd)
#     save_results : default 0/FALSE
#     min_n        : skip taxa smaller than this (default 1)
#     threads      : Julia EM threads (also set JULIA_NUM_THREADS before startup)
#
# Resource profiling: run ONE taxon under /usr/bin/time -l (see sweep_taxa_50K.sh).

suppressMessages({library(RTIGER); library(readr); library(dplyr)})

JULIA_HOME <- Sys.getenv("JULIA_HOME", "")
if (!nzchar(JULIA_HOME) || !file.exists(file.path(JULIA_HOME, "julia"))) {
  launcher <- path.expand("~/.juliaup/bin/julia")
  if (file.exists(launcher))
    JULIA_HOME <- system(paste(shQuote(launcher), "--startup-file=no -e",
                               shQuote("print(Sys.BINDIR)")), intern = TRUE)
}
setupJulia(JULIA_HOME = JULIA_HOME); sourceJulia()

args       <- commandArgs(trailingOnly = TRUE)
taxon_only <- if (length(args) >= 1 && nzchar(args[1])) args[1] else NA_character_
r          <- if (length(args) >= 2) as.integer(args[2]) else 8L
save_it    <- if (length(args) >= 3) as.logical(as.integer(args[3])) else FALSE
min_n      <- if (length(args) >= 4) as.integer(args[4]) else 1L
threads    <- if (length(args) >= 5) as.integer(args[5]) else 1L

root   <- "data/rtiger_50K"
ed_all <- read_csv(file.path(root, "expDesign_all.csv"), show_col_types = FALSE) |>
  mutate(taxon = sub("[.].*", "", donor))
sl     <- read_csv(file.path(root, "seqlengths.csv"), show_col_types = FALSE)
seqv   <- setNames(sl$len, sl$chr_label)

# Coverage QC: RTIGER aborts if any sample has a chromosome with < 2*rigidity covered
# markers. Drop failed/empty libraries (qc_coverage_50K.sh).
qc_path <- "agent/coverage_qc.csv"
if (file.exists(qc_path)) {
  qc <- read_csv(qc_path, show_col_types = FALSE)
  pass <- qc$n_chr_cov == 10 & qc$min_chr_cov >= 2L * r
  n0 <- nrow(ed_all); ed_all <- ed_all[ed_all$name %in% qc$sample[pass], ]
  cat(sprintf("Coverage QC (r=%d, need >=%d covered/chr): dropped %d/%d libraries.\n",
              r, 2L * r, n0 - nrow(ed_all), n0))
} else cat("WARNING: agent/coverage_qc.csv missing — fitting WITHOUT QC.\n")

sizes <- ed_all |> count(taxon, name = "n") |> arrange(desc(n)) |> filter(n >= min_n)
if (!is.na(taxon_only)) sizes <- sizes |> filter(taxon == taxon_only)
if (nrow(sizes) == 0) stop("No taxon groups to fit (check args / min_n).")

fits_dir <- file.path(root, "fits_taxa"); dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Per-taxon RTIGER: %d taxa, rigidity=%d, threads=%d, save.results=%s\n",
            nrow(sizes), r, threads, save_it))
cat(sprintf("Largest group: %s (n=%d) — the memory/runtime stress fit.\n", sizes$taxon[1], sizes$n[1]))

# Merge each taxon's row into run_summary.csv, REPLACING only its own prior row and
# keeping rows for taxa not in this invocation. This makes per-taxon runs (taxon arg)
# accumulate into one summary over ALL taxa instead of clobbering it down to the last run.
summary_path <- file.path(fits_dir, "run_summary.csv")
for (i in seq_len(nrow(sizes))) {
  tx <- sizes$taxon[i]; n <- sizes$n[i]
  ed <- ed_all |> filter(taxon == tx) |> select(files, name) |> as.data.frame()
  outdir <- file.path(fits_dir, tx); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prog <- file.path(outdir, "fit_progress.log"); if (file.exists(prog)) file.remove(prog)
  cat(sprintf("\n[%d/%d] %s  n=%d  -> %s\n", i, nrow(sizes), tx, n, outdir))
  t0 <- Sys.time(); status <- "ok"
  res <- tryCatch(
    RTIGER(expDesign = ed, outputdir = outdir, seqlengths = seqv, rigidity = r,
           autotune = FALSE, progress_log = prog, save.results = save_it,
           threads = threads, verbose = TRUE),
    error = function(e) { status <<- paste0("ERROR: ", conditionMessage(e)); NULL })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(res)) {
    saveRDS(res, file.path(outdir, "rtiger_result.rds"))
    cat(sprintf("    done in %.1fs (%d samples)\n", secs, length(res@Viterbi)))
  } else cat(sprintf("    FAILED in %.1fs: %s\n", secs, status))
  this_row <- data.frame(taxon = tx, n = n, seconds = round(secs, 1), status = status)
  prior    <- if (file.exists(summary_path))
                read_csv(summary_path, show_col_types = FALSE, col_types = NULL) else this_row[0, ]
  merged   <- bind_rows(prior[!prior$taxon %in% tx, , drop = FALSE], this_row)
  merged   <- merged[order(-merged$n), ]
  write_csv(merged, summary_path)
}
cat(sprintf("\nALL DONE. %d taxa; summary -> %s\n", nrow(sizes), file.path(fits_dir, "run_summary.csv")))
