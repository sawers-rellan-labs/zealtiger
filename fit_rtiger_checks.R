#!/usr/bin/env Rscript
# Fit RTIGER on the CHECK samples (B73, Purple) the production per-taxon run skips.
#
# The checks have no teosinte donor, so they don't belong to any taxon group.
# We pool by CHECK GROUP instead (Purple together, B73 together) — Purples share
# one emission regime (~90% non-B73), B73s another (pure recurrent). Same RTIGER
# settings as fit_rtiger_by_taxa.R: autotune=FALSE, fixed rigidity.
#
# Input : data/rtiger_50K/{expDesign_checks.csv (files,name,group), seqlengths.csv}
#         from extract_check_counts_50K.sh.
# Output: data/rtiger_50K/fits_checks/<group>/rtiger_result.rds
# Then  : Rscript rtiger_to_common_schema.R data/rtiger_50K/fits_checks \
#               data/rtiger_50K/calls_checks_common_schema.csv
#
# Usage: Rscript fit_rtiger_checks.R [rigidity=5] [group="" (all)] [threads=1]

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
r          <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 5L
group_only <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NA_character_
threads    <- if (length(args) >= 3) as.integer(args[3]) else 1L

root   <- "data/rtiger_50K"
ed_all <- read_csv(file.path(root, "expDesign_checks.csv"), show_col_types = FALSE)
sl     <- read_csv(file.path(root, "seqlengths.csv"), show_col_types = FALSE)
seqv   <- setNames(sl$len, sl$chr_label)

# Coverage QC: RTIGER aborts if any sample has a chromosome with < 2*rigidity
# covered markers. Drop libraries that fail (same guard as the taxa fit).
covered_ok <- function(file) {
  d <- suppressMessages(read_tsv(file, col_names = FALSE, show_col_types = FALSE))
  cov <- d[d$X4 + d$X6 > 0, ]
  per_chr <- table(cov$X1)
  length(per_chr) == 10 && min(per_chr) >= 2L * r
}
keep <- vapply(ed_all$files, covered_ok, logical(1))
if (any(!keep))
  cat(sprintf("Coverage QC (r=%d): dropping %d/%d check libraries: %s\n",
              r, sum(!keep), nrow(ed_all), paste(ed_all$name[!keep], collapse = ", ")))
ed_all <- ed_all[keep, ]

groups <- ed_all |> count(group, name = "n") |> arrange(desc(n))
if (!is.na(group_only)) groups <- groups |> filter(group == group_only)
if (nrow(groups) == 0) stop("No check groups to fit.")

fits_dir <- file.path(root, "fits_checks")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Per-group RTIGER on checks: %d group(s), rigidity=%d, threads=%d\n",
            nrow(groups), r, threads))

for (i in seq_len(nrow(groups))) {
  g <- groups$group[i]; n <- groups$n[i]
  ed <- ed_all |> filter(group == g) |> select(files, name) |> as.data.frame()
  outdir <- file.path(fits_dir, g); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prog <- file.path(outdir, "fit_progress.log"); if (file.exists(prog)) file.remove(prog)
  cat(sprintf("\n[%d/%d] %s  n=%d  -> %s\n", i, nrow(groups), g, n, outdir))
  t0 <- Sys.time(); status <- "ok"
  res <- tryCatch(
    RTIGER(expDesign = ed, outputdir = outdir, seqlengths = seqv, rigidity = r,
           autotune = FALSE, progress_log = prog, save.results = FALSE,
           threads = threads, verbose = TRUE),
    error = function(e) { status <<- paste0("ERROR: ", conditionMessage(e)); NULL })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(res)) {
    saveRDS(res, file.path(outdir, "rtiger_result.rds"))
    cat(sprintf("    done in %.1fs (%d samples)\n", secs, length(res@Viterbi)))
  } else cat(sprintf("    FAILED in %.1fs: %s\n", secs, status))
}
cat("\nALL DONE. Convert with:\n  Rscript rtiger_to_common_schema.R",
    file.path(root, "fits_checks"), file.path(root, "calls_checks_common_schema.csv"), "\n")
