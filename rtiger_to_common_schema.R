#!/usr/bin/env Rscript
# Convert the per-donor RTIGER fits into the common cross-source call schema.
#
# Reads every  data/rtiger_50K/fits/<donor>/rtiger_result.rds  and emits one
# tidy table of ancestry segments in the schema shared by all four call sources
# (simulation truth, RTIGER, method-B GMM/HMM, BRBseq):
#
#   source, donor, name, chr, start_bp, end_bp, state
#     state: 0 = homozygous recurrent (B73/AA), 1 = het (AB), 2 = homozygous donor (BB/teo)
#
# Reuses the tested extractor segments_from_rtiger_object() from R/06_score_rtiger.R
# (handles the in-memory object, so it works with save.results=FALSE — no BEDs).
# Idempotent and incremental: run it now over the fits done so far, and again
# when all 82 finish.
#
# Usage:  Rscript rtiger_to_common_schema.R [fits_dir] [out_csv]

suppressMessages({
  library(tidyverse); library(GenomicRanges); library(S4Vectors)
})
source("R/06_score_rtiger.R")   # defines segments_from_rtiger_object(), .vit_map

args     <- commandArgs(trailingOnly = TRUE)
fits_dir <- if (length(args) >= 1) args[1] else "data/rtiger_50K/fits"
out_csv  <- if (length(args) >= 2) args[2] else "data/rtiger_50K/calls_common_schema.csv"

rds <- Sys.glob(file.path(fits_dir, "*", "rtiger_result.rds"))
if (length(rds) == 0) stop("No rtiger_result.rds under ", fits_dir)
cat(sprintf("Converting %d donor fit(s) from %s\n", length(rds), fits_dir))

calls <- purrr::map(rds, function(f) {
  donor <- basename(dirname(f))
  res   <- readRDS(f)
  seg   <- tryCatch(segments_from_rtiger_object(res),
                    error = function(e) { warning(donor, ": ", conditionMessage(e)); NULL })
  if (is.null(seg) || nrow(seg) == 0) return(NULL)
  dplyr::mutate(seg, source = "RTIGER_SNP50K", donor = donor, .before = 1)
}) |> purrr::list_rbind() |>
  dplyr::arrange(donor, name, chr, start_bp)

readr::write_csv(calls, out_csv)

# ---- summary + polarity sanity check --------------------------------------
# BC2S3 expectation: mean donor dosage ~ 0.125 (state contributes state/2 per bp).
# A flipped REF/ALT polarity would show ~0.875 instead.
per_sample <- calls |>
  dplyr::mutate(len = end_bp - start_bp + 1) |>
  dplyr::group_by(name) |>
  dplyr::summarise(dosage = sum(len * state / 2) / sum(len), .groups = "drop")

cat(sprintf("\nWrote %s\n", out_csv))
cat(sprintf("  donors: %d | samples: %d | segments: %d\n",
            dplyr::n_distinct(calls$donor), dplyr::n_distinct(calls$name), nrow(calls)))
cat(sprintf("  state mix (segment bp): AA=%.3f  AB=%.3f  BB=%.3f\n",
            { l <- calls$end_bp - calls$start_bp + 1
              sum(l[calls$state==0])/sum(l) }, { l <- calls$end_bp - calls$start_bp + 1
              sum(l[calls$state==1])/sum(l) }, { l <- calls$end_bp - calls$start_bp + 1
              sum(l[calls$state==2])/sum(l) }))
cat(sprintf("  mean donor dosage = %.3f  (BC2S3 expectation ~0.125; ~0.875 would mean flipped polarity)\n",
            mean(per_sample$dosage)))
