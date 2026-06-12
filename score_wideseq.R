#!/usr/bin/env Rscript
# Score the wideseq (method B) 3-state caller against the per-Mb simulated truth.
#
# Runs the recommended call_ancestry.Rmd path (Kgmm_HMM): on each
# <name>_bin_genotypes.tsv, bins with ALT_FREQ==0 are REF by rule; non-zero bins are
# clustered k=3 by a Gaussian mixture (rebmix; kmeans fallback) on ALT_FREQ and
# relabelled REF/HET/ALT by cluster mean; then a 3-state HMM (HMM pkg, BC2S3 Mendelian
# priors, fixed emissions, high self-transition) smooths each chromosome.
#
# Scores per-bin calls vs results/wideseq_benchmark/truth_bins.csv using the SAME
# R/06 functions as score_rtiger.R (assign_called_state / score_markers /
# score_segments / score_cos) so wideseq and RTIGER outputs are directly comparable
# on the shared simcross truth.
#
# Usage:
#   Rscript score_wideseq.R [--bins DIR]
#     --bins : dir of <name>_bin_genotypes.tsv (default results/wideseq_benchmark/bins)
# Writes results/wideseq_benchmark/scoring/{marker,segment,crossover}_scores.csv.
#
# Caller is sourced from R/wideseq_caller.R (shared with the joint comparison); no hard
# deps beyond tidyverse/data.table (Viterbi inlined there), rebmix used if installed.

suppressMessages({library(tidyverse); library(data.table)})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args     <- commandArgs(trailingOnly = TRUE)
bins_dir <- if ("--bins" %in% args) args[which(args == "--bins") + 1] else
            "results/wideseq_benchmark/bins"
out_root <- "results/wideseq_benchmark"
out_dir  <- file.path(out_root, "scoring")
fs::dir_create(out_dir)
bc2s3 <- calculate_nil_frequencies(bc = 2, s = 3)   # BC2S3 HMM start priors (sourced)

# ---- assemble call sets across all samples ---------------------------------
truth_bins <- readr::read_csv(file.path(out_root, "truth_bins.csv"), show_col_types = FALSE)
truth_seg  <- readr::read_csv(file.path(out_root, "truth_segments.csv"), show_col_types = FALSE)
BIN <- 1e6
geno2num <- c(REF = 0L, HET = 1L, ALT = 2L)

bin_files <- fs::dir_ls(bins_dir, glob = "*_bin_genotypes.tsv")
stopifnot(length(bin_files) > 0)
cat(sprintf("Scoring %d samples (caller=%s)\n", length(bin_files),
            if (.ws_has_rebmix()) "Kgmm_HMM (rebmix)" else "Kkmeans_HMM (rebmix missing)"))

called_list <- vector("list", length(bin_files))
for (i in seq_along(bin_files)) {
  obs  <- readr::read_tsv(bin_files[i], show_col_types = FALSE)
  nm   <- obs$SAMPLE[1]
  cs   <- ws_call_sample(obs, bc2s3)
  called_list[[i]] <- cs %>%
    transmute(name = nm, chr = as.integer(sub("chr", "", CONTIG)), BIN_POS,
              state = geno2num[as.character(call_hmm)])
}
called_bins <- purrr::list_rbind(called_list)

# per-bin "markers" (bin midpoints) carrying the TRUE state, for score_markers()
truth_markers <- truth_bins %>%
  transmute(name, chr = as.integer(sub("chr", "", CONTIG)),
            bp = as.integer((BIN_POS - 0.5) * BIN), dosage = true_state)

# called segments = RLE of per-bin calls per (name, chr)
called_seg <- called_bins %>% arrange(name, chr, BIN_POS) %>%
  group_by(name, chr) %>%
  mutate(run = cumsum(c(TRUE, diff(state) != 0))) %>%
  group_by(name, chr, run) %>%
  summarise(start_bp = (min(BIN_POS) - 1L) * BIN, end_bp = max(BIN_POS) * BIN,
            state = first(state), .groups = "drop") %>%
  select(name, chr, start_bp, end_bp, state)
readr::write_csv(called_seg, file.path(out_dir, "called_segments.csv"))  # for size-dist comparison

# ---- score (same R/06 functions as score_rtiger.R) -------------------------
mk <- assign_called_state(called_seg, truth_markers)
sm <- score_markers(mk)

seg_classes <- list(
  homozygous_donor = 2L, heterozygous = 1L, donor_present = c(1L, 2L))
seg_tbl <- purrr::imap_dfr(seg_classes, function(tgt, nm)
  score_segments(called_seg, truth_seg, target = tgt, label = nm)$summary)

co <- score_cos(called_seg, truth_seg) %>%
  summarise(mean_true_co = mean(true_co), mean_called_co = mean(called_co),
            co_bias = mean(called_co - true_co),
            co_rmse = sqrt(mean((called_co - true_co)^2)))

marker_tbl <- sm$per_state %>% mutate(accuracy = sm$accuracy, uncovered = sm$uncovered)
readr::write_csv(marker_tbl, file.path(out_dir, "marker_scores.csv"))
readr::write_csv(seg_tbl,    file.path(out_dir, "segment_scores.csv"))
readr::write_csv(co,         file.path(out_dir, "crossover_scores.csv"))

cat("\n================ BIN-LEVEL (per-state recall/precision) ============\n")
print(as.data.frame(marker_tbl), digits = 3)
cat("\n================ SEGMENT-LEVEL (reciprocal-overlap >= 0.5) =========\n")
print(as.data.frame(seg_tbl), digits = 3)
cat("\n================ CROSSOVERS ========================================\n")
print(as.data.frame(co), digits = 3)
cat("\nConfusion (rows=true, cols=called):\n"); print(sm$confusion)
cat(sprintf("\nWrote marker_scores.csv, segment_scores.csv, crossover_scores.csv in %s\n", out_dir))
