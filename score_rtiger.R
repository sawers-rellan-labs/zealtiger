#!/usr/bin/env Rscript
# Score RTIGER calls against the simulated ground truth.
#
# Usage:
#   Rscript score_rtiger.R --rtiger <dir-with-CompleteBlock-beds | results.rds>
#   Rscript score_rtiger.R --demo            # self-test: perfect + degraded truth
#
# Reads truth from results/rtiger_benchmark/ (truth_segments.csv, truth_markers.rds).
# Writes scoring tables + plots to results/rtiger_benchmark/scoring/.

suppressMessages({
  library(tidyverse); library(data.table)
})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args      <- commandArgs(trailingOnly = TRUE)
demo      <- "--demo" %in% args
rtiger_in <- if ("--rtiger" %in% args) args[which(args == "--rtiger") + 1] else NULL
truth_dir <- "results/rtiger_benchmark"
out_dir   <- file.path(truth_dir, "scoring")
fs::dir_create(out_dir)

truth_seg     <- readr::read_csv(file.path(truth_dir, "truth_segments.csv"),
                                 show_col_types = FALSE)
truth_markers <- readRDS(file.path(truth_dir, "truth_markers.rds"))

# ---- assemble the call sets to score --------------------------------------
call_sets <- list()
if (demo || is.null(rtiger_in)) {
  if (is.null(rtiger_in)) cat("No --rtiger given; running --demo self-test.\n")
  call_sets[["perfect (truth)"]] <-
    dplyr::select(truth_seg, name, chr, start_bp, end_bp, state)
  call_sets[["degraded mock"]] <- degrade_truth(truth_seg)
} else {
  call_sets[["RTIGER"]] <- read_rtiger_segments(rtiger_in)
}

score_one <- function(called, tag) {
  mk  <- assign_called_state(called, truth_markers)
  sm  <- score_markers(mk)
  # zygosity classes first, then the donor-present TOTAL (states 1+2) as the
  # aggregate row (its n_truth = homozygous_donor + heterozygous).
  seg <- dplyr::bind_rows(
    score_segments(called, truth_seg, target = 2L,        label = "homozygous_donor")$summary,
    score_segments(called, truth_seg, target = 1L,        label = "heterozygous")$summary,
    score_segments(called, truth_seg, target = c(1L, 2L), label = "donor_present")$summary
  )
  co  <- score_cos(called, truth_seg)
  list(tag = tag, markers = sm, segments = seg, cos = co,
       confusion = sm$confusion)
}

results <- purrr::imap(call_sets, ~ score_one(.x, .y))

# ---- write + report --------------------------------------------------------
seg_tbl <- purrr::map_dfr(results, ~ dplyr::mutate(.x$segments, call_set = .x$tag,
                                                   .before = 1))
mk_tbl  <- purrr::map_dfr(results, ~ dplyr::mutate(.x$markers$per_state,
             call_set = .x$tag, accuracy = .x$markers$accuracy,
             uncovered = .x$markers$uncovered, .before = 1))
co_tbl  <- purrr::map_dfr(results, ~ dplyr::summarise(.x$cos,
             call_set = .x$tag,
             mean_true_co = mean(true_co), mean_called_co = mean(called_co),
             co_bias = mean(called_co - true_co),
             co_rmse = sqrt(mean((called_co - true_co)^2))))

readr::write_csv(seg_tbl, file.path(out_dir, "segment_scores.csv"))
readr::write_csv(mk_tbl,  file.path(out_dir, "marker_scores.csv"))
readr::write_csv(co_tbl,  file.path(out_dir, "crossover_scores.csv"))

# CO scatter (true vs called) for the last/primary call set
primary <- results[[length(results)]]
p_co <- ggplot(primary$cos, aes(true_co, called_co)) +
  geom_abline(slope = 1, linetype = 2, colour = "grey50") +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.5) +
  labs(x = "True detectable COs / line", y = "RTIGER-called COs / line",
       title = paste0("Crossover recovery — ", primary$tag)) + theme_bw()
ggsave(file.path(out_dir, "crossover_recovery.png"), p_co, width = 6, height = 5,
       dpi = 150)

cat("\n================ MARKER-LEVEL (per-state recall/precision) ============\n")
print(as.data.frame(mk_tbl), row.names = FALSE, digits = 3)
cat("\n================ SEGMENT-LEVEL (reciprocal-overlap >= 0.5) ============\n")
print(as.data.frame(seg_tbl), row.names = FALSE, digits = 3)
cat("\n================ CROSSOVERS ===========================================\n")
print(as.data.frame(co_tbl), row.names = FALSE, digits = 3)
cat("\nConfusion matrix (rows=true, cols=called) for '", primary$tag, "':\n", sep = "")
print(primary$confusion)
cat("\nWrote segment_scores.csv, marker_scores.csv, crossover_scores.csv,",
    "crossover_recovery.png in", out_dir, "\n")
