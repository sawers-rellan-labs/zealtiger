#!/usr/bin/env Rscript
# Score BOTH callers against ONE canonical truth (results/joint_benchmark), produced by
# make_joint_benchmark.R from a single set of BC2S3 mosaics. RTIGER (r=8) ran on the skim
# input; wideseq (Kgmm_HMM, shared caller R/wideseq_caller.R) runs here on the 1Mb input.
# Both sets of called SEGMENTS are scored vs the canonical truth_markers / truth_segments
# with the same R/06 functions, so segment recall and crossover counts share ONE denominator.
#
# Usage: Rscript score_joint.R     (needs results/joint_benchmark/rtiger/rtiger_out_r8 fitted)
# Writes results/joint_benchmark/scoring/{marker,segment,crossover}_scores.csv and the
# introgression-size table docs/rtiger_vs_wideseq_blocks.csv (for the comparison notebook).

suppressMessages({library(tidyverse); library(data.table)})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

root <- "results/joint_benchmark"
out  <- file.path(root, "scoring"); fs::dir_create(out)
BIN  <- 1e6
truth_markers <- readRDS(file.path(root, "truth_markers.rds"))
truth_seg     <- read_csv(file.path(root, "truth_segments.csv"), show_col_types = FALSE)
bc2s3 <- calculate_nil_frequencies(2, 3)

# --- RTIGER called segments (from the saved r=8 fit) ---
rtg_seg <- segments_from_rtiger_object(
  readRDS(file.path(root, "rtiger", "rtiger_out_r8", "rtiger_result.rds")))

# --- wideseq called segments (shared Kgmm_HMM caller on the 1Mb bins) ---
geno2num <- c(REF = 0L, HET = 1L, ALT = 2L)
bin_files <- fs::dir_ls(file.path(root, "wideseq", "bins"), glob = "*_bin_genotypes.tsv")
ws_seg <- purrr::map_dfr(bin_files, function(f) {
  obs <- read_tsv(f, show_col_types = FALSE)
  ws_call_sample(obs, bc2s3) |>
    transmute(name = obs$SAMPLE[1], chr = as.integer(sub("chr", "", CONTIG)), BIN_POS,
              state = geno2num[as.character(call_hmm)]) |>
    arrange(name, chr, BIN_POS) |> group_by(name, chr) |>
    mutate(run = cumsum(c(TRUE, diff(state) != 0))) |> group_by(name, chr, run) |>
    summarise(start_bp = (min(BIN_POS) - 1L) * BIN, end_bp = max(BIN_POS) * BIN,
              state = first(state), .groups = "drop") |>
    select(name, chr, start_bp, end_bp, state)
})

score_one <- function(seg, tag) {
  sm <- score_markers(assign_called_state(seg, truth_markers))
  classes <- list(homozygous_donor = 2L, heterozygous = 1L, donor_present = c(1L, 2L))
  list(
    marker = sm$per_state |> mutate(method = tag, accuracy = sm$accuracy, uncovered = sm$uncovered),
    seg = purrr::imap_dfr(classes, ~ score_segments(seg, truth_seg, target = .x, label = .y)$summary) |>
            mutate(method = tag, .before = 1),
    co = score_cos(seg, truth_seg) |>
            summarise(method = tag, mean_true_co = mean(true_co), mean_called_co = mean(called_co),
                      co_bias = mean(called_co - true_co), co_rmse = sqrt(mean((called_co - true_co)^2))))
}
R <- score_one(rtg_seg, "RTIGER (r=8)")
W <- score_one(ws_seg,  "wideseq (Kgmm_HMM)")
marker <- bind_rows(R$marker, W$marker)
seg    <- bind_rows(R$seg, W$seg)
co     <- bind_rows(R$co, W$co)
write_csv(marker, file.path(out, "marker_scores.csv"))
write_csv(seg,    file.path(out, "segment_scores.csv"))
write_csv(co,     file.path(out, "crossover_scores.csv"))

# --- introgression-size distribution: donor-present blocks vs the SAME truth ---
donor_blocks <- function(s) s |> mutate(chr = as.integer(chr)) |> arrange(name, chr, start_bp) |>
  group_by(name, chr) |>
  mutate(donor = as.integer(state >= 1), run = cumsum(c(TRUE, diff(donor) != 0))) |>
  group_by(name, chr, run) |>
  summarise(donor = first(donor), start_bp = min(start_bp), end_bp = max(end_bp), .groups = "drop") |>
  filter(donor == 1) |> transmute(name, chr, size_mb = (end_bp - start_bp) / 1e6)
blocks <- bind_rows(
  donor_blocks(truth_seg) |> mutate(src = "simulation (truth)"),
  donor_blocks(rtg_seg)   |> mutate(src = "RTIGER (r=8)"),
  donor_blocks(ws_seg)    |> mutate(src = "wideseq (Kgmm_HMM)"))
write_csv(blocks, "docs/rtiger_vs_wideseq_blocks.csv")

cat("================ MARKER-LEVEL (vs canonical fine truth) ================\n")
print(as.data.frame(marker), digits = 3)
cat("\n================ SEGMENT-LEVEL (reciprocal overlap >= 0.5) ============\n")
print(as.data.frame(seg[, c("method","class","n_truth","n_called","precision","recall","F1","FDR")]), digits = 3)
cat("\n================ CROSSOVERS (vs the SAME canonical truth) =============\n")
print(as.data.frame(co), digits = 3)
cat("\n================ INTROGRESSION-SIZE (donor-present blocks) ============\n")
blocks |> group_by(src) |>
  summarise(blocks_per_sample = round(n() / n_distinct(name), 1),
            median_Mb = round(median(size_mb), 1), mean_Mb = round(mean(size_mb), 1),
            pct_under_1Mb = round(100 * mean(size_mb < 1), 1),
            donor_Mb_per_sample = round(sum(size_mb) / n_distinct(name), 0), .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)
cat("\nWrote scoring CSVs + docs/rtiger_vs_wideseq_blocks.csv\n")
