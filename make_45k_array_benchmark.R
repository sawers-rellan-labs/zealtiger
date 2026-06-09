#!/usr/bin/env Rscript
# Generate a 45K-ARRAY RTIGER benchmark to SELECT the rigidity for the MolBreeding run.
#
# Unlike the 50K skim benchmark (coverage exp-floor, het-undercalling), this models
# the MolBreeding array regime:
#   * markers sit at the REAL lifted v5 positions RTIGER will use (the ~21k
#     informative sites, read from data/molbreeding_45k/counts/), so spacing is faithful;
#   * emission is CLEAN genotype → 0/1/2 pseudo-counts (RefRef=(D,0), het=(D/2,D/2),
#     AltAlt=(0,D)) with ~1.3% MCAR dropout (measured) + a small miscall rate — NO
#     coverage exp-floor, because an array has no per-sample coverage.
#
# Ground truth = simulated BC2S3 NIL dosage mosaics (same simcross machinery as the
# 50K benchmark). Output mirrors results/rtiger_benchmark/ layout.
#
# Usage: Rscript make_45k_array_benchmark.R [n_lines] [pi_miss] [depth] [error]

suppressMessages({library(tidyverse); library(simcross)})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args    <- commandArgs(trailingOnly = TRUE)
n_lines <- if (length(args) >= 1) as.integer(args[1]) else 100L
pi_miss <- if (length(args) >= 2) as.numeric(args[2]) else 0.013   # measured array missingness
D       <- if (length(args) >= 3) as.integer(args[3]) else 20L
error   <- if (length(args) >= 4) as.numeric(args[4]) else 0.002   # array miscall rate (small)
Dh      <- as.integer(round(D / 2))
set.seed(20260609L)

out_dir   <- "results/rtiger_45k_benchmark"
count_dir <- file.path(out_dir, "counts")
fs::dir_create(count_dir)

# --- marker grid = the real informative v5 sites the MolBreeding run uses ------
cf <- readr::read_tsv(list.files("data/molbreeding_45k/counts", full.names = TRUE)[1],
                      col_names = c("chr_label","bp","rb","rc","ab","ac"), show_col_types = FALSE)
anchors <- readRDS("results/maize_map_v5_clean.rds")
chr_cm_lengths <- anchors |> group_by(chr) |> summarise(L = max(cm), .groups = "drop") |> arrange(chr)
interp <- split(anchors, anchors$chr) |> map(make_bp_to_cm)

markers <- cf |>
  transmute(chr = as.integer(sub("chr", "", chr_label)), chr_label, bp) |>
  arrange(chr, bp) |>
  mutate(ref_base = "A", alt_base = "C", cm = NA_real_)
for (c in sort(unique(markers$chr))) {
  idx <- markers$chr == c
  markers$cm[idx] <- interp[[as.character(c)]](markers$bp[idx])
}
cat(sprintf("45K array benchmark: %d markers, %d lines, pi_miss=%.3f, D=%d, error=%.3f\n",
            nrow(markers), n_lines, pi_miss, D, error))

# --- clean-array emission: genotype → pseudo-counts + MCAR dropout ------------
draw_array_counts <- function(dosage, pi_miss, D, Dh, error) {
  n <- length(dosage)
  st <- dosage
  if (error > 0) { fl <- runif(n) < error; if (any(fl)) st[fl] <- sample(0:2, sum(fl), replace = TRUE) }
  ref <- ifelse(st == 0L, D, ifelse(st == 1L, Dh, 0L))
  alt <- ifelse(st == 2L, D, ifelse(st == 1L, Dh, 0L))
  miss <- runif(n) < pi_miss
  ref[miss] <- 0L; alt[miss] <- 0L
  list(ref = as.integer(ref), alt = as.integer(alt))
}

ped <- bc2s3_pedigree(); stopifnot(check_pedigree(ped, ignore_sex = TRUE))
names_vec <- sprintf("sim_%04d", seq_len(n_lines))
truth_seg <- vector("list", n_lines); truth_mark <- vector("list", n_lines)

for (i in seq_len(n_lines)) {
  dosage <- simulate_sample_dosage(ped, chr_cm_lengths, markers, m = 10, p = 0, donor_allele = 2L)
  cts <- draw_array_counts(dosage, pi_miss, D, Dh, error)
  write_rtiger_sample(markers, cts$ref, cts$alt, file.path(count_dir, paste0(names_vec[i], ".tsv")))
  seg <- truth_segments_from_dosage(markers, dosage)
  truth_seg[[i]]  <- mutate(seg, name = names_vec[i], .before = 1)
  truth_mark[[i]] <- tibble(name = names_vec[i], chr = markers$chr, bp = markers$bp, dosage = dosage)
}

expDesign <- data.frame(files = fs::path_abs(file.path(count_dir, paste0(names_vec, ".tsv"))),
                        name = names_vec)
readr::write_csv(expDesign, file.path(out_dir, "expDesign.csv"))
file.copy("data/molbreeding_45k/seqlengths.csv", file.path(out_dir, "seqlengths.csv"), overwrite = TRUE)
readr::write_csv(list_rbind(truth_seg), file.path(out_dir, "truth_segments.csv"))
saveRDS(list_rbind(truth_mark), file.path(out_dir, "truth_markers.rds"))
cat(sprintf("Wrote benchmark to %s (mean donor frac %.4f)\n", out_dir,
            mean(map_dbl(truth_mark, ~ mean(.x$dosage) / 2))))
