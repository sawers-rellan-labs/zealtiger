#!/usr/bin/env Rscript
# Generate a per-MEGABASE wideseq benchmark from simulated BC2S3 NILs.
#
# wideseq (method B; BzeaSeq docs/call_ancestry.Rmd, rpubs 1306822) is a 3-state
# REF/HET/ALT (0/1/2) caller: it bins allele counts to 1 Mb, computes the binned ALT
# frequency, clusters non-zero bins (k=3) and smooths with a 3-state HMM (BC2S3 priors).
# Its real input is the ~27.6M teosinte-vs-B73 variants from Schnable 2023 (~12,950/Mb,
# B73 NAM v5), genotyped at ~0.59x — so each 1 Mb bin aggregates THOUSANDS of informative
# variants and ALT_FREQ is well resolved despite low per-variant coverage.
#
# We simulate the ~2,300 per-Mb bins per plant DIRECTLY — never materialising 27M
# per-variant rows (that is the hundreds-of-GB server intermediate we want to avoid).
# Per bin we draw INFORMATIVE_VARIANT_COUNT, DEPTH_SUM and ALT_COUNT analytically from
# the bin's true ALT fraction + the WIDESEQ exp-floor coverage model:
#   present_prob = (1-pi)(1-e^{-k*lambda});  cond_mean = lambda/present_prob
#   INFORMATIVE_VARIANT_COUNT ~ Poisson(sites_per_mb * present_prob)
#   DEPTH_SUM = IVC + Poisson(IVC*(cond_mean-1))           # present site depth >= 1
#   ALT_COUNT ~ Binomial(DEPTH_SUM, p_eff(true_alt_frac))  # bin-level alt read prob
# The true per-bin ALT fraction comes from a coarse simcross sub-grid (points_per_bin),
# which also gives the bin's dominant true state. Same simcross mosaics/seed as
# make_rtiger_benchmark.R so RTIGER and wideseq share one ground truth.
#
# Output (results/wideseq_benchmark/):
#   bins/<name>_bin_genotypes.tsv : SAMPLE CONTIG BIN_POS INFORMATIVE_VARIANT_COUNT
#                                   DEPTH_SUM ALT_COUNT ALT_FREQ   (call_ancestry input)
#   truth_bins.csv                : name CONTIG BIN_POS true_state(0/1/2) true_alt_frac
#   truth_segments.csv            : name chr start_bp end_bp state  (RLE of true_state)
#   sample_stats.csv, params.json
#
# Usage:
#   Rscript make_wideseq_benchmark.R [n_lines] [sites_per_mb] [lambda] [points_per_bin]
#     n_lines        : default 100
#     sites_per_mb   : teosinte-panel variant density (default 12948 = 1294.82/100kb,
#                      the Schnable2023 27.6M-SNP density)
#     lambda         : mean per-variant coverage (default 0.59 = wideseq)
#     points_per_bin : simcross truth sub-grid points per Mb (default 40 -> 25kb breakpts)
#
# Requires results/maize_map_v5_clean.rds (run run_all.R Stage 1 first).

suppressMessages({
  library(tidyverse)
  library(simcross)
})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args <- commandArgs(trailingOnly = TRUE)
params <- list(
  n_lines        = if (length(args) >= 1) as.integer(args[1]) else 100L,
  sites_per_mb   = if (length(args) >= 2) as.numeric(args[2]) else 12948,
  lambda_mean    = if (length(args) >= 3) as.numeric(args[3]) else 0.590,
  points_per_bin = if (length(args) >= 4) as.integer(args[4]) else 40L,
  bin_size       = 1e6,
  # WIDESEQ exp-floor missing-data model (NOT SNP50K): missing(l)=pi+(1-pi)e^{-k l}
  pi_floor       = 0.346,
  k_decay        = 1.256,
  # Per-sample lambda ~ NORMAL, fit to 1434 BZea wideseq samples (agent/fit_wideseq_lambda.R):
  # mean 0.590, sd 0.226 (CV 0.38, skew 0.34) -- Normal beats Gamma/lognormal by AIC.
  lambda_sd      = 0.226,
  lambda_min     = 0.05,    # positivity floor (real population has a small near-0 failed-lib tail)
  error          = 0.005,
  m              = 10, p = 0,
  donor_allele   = 2L,
  seed           = 20260609L  # SAME seed as make_rtiger_benchmark.R -> shared mosaics
)
set.seed(params$seed)

stopifnot(fs::file_exists("results/maize_map_v5_clean.rds"))
anchors_clean <- readRDS("results/maize_map_v5_clean.rds")
chr_cm_lengths <- anchors_clean %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(L = max(cm), .groups = "drop") %>% dplyr::arrange(chr)

genome_mb <- anchors_clean %>% group_by(chr) %>%
  summarise(span = max(bp) - min(bp), .groups = "drop") %>%
  summarise(mb = sum(span) / 1e6) %>% pull(mb)

# Coarse TRUTH sub-grid only (cheap): points_per_bin per Mb, NOT the 27M variants.
n_truth <- as.integer(round(genome_mb * params$points_per_bin))
markers <- build_marker_grid(anchors_clean, n_truth, chr_prefix = "chr")
markers <- dplyr::mutate(markers, BIN_POS = as.integer(ceiling(bp / params$bin_size)))

out_dir <- "results/wideseq_benchmark"
bin_dir <- file.path(out_dir, "bins")
fs::dir_create(bin_dir)
cat(sprintf("wideseq benchmark: %d lines, ~%.0f Mb, sites/Mb=%.0f, lambda=%.2f, truth subgrid=%d/Mb\n",
            params$n_lines, genome_mb, params$sites_per_mb, params$lambda_mean, params$points_per_bin))

# per-sample coverage ~ Normal(mean, sd) fit to the BZea wideseq population, floored
# at lambda_min for positivity (sd=0 -> constant lambda).
draw_lambda <- function(n) {
  if (params$lambda_sd <= 0) return(rep(params$lambda_mean, n))
  pmax(stats::rnorm(n, params$lambda_mean, params$lambda_sd), params$lambda_min)
}
lambdas <- draw_lambda(params$n_lines)

ped <- bc2s3_pedigree()
stopifnot(check_pedigree(ped, ignore_sex = TRUE))

dominant_state <- function(d) (which.max(tabulate(d + 1L, nbins = 3L)) - 1L)

names_vec  <- sprintf("sim_%04d", seq_len(params$n_lines))
truth_bin  <- vector("list", params$n_lines)
truth_seg  <- vector("list", params$n_lines)
stats_list <- vector("list", params$n_lines)

for (i in seq_len(params$n_lines)) {
  dosage <- simulate_sample_dosage(ped, chr_cm_lengths, markers,
                                   m = params$m, p = params$p,
                                   donor_allele = params$donor_allele)

  # per-bin TRUE state + alt fraction from the coarse sub-grid
  tb <- tibble::tibble(CONTIG = markers$chr_label, BIN_POS = markers$BIN_POS, dosage) %>%
    dplyr::group_by(CONTIG, BIN_POS) %>%
    dplyr::summarise(true_state = dominant_state(dosage),
                     true_alt_frac = mean(dosage) / 2, .groups = "drop") %>%
    dplyr::mutate(CONTIG = factor(CONTIG, levels = paste0("chr", 1:10))) %>%
    dplyr::arrange(CONTIG, BIN_POS)
  nb <- nrow(tb)

  # ---- DIRECT per-bin observation draw (no per-variant intermediate) -------
  # lambda = DEPTH_SUM/VARIANT_COUNT = mean reads per wideseq site (incl. uncovered);
  # present_prob = 1 - missing; a COVERED site carries cond_mean = lambda/present ~1.7
  # reads (low coverage, like the skim) -> the per-bin ALT_FREQ is well resolved only
  # because thousands of sites are pooled, NOT because any site is deep.
  lam          <- lambdas[i]
  present_prob <- (1 - params$pi_floor) * (1 - exp(-params$k_decay * lam))
  cond_mean    <- lam / present_prob
  variant_count <- stats::rpois(nb, params$sites_per_mb)              # panel sites in bin
  ivc          <- stats::rbinom(nb, variant_count, present_prob)      # covered (informative) sites
  depth_sum    <- ivc + stats::rpois(nb, ivc * max(cond_mean - 1, 0)) # covered-site depth >= 1
  p_eff        <- tb$true_alt_frac * (1 - params$error) + (1 - tb$true_alt_frac) * params$error
  alt_count    <- stats::rbinom(nb, depth_sum, p_eff)
  alt_freq     <- ifelse(depth_sum > 0, alt_count / depth_sum, 0)

  obs <- tibble::tibble(
    SAMPLE = names_vec[i], CONTIG = tb$CONTIG, BIN_POS = tb$BIN_POS,
    VARIANT_COUNT = variant_count, INFORMATIVE_VARIANT_COUNT = ivc, DEPTH_SUM = depth_sum,
    ALT_COUNT = alt_count, ALT_FREQ = alt_freq)
  readr::write_tsv(obs, file.path(bin_dir, paste0(names_vec[i], "_bin_genotypes.tsv")))

  truth_bin[[i]] <- tb %>% dplyr::transmute(name = names_vec[i], CONTIG, BIN_POS,
                                            true_state, true_alt_frac)
  truth_seg[[i]] <- tb %>%
    dplyr::group_by(CONTIG) %>%
    dplyr::mutate(run = cumsum(c(TRUE, diff(true_state) != 0))) %>%
    dplyr::group_by(CONTIG, run, true_state) %>%
    dplyr::summarise(start_bp = (min(BIN_POS) - 1L) * params$bin_size,
                     end_bp   = max(BIN_POS) * params$bin_size, .groups = "drop") %>%
    dplyr::transmute(name = names_vec[i], chr = as.integer(sub("chr", "", CONTIG)),
                     start_bp, end_bp, state = true_state) %>%
    dplyr::arrange(chr, start_bp)

  n_co <- sum(unlist(tapply(tb$true_state, tb$CONTIG, function(s) sum(diff(s) != 0))))
  stats_list[[i]] <- tibble::tibble(
    name = names_vec[i], lambda = lam,
    lambda_realized  = sum(depth_sum) / sum(variant_count),   # == DEPTH_SUM/VARIANT_COUNT
    missing_realized = 1 - sum(ivc) / sum(variant_count),     # == 1 - INFORMATIVE/VARIANT
    mean_depth_sum = mean(depth_sum), mean_ivc = mean(ivc),
    n_bins = nb, n_co_detectable = n_co, donor_frac = mean(tb$true_alt_frac))
}

readr::write_csv(purrr::list_rbind(truth_bin), file.path(out_dir, "truth_bins.csv"))
readr::write_csv(purrr::list_rbind(truth_seg), file.path(out_dir, "truth_segments.csv"))
sample_stats <- purrr::list_rbind(stats_list)
readr::write_csv(sample_stats, file.path(out_dir, "sample_stats.csv"))

jsonlite::write_json(
  c(params, list(
    genome_mb = genome_mb, n_truth_subgrid = nrow(markers),
    realized_mean_depth_sum = mean(sample_stats$mean_depth_sum),
    realized_mean_ivc       = mean(sample_stats$mean_ivc),
    mean_co_detectable      = mean(sample_stats$n_co_detectable),
    simcross = as.character(packageVersion("simcross")),
    date = as.character(Sys.time())
  )),
  file.path(out_dir, "params.json"), pretty = TRUE, auto_unbox = TRUE)

cat("\n================ WIDESEQ BENCHMARK DATASET ================\n")
cat(sprintf("Lines: %d | %d bins/line | files in %s\n",
            params$n_lines, round(mean(sample_stats$n_bins)), bin_dir))
cat(sprintf("Realized: lambda (DEPTH_SUM/VARIANT_COUNT) %.3f, missingness %.3f (targets 0.59 / 0.67)\n",
            mean(sample_stats$lambda_realized), mean(sample_stats$missing_realized)))
cat(sprintf("Realized per-bin: mean DEPTH_SUM %.0f, mean informative variants %.0f (~%.1f reads/covered site)\n",
            mean(sample_stats$mean_depth_sum), mean(sample_stats$mean_ivc),
            mean(sample_stats$mean_depth_sum) / mean(sample_stats$mean_ivc)))
cat(sprintf("Mean detectable COs/line (1Mb res): %.1f | mean donor frac: %.4f\n",
            mean(sample_stats$n_co_detectable), mean(sample_stats$donor_frac)))
cat("Wrote bins/, truth_bins.csv, truth_segments.csv, sample_stats.csv, params.json\n")
