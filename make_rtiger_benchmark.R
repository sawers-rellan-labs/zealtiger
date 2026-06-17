#!/usr/bin/env Rscript
# Generate an RTIGER benchmark dataset from simulated BC2S3 NILs.
#
# Emits, under results/rtiger_benchmark/:
#   counts/<name>.tsv     one RTIGER allele-count file per simulated NIL
#   expDesign.csv         files + name table for RTIGER(expDesign = ...)
#   seqlengths.csv        named chromosome lengths (bp) for RTIGER seqlengths
#   truth_segments.csv    ground-truth dosage segments (for precision/recall)
#   truth_markers.rds     per-marker true dosage (for marker-level scoring)
#   sample_stats.csv      per-sample lambda, realized coverage/missingness, COs
#   params.json           run parameters
#
# Usage:
#   Rscript make_rtiger_benchmark.R [n_lines]   # default 100
#
# Requires results/maize_map_v5_clean.rds (run run_all.R Stage 1 first).

suppressMessages({
  library(tidyverse)
  library(simcross)
})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args <- commandArgs(trailingOnly = TRUE)
# Coverage/grid defaults model the SNP50K skim; override via env vars to benchmark
# other regimes (e.g. seqcapture capped to ~20x: LAMBDA_MEAN=20 PI_FLOOR=0.002
# OUT_DIR=results/rtiger_benchmark_20x). pi_floor is the structural missingness
# (SNP50K ~0.16; GBTS seqcapture ~0.002 — near-complete).
env_num <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }
env_chr <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

# ---- data-driven coverage calibration from the real skim -------------------
# lambda mean + Gamma shape are ESTIMATED from the real skim per-sample depth
# (QC'd to >= floor), not hardcoded. Cached (real_skim_lambda.csv) so we don't
# re-read every count file each run; computed from the count files if the cache
# is absent. Env LAMBDA_MEAN/LAMBDA_SHAPE override; the literals below are only a
# fallback when the real skim mount is absent (e.g. a clean clone).
lambda_floor    <- env_num("LAMBDA_MIN", 0.10)
real_lambda_csv <- "results/sim_calibration/real_skim_lambda.csv"
skim_design     <- "data/rtiger_50K/expDesign_all.csv"
real_lambda <- if (fs::file_exists(real_lambda_csv)) {
  readr::read_csv(real_lambda_csv, show_col_types = FALSE)$lambda
} else if (fs::file_exists(skim_design)) {
  ed <- readr::read_csv(skim_design, show_col_types = FALSE)
  lam <- vapply(ed$files, function(p)
    mean(Reduce(`+`, data.table::fread(p, header = FALSE, select = c(4, 6)))), numeric(1))
  fs::dir_create(dirname(real_lambda_csv))
  readr::write_csv(tibble::tibble(name = ed$name, lambda = lam), real_lambda_csv)
  lam
} else NULL
cov <- if (!is.null(real_lambda)) fit_lambda_gamma(real_lambda, floor = lambda_floor) else
  list(mean = 0.43, shape = 8, n = NA_integer_, source = "fallback literal (no real skim)")
cat(sprintf("Coverage calibration [%s]: lambda_mean %.3f, Gamma shape %.2f (n=%s)\n",
            cov$source, cov$mean, cov$shape, as.character(cov$n)))

params <- list(
  n_lines      = if (length(args) >= 1) as.integer(args[1]) else 100L,
  n_markers    = as.integer(env_num("N_MARKERS", 50000)),
  lambda_mean  = env_num("LAMBDA_MEAN", cov$mean),     # estimated from real skim
  lambda_shape = env_num("LAMBDA_SHAPE", cov$shape),   # estimated Gamma shape (~8)
  lambda_min   = lambda_floor,                          # low-coverage floor
  pi_floor     = env_num("PI_FLOOR", 0.161),      # exp-floor structural missingness
  k_decay      = env_num("K_DECAY", 1.042),       # exp-floor coverage decay
  error        = env_num("ERROR", 0.005),         # per-base sequencing error
  m            = 10, p = 0,
  donor_allele = 2L,
  seed         = as.integer(env_num("SEED", 20260609))
)
out_dir   <- env_chr("OUT_DIR", "results/rtiger_benchmark")
set.seed(params$seed)

stopifnot(fs::file_exists("results/maize_map_v5_clean.rds"))
anchors_clean <- readRDS("results/maize_map_v5_clean.rds")
chr_cm_lengths <- anchors_clean %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(L = max(cm), .groups = "drop") %>%
  dplyr::arrange(chr)

count_dir <- file.path(out_dir, "counts")
fs::dir_create(count_dir)

markers <- build_marker_grid(anchors_clean, params$n_markers, chr_prefix = "chr")
seqlengths <- markers %>%
  dplyr::group_by(chr_label) %>%
  dplyr::summarise(len = max(bp) + 1L, .groups = "drop")

# per-sample mean coverage: Gamma(mean = lambda_mean) unless shape is Inf
draw_lambda <- function(n) {
  if (!is.finite(params$lambda_shape)) return(rep(params$lambda_mean, n))
  stats::rgamma(n, shape = params$lambda_shape,
                rate = params$lambda_shape / params$lambda_mean)
}
lambdas <- draw_lambda(params$n_lines)
lambdas <- lambdas * (params$lambda_mean / mean(lambdas))  # pin mean to target
lambdas <- pmax(lambdas, params$lambda_min)                # realistic coverage floor

ped <- bc2s3_pedigree()
stopifnot(check_pedigree(ped, ignore_sex = TRUE))

names_vec  <- sprintf("sim_%04d", seq_len(params$n_lines))
truth_seg  <- vector("list", params$n_lines)
truth_mark <- vector("list", params$n_lines)
stats_list <- vector("list", params$n_lines)

for (i in seq_len(params$n_lines)) {
  dosage <- simulate_sample_dosage(ped, chr_cm_lengths, markers,
                                   m = params$m, p = params$p,
                                   donor_allele = params$donor_allele)
  cts <- draw_allele_counts(dosage, lambdas[i], params$pi_floor,
                            params$k_decay, params$error)
  depth <- cts$ref + cts$alt
  path  <- file.path(count_dir, paste0(names_vec[i], ".tsv"))
  write_rtiger_sample(markers, cts$ref, cts$alt, path)

  seg <- truth_segments_from_dosage(markers, dosage)
  truth_seg[[i]]  <- dplyr::mutate(seg, name = names_vec[i], .before = 1)
  truth_mark[[i]] <- tibble::tibble(name = names_vec[i], chr = markers$chr,
                                    bp = markers$bp, dosage = dosage)
  # detectable COs at this marker resolution = dosage-state transitions
  n_co <- markers %>% dplyr::mutate(state = dosage) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(co = sum(diff(state) != 0), .groups = "drop") %>%
    dplyr::summarise(n = sum(co)) %>% dplyr::pull(n)
  stats_list[[i]] <- tibble::tibble(
    name = names_vec[i], lambda = lambdas[i],
    mean_depth = mean(depth), missing_obs = mean(depth == 0),
    n_co_detectable = n_co,
    donor_markers = sum(dosage > 0), donor_frac = mean(dosage) / 2
  )
}

# ---- write design + truth + metadata --------------------------------------
expDesign <- data.frame(
  files = fs::path_abs(file.path(count_dir, paste0(names_vec, ".tsv"))),
  name  = names_vec
)
readr::write_csv(expDesign, file.path(out_dir, "expDesign.csv"))
readr::write_csv(seqlengths, file.path(out_dir, "seqlengths.csv"))

truth_segments <- purrr::list_rbind(truth_seg)
readr::write_csv(truth_segments, file.path(out_dir, "truth_segments.csv"))
saveRDS(purrr::list_rbind(truth_mark), file.path(out_dir, "truth_markers.rds"))

sample_stats <- purrr::list_rbind(stats_list)
readr::write_csv(sample_stats, file.path(out_dir, "sample_stats.csv"))

jsonlite::write_json(
  c(params, list(
    seqnames = seqlengths$chr_label,
    realized_mean_coverage   = mean(sample_stats$mean_depth),
    realized_mean_missingness = mean(sample_stats$missing_obs),
    mean_co_detectable        = mean(sample_stats$n_co_detectable),
    simcross = as.character(packageVersion("simcross")),
    date = as.character(Sys.time())
  )),
  file.path(out_dir, "params.json"), pretty = TRUE, auto_unbox = TRUE
)

cat("\n================ RTIGER BENCHMARK DATASET ================\n")
cat(sprintf("Lines: %d | markers/line: %d | files in %s\n",
            params$n_lines, nrow(markers), count_dir))
cat(sprintf("Target SNP50K:   lambda %.2f, missingness %.3f\n",
            params$lambda_mean, params$pi_floor +
              (1 - params$pi_floor) * exp(-params$k_decay * params$lambda_mean)))
cat(sprintf("Realized:        mean coverage %.3f, mean missingness %.3f\n",
            mean(sample_stats$mean_depth), mean(sample_stats$missing_obs)))
cat(sprintf("Mean detectable COs/line (at 50K res): %.1f | mean donor frac: %.4f\n",
            mean(sample_stats$n_co_detectable), mean(sample_stats$donor_frac)))
cat(sprintf("Truth: %d segments across %d lines\n",
            nrow(truth_segments), params$n_lines))
cat("\nWrote expDesign.csv, seqlengths.csv, truth_segments.csv,",
    "truth_markers.rds, sample_stats.csv, params.json\n")
