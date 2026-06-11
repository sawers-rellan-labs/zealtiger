#!/usr/bin/env Rscript
# Generate a per-MEGABASE wideseq benchmark from simulated BC2S3 NILs.
#
# wideseq (method B; docs/call_ancestry.Rmd in the BzeaSeq repo, rpubs 1306822) is a
# 3-state REF/HET/ALT (0/1/2) caller: it bins allele counts to 1 Mb, computes the
# binned ALT frequency, clusters non-zero bins (k=3: Ckmeans / arcsin / GMM) and
# smooths with a 3-state HMM (BC2S3 Mendelian priors). So the natural benchmark is
# per-Mb, not a 2.5 Gb genome at 0.4x.
#
# We reuse the SAME simcross mosaics and the SAME SNP50K coverage model as
# make_rtiger_benchmark.R (build_marker_grid + simulate_sample_dosage +
# draw_allele_counts), then aggregate per 1 Mb EXACTLY as BzeaSeq's
# ancestry_processing.R does:
#   BIN_POS = ceiling(POSITION / bin_size)
#   INFORMATIVE_VARIANT_COUNT = sum(REF+ALT > 0)
#   DEPTH_SUM = sum(REF+ALT);  ALT_FREQ = sum(ALT)/DEPTH_SUM
# Because both benchmarks derive from the same per-variant draws, RTIGER (per marker)
# and wideseq (per Mb) can be scored against ONE shared ground truth.
#
# Emits, under results/wideseq_benchmark/:
#   bins/<name>_bin_genotypes.tsv : SAMPLE CONTIG BIN_POS INFORMATIVE_VARIANT_COUNT
#                                   DEPTH_SUM ALT_COUNT ALT_FREQ   (the call_ancestry input)
#   truth_bins.csv                : name CONTIG BIN_POS true_state(0/1/2) true_alt_frac n_var
#   truth_segments.csv            : name chr start_bp end_bp state  (RLE of true_state)
#   sample_stats.csv, params.json
#
# Usage:
#   Rscript make_wideseq_benchmark.R [n_lines] [variants_per_mb] [lambda]
#     n_lines         : default 100
#     variants_per_mb : informative wideseq sites per Mb (default 200; set to match the
#                       real panel density — count wideseq_all.pos / genome Mb)
#     lambda          : mean per-variant coverage (default 0.43 = SNP50K skim; raise for
#                       higher-coverage wideseq WGS)
#
# Requires results/maize_map_v5_clean.rds (run run_all.R Stage 1 first).

suppressMessages({
  library(tidyverse)
  library(simcross)
})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args <- commandArgs(trailingOnly = TRUE)
params <- list(
  n_lines         = if (length(args) >= 1) as.integer(args[1]) else 100L,
  variants_per_mb = if (length(args) >= 2) as.integer(args[2]) else 200L,
  lambda_mean     = if (length(args) >= 3) as.numeric(args[3]) else 0.43,
  bin_size        = 1e6,
  lambda_shape    = 2.5,      # Gamma shape for per-sample lambda (Inf = constant)
  lambda_min      = 0.15,
  pi_floor        = 0.161,    # SNP50K exp-floor structural missingness
  k_decay         = 1.042,    # SNP50K exp-floor coverage decay
  error           = 0.005,
  m               = 10, p = 0,
  donor_allele    = 2L,
  seed            = 20260609L # SAME seed as make_rtiger_benchmark.R -> shared mosaics
)
set.seed(params$seed)

stopifnot(fs::file_exists("results/maize_map_v5_clean.rds"))
anchors_clean <- readRDS("results/maize_map_v5_clean.rds")
chr_cm_lengths <- anchors_clean %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(L = max(cm), .groups = "drop") %>%
  dplyr::arrange(chr)

# Genome size in Mb -> total variant grid (cheap: ~genome_Mb * variants_per_mb markers,
# far below per-base simulation). build_marker_grid spreads markers proportional to span.
genome_mb  <- anchors_clean %>% group_by(chr) %>%
  summarise(span = max(bp) - min(bp), .groups = "drop") %>%
  summarise(mb = sum(span) / 1e6) %>% pull(mb)
n_markers  <- as.integer(round(genome_mb * params$variants_per_mb))

out_dir <- "results/wideseq_benchmark"
bin_dir <- file.path(out_dir, "bins")
fs::dir_create(bin_dir)

markers <- build_marker_grid(anchors_clean, n_markers, chr_prefix = "chr")
markers <- dplyr::mutate(markers, BIN_POS = as.integer(ceiling(bp / params$bin_size)))
cat(sprintf("wideseq benchmark: %d variant sites (~%d/Mb over %.0f Mb), %d lines, lambda=%.2f\n",
            nrow(markers), params$variants_per_mb, genome_mb, params$n_lines, params$lambda_mean))

# per-sample mean coverage: Gamma(mean = lambda_mean), pinned, floored (as in the skim)
draw_lambda <- function(n) {
  if (!is.finite(params$lambda_shape)) return(rep(params$lambda_mean, n))
  stats::rgamma(n, shape = params$lambda_shape,
                rate = params$lambda_shape / params$lambda_mean)
}
lambdas <- draw_lambda(params$n_lines)
lambdas <- lambdas * (params$lambda_mean / mean(lambdas))
lambdas <- pmax(lambdas, params$lambda_min)

ped <- bc2s3_pedigree()
stopifnot(check_pedigree(ped, ignore_sex = TRUE))

names_vec  <- sprintf("sim_%04d", seq_len(params$n_lines))
truth_bin  <- vector("list", params$n_lines)
truth_seg  <- vector("list", params$n_lines)
stats_list <- vector("list", params$n_lines)

# dominant 0/1/2 state in a bin (ties -> lower state, conservative toward REF)
dominant_state <- function(d) {
  tab <- tabulate(d + 1L, nbins = 3L)   # counts of states 0,1,2
  (which.max(tab) - 1L)
}

for (i in seq_len(params$n_lines)) {
  dosage <- simulate_sample_dosage(ped, chr_cm_lengths, markers,
                                   m = params$m, p = params$p,
                                   donor_allele = params$donor_allele)
  cts   <- draw_allele_counts(dosage, lambdas[i], params$pi_floor,
                              params$k_decay, params$error)
  depth <- cts$ref + cts$alt

  # ---- bin to 1 Mb, exactly as ancestry_processing.R ----------------------
  binned <- tibble::tibble(
    CONTIG = markers$chr_label, BIN_POS = markers$BIN_POS,
    ref = cts$ref, alt = cts$alt, depth = depth, dosage = dosage
  ) %>%
    dplyr::group_by(CONTIG, BIN_POS) %>%
    dplyr::summarise(
      INFORMATIVE_VARIANT_COUNT = sum(depth > 0),
      DEPTH_SUM                 = sum(depth),
      ALT_COUNT                 = sum(alt),
      ALT_FREQ                  = ifelse(sum(depth) > 0, sum(alt) / sum(depth), 0),
      true_state                = dominant_state(dosage),
      true_alt_frac             = mean(dosage) / 2,
      n_var                     = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(CONTIG = factor(CONTIG, levels = paste0("chr", 1:10))) %>%
    dplyr::arrange(CONTIG, BIN_POS)

  # observed input file the wideseq caller consumes
  obs <- binned %>%
    dplyr::transmute(SAMPLE = names_vec[i], CONTIG, BIN_POS,
                     INFORMATIVE_VARIANT_COUNT, DEPTH_SUM, ALT_COUNT, ALT_FREQ)
  readr::write_tsv(obs, file.path(bin_dir, paste0(names_vec[i], "_bin_genotypes.tsv")))

  truth_bin[[i]] <- binned %>%
    dplyr::transmute(name = names_vec[i], CONTIG, BIN_POS, true_state, true_alt_frac, n_var)

  # truth segments = RLE of the per-bin true_state along each chromosome
  truth_seg[[i]] <- binned %>%
    dplyr::group_by(CONTIG) %>%
    dplyr::mutate(run = cumsum(c(TRUE, diff(true_state) != 0))) %>%
    dplyr::group_by(CONTIG, run, true_state) %>%
    dplyr::summarise(start_bp = (min(BIN_POS) - 1L) * params$bin_size,
                     end_bp   = max(BIN_POS) * params$bin_size,
                     .groups = "drop") %>%
    dplyr::transmute(name = names_vec[i], chr = as.integer(sub("chr", "", CONTIG)),
                     start_bp, end_bp, state = true_state) %>%
    dplyr::arrange(chr, start_bp)

  n_co <- sum(unlist(tapply(binned$true_state, binned$CONTIG,
                            function(s) sum(diff(s) != 0))))
  stats_list[[i]] <- tibble::tibble(
    name = names_vec[i], lambda = lambdas[i],
    mean_depth_sum = mean(binned$DEPTH_SUM),
    mean_ivc       = mean(binned$INFORMATIVE_VARIANT_COUNT),
    n_bins         = nrow(binned),
    n_co_detectable = n_co,
    donor_frac     = mean(binned$true_alt_frac)
  )
}

readr::write_csv(purrr::list_rbind(truth_bin), file.path(out_dir, "truth_bins.csv"))
readr::write_csv(purrr::list_rbind(truth_seg), file.path(out_dir, "truth_segments.csv"))
sample_stats <- purrr::list_rbind(stats_list)
readr::write_csv(sample_stats, file.path(out_dir, "sample_stats.csv"))

jsonlite::write_json(
  c(params, list(
    n_markers = nrow(markers), genome_mb = genome_mb,
    realized_mean_depth_sum = mean(sample_stats$mean_depth_sum),
    realized_mean_ivc       = mean(sample_stats$mean_ivc),
    mean_co_detectable      = mean(sample_stats$n_co_detectable),
    simcross = as.character(packageVersion("simcross")),
    date = as.character(Sys.time())
  )),
  file.path(out_dir, "params.json"), pretty = TRUE, auto_unbox = TRUE
)

cat("\n================ WIDESEQ BENCHMARK DATASET ================\n")
cat(sprintf("Lines: %d | %d bins/line (~%d/Mb variants) | files in %s\n",
            params$n_lines, round(mean(sample_stats$n_bins)), params$variants_per_mb, bin_dir))
cat(sprintf("Realized per-bin: mean DEPTH_SUM %.1f, mean informative variants %.1f\n",
            mean(sample_stats$mean_depth_sum), mean(sample_stats$mean_ivc)))
cat(sprintf("Mean detectable COs/line (1Mb res): %.1f | mean donor frac: %.4f\n",
            mean(sample_stats$n_co_detectable), mean(sample_stats$donor_frac)))
cat("Wrote bins/<name>_bin_genotypes.tsv, truth_bins.csv, truth_segments.csv,",
    "sample_stats.csv, params.json\n")
cat("NOTE: tune variants_per_mb / lambda so DEPTH_SUM & INFORMATIVE_VARIANT_COUNT match",
    "a real <sample>_bin_genotypes.tsv before trusting the comparison.\n")
