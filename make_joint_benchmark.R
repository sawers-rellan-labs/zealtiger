#!/usr/bin/env Rscript
# ONE genome simulation, assessed by BOTH callers — decouples the genome/segment
# simulation from the chromosome-painting method.
#
# Per sample we draw the BC2S3 mosaic ONCE (simcross), then derive from that single
# mosaic: (a) a CANONICAL high-resolution truth, (b) the RTIGER input (SNP50K skim:
# per-marker allele counts at ~50K markers, lambda~0.43), (c) the wideseq input
# (teosinte panel: per-1Mb-bin ALT_FREQ, lambda~0.59). Both callers are later scored
# against the SAME canonical truth, so segment recall and crossover counts share one
# denominator. (Coverage differs by platform — a property of the assay, not the genome.)
#
# Output (results/joint_benchmark/):
#   truth_markers.rds   canonical per-marker dosage on a FINE grid (name,chr,bp,dosage)
#   truth_segments.csv  canonical RLE segments (name,chr,start_bp,end_bp,state)
#   rtiger/{expDesign.csv,seqlengths.csv,counts/<name>.tsv}   RTIGER input (skim)
#   wideseq/bins/<name>_bin_genotypes.tsv                     wideseq input (1Mb)
#   sample_stats.csv, params.json
#
# Usage: Rscript make_joint_benchmark.R [n_lines]
# Requires results/maize_map_v5_clean.rds.

suppressMessages({library(tidyverse); library(simcross)})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

args <- commandArgs(trailingOnly = TRUE)
P <- list(
  n_lines        = if (length(args) >= 1) as.integer(args[1]) else 100L,
  truth_per_mb   = 200L,     # canonical truth grid (~5 kb); ~ exact CO count
  n_markers_skim = 50000L,   # RTIGER marker grid
  bin_size       = 1e6,
  sites_per_mb   = 12948,    # wideseq teosinte-panel density
  m = 10, p = 0, donor_allele = 2L, error = 0.005, seed = 20260609L,
  # RTIGER / SNP50K skim coverage
  skim_lambda = 0.43, skim_shape = 2.5, skim_min = 0.15, skim_pi = 0.161, skim_k = 1.042,
  # wideseq coverage (Normal per-sample lambda fit to BZea; wideseq exp-floor)
  ws_lambda = 0.590, ws_sd = 0.226, ws_min = 0.05, ws_pi = 0.346, ws_k = 1.256
)
set.seed(P$seed)

stopifnot(fs::file_exists("results/maize_map_v5_clean.rds"))
anchors <- readRDS("results/maize_map_v5_clean.rds")
chr_cm <- anchors %>% group_by(chr) %>% summarise(L = max(cm), .groups = "drop") %>% arrange(chr)
genome_mb <- anchors %>% group_by(chr) %>% summarise(s = max(bp) - min(bp), .groups = "drop") %>%
  summarise(mb = sum(s) / 1e6) %>% pull(mb)

# Build BOTH marker grids up front (deterministic given the seed).
markers_fine <- build_marker_grid(anchors, as.integer(genome_mb * P$truth_per_mb), "chr") %>%
  mutate(BIN_POS = as.integer(ceiling(bp / P$bin_size)))
markers_skim <- build_marker_grid(anchors, P$n_markers_skim, "chr")
seqlengths_skim <- markers_skim %>% group_by(chr_label) %>%
  summarise(len = max(bp) + 1L, .groups = "drop")

out <- "results/joint_benchmark"
fs::dir_create(file.path(out, "rtiger", "counts"))
fs::dir_create(file.path(out, "wideseq", "bins"))

# dosage 0/1/2 from ONE sim object at an arbitrary marker grid (chr, cm)
dosage_at <- function(sim, mk) {
  d <- integer(nrow(mk))
  for (ci in seq_along(chr_cm$L)) {
    ch <- chr_cm$chr[ci]; rows <- which(mk$chr == ch); if (!length(rows)) next
    nil <- sim[[ci]][["8"]]; cm <- mk$cm[rows]
    d[rows] <- (hap_allele_at(nil$mat, cm) == P$donor_allele) +
               (hap_allele_at(nil$pat, cm) == P$donor_allele)
  }
  d
}
draw_norm <- function(n, mean, sd, lo) pmax(stats::rnorm(n, mean, sd), lo)
draw_gam  <- function(n, mean, shape, lo) {
  g <- stats::rgamma(n, shape = shape, rate = shape / mean); pmax(g * (mean / mean(g)), lo)
}
dominant <- function(d) (which.max(tabulate(d + 1L, nbins = 3L)) - 1L)

ped <- bc2s3_pedigree(); stopifnot(check_pedigree(ped, ignore_sex = TRUE))
lam_skim <- draw_gam(P$n_lines, P$skim_lambda, P$skim_shape, P$skim_min)
lam_ws   <- draw_norm(P$n_lines, P$ws_lambda, P$ws_sd, P$ws_min)

nm <- sprintf("sim_%04d", seq_len(P$n_lines))
truth_mark <- truth_seg <- stats_list <- vector("list", P$n_lines)
cat(sprintf("Joint benchmark: %d lines | truth %d/Mb, RTIGER %d markers (skim), wideseq %.0f sites/Mb (1Mb)\n",
            P$n_lines, P$truth_per_mb, P$n_markers_skim, P$sites_per_mb))

for (i in seq_len(P$n_lines)) {
  sim <- simcross::sim_from_pedigree(ped, L = chr_cm$L, m = P$m, p = P$p)  # ONE mosaic
  d_fine <- dosage_at(sim, markers_fine)
  d_skim <- dosage_at(sim, markers_skim)

  # ---- canonical truth (fine) ----
  truth_mark[[i]] <- tibble(name = nm[i], chr = markers_fine$chr, bp = markers_fine$bp, dosage = d_fine)
  truth_seg[[i]]  <- truth_segments_from_dosage(markers_fine, d_fine) %>% mutate(name = nm[i], .before = 1)

  # ---- RTIGER input: SNP50K skim allele counts ----
  cts <- draw_allele_counts(d_skim, lam_skim[i], P$skim_pi, P$skim_k, P$error)
  write_rtiger_sample(markers_skim, cts$ref, cts$alt,
                      file.path(out, "rtiger", "counts", paste0(nm[i], ".tsv")))

  # ---- wideseq input: per-1Mb bins from the SAME mosaic (true alt frac from fine grid) ----
  tb <- tibble(BIN_POS = markers_fine$BIN_POS, CONTIG = markers_fine$chr_label, dosage = d_fine) %>%
    group_by(CONTIG, BIN_POS) %>%
    summarise(true_alt_frac = mean(dosage) / 2, .groups = "drop") %>%
    mutate(CONTIG = factor(CONTIG, levels = paste0("chr", 1:10))) %>% arrange(CONTIG, BIN_POS)
  nb <- nrow(tb)
  present <- (1 - P$ws_pi) * (1 - exp(-P$ws_k * lam_ws[i])); cond <- lam_ws[i] / present
  vc  <- stats::rpois(nb, P$sites_per_mb)
  ivc <- stats::rbinom(nb, vc, present)
  dep <- ivc + stats::rpois(nb, ivc * max(cond - 1, 0))
  peff <- tb$true_alt_frac * (1 - P$error) + (1 - tb$true_alt_frac) * P$error
  altc <- stats::rbinom(nb, dep, peff)
  tibble(SAMPLE = nm[i], CONTIG = tb$CONTIG, BIN_POS = tb$BIN_POS,
         VARIANT_COUNT = vc, INFORMATIVE_VARIANT_COUNT = ivc, DEPTH_SUM = dep,
         ALT_COUNT = altc, ALT_FREQ = ifelse(dep > 0, altc / dep, 0)) %>%
    write_tsv(file.path(out, "wideseq", "bins", paste0(nm[i], "_bin_genotypes.tsv")))

  stats_list[[i]] <- tibble(name = nm[i], lam_skim = lam_skim[i], lam_ws = lam_ws[i],
    true_co = sum(unlist(tapply(d_fine, markers_fine$chr, function(s) sum(diff(s) != 0)))),
    donor_frac = mean(d_fine) / 2)
}

truth_markers <- purrr::list_rbind(truth_mark)
saveRDS(truth_markers, file.path(out, "truth_markers.rds"))
write_csv(purrr::list_rbind(truth_seg), file.path(out, "truth_segments.csv"))
ed <- data.frame(files = fs::path_abs(file.path(out, "rtiger", "counts", paste0(nm, ".tsv"))), name = nm)
write_csv(ed, file.path(out, "rtiger", "expDesign.csv"))
write_csv(seqlengths_skim, file.path(out, "rtiger", "seqlengths.csv"))
ss <- purrr::list_rbind(stats_list); write_csv(ss, file.path(out, "sample_stats.csv"))
jsonlite::write_json(c(P, list(genome_mb = genome_mb, mean_true_co = mean(ss$true_co),
  mean_donor_frac = mean(ss$donor_frac), date = as.character(Sys.time()))),
  file.path(out, "params.json"), pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("\nDONE. canonical truth: %d markers/line (%d/Mb), mean true COs/line %.1f, donor frac %.3f\n",
            nrow(markers_fine), P$truth_per_mb, mean(ss$true_co), mean(ss$donor_frac)))
cat("  RTIGER input -> rtiger/ (skim) ; wideseq input -> wideseq/bins/ (1Mb)\n")
cat("Next: fit RTIGER on rtiger/expDesign.csv, run wideseq on wideseq/bins/, score both vs truth_segments.csv\n")
