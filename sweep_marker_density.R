#!/usr/bin/env Rscript
# Marker-density sweep for the single-locus validation: Cohen's d and z of the
# Mb-fraction deviation from the BC2S3 expectation, vs marker count. Controlled
# design (same simulated NILs re-measured on EVERY grid, so the trend isolates
# marker discretization, not run-to-run noise), replicated R times to bound the
# estimates (SE of d across reps). See docs/13-single-locus-validation.md.
#
# Effect-size focus: Cohen's d = (mean - expected)/SD is N-independent in
# expectation; z = (mean - expected)/SE = d*sqrt(N) is reported too, but with
# large N even a negligible d is "significant" -- judge by d.
#
# Usage:
#   Rscript sweep_marker_density.R [--n=1500] [--reps=5] [--m=10] [--p=0]
#       [--grids=1000,2000,4000,8000,20000,50000,100000,200000]
#       [--map=results/maize_map_v5_clean.rds] [--seed=100]
#       [--out=results/single_locus_check]

suppressMessages({ library(dplyr); library(readr) })
getopt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
N      <- as.integer(getopt("--n", "1500"))
reps   <- as.integer(getopt("--reps", "5"))
m_int  <- as.integer(getopt("--m", "10"))
p_esc  <- as.numeric(getopt("--p", "0"))
grids  <- as.integer(strsplit(getopt("--grids",
            "1000,2000,4000,8000,20000,50000,100000,200000"), ",")[[1]])
map_rds<- getopt("--map", "results/maize_map_v5_clean.rds")
base   <- as.integer(getopt("--seed", "100"))
outdir <- getopt("--out", "results/single_locus_check")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source("R/02_simulate.R"); source("R/05_make_rtiger_input.R"); source("R/07_single_locus.R")
anchors <- readRDS(map_rds)
ped <- bc2s3_pedigree()
ccl <- anchors %>% group_by(chr) %>% summarise(L = max(cm), .groups = "drop") %>% arrange(chr)
p0  <- single_locus_p0()
mk  <- lapply(grids, function(M) build_marker_grid(anchors, M))
genome_mb <- sum(tapply(anchors$bp, anchors$chr, function(b) diff(range(b)))) / 1e6

dose_on <- function(sim, m) {
  d <- integer(nrow(m))
  for (ci in seq_along(ccl$L)) {
    ch <- ccl$chr[ci]; r <- which(m$chr == ch); if (!length(r)) next
    nil <- sim[[ci]][["8"]]; cm <- m$cm[r]
    d[r] <- (hap_allele_at(nil$mat, cm) == 2L) + (hap_allele_at(nil$pat, cm) == 2L)
  }
  d
}
cohen_d <- function(x, mu) (mean(x) - mu) / stats::sd(x)
zval    <- function(x, mu) (mean(x) - mu) / (stats::sd(x) / sqrt(length(x)))

rows <- list()
for (rep in seq_len(reps)) {
  set.seed(base + rep)
  Fs <- lapply(grids, function(x) matrix(NA_real_, N, 3))
  for (i in seq_len(N)) {
    sim <- simcross::sim_from_pedigree(ped, L = ccl$L, m = m_int, p = p_esc)
    for (g in seq_along(grids)) Fs[[g]][i, ] <- state_fractions_mb(mk[[g]], dose_on(sim, mk[[g]]))
  }
  for (g in seq_along(grids)) {
    Fg <- Fs[[g]]
    rows[[length(rows) + 1]] <- data.frame(
      rep = rep, markers = grids[g], spacing_Mb = genome_mb / grids[g],
      het_est = mean(Fg[, 2]),
      het_d = cohen_d(Fg[, 2], p0["Het"]), het_z = zval(Fg[, 2], p0["Het"]),
      alt_d = cohen_d(Fg[, 3], p0["ALT"]), alt_z = zval(Fg[, 3], p0["ALT"]),
      ref_d = cohen_d(Fg[, 1], p0["REF"]), ref_z = zval(Fg[, 1], p0["REF"]))
  }
  cat(sprintf("rep %d/%d done\n", rep, reps))
}
per_rep <- bind_rows(rows)
write_csv(per_rep, file.path(outdir, "marker_sweep_per_rep.csv"))

summ <- per_rep %>% group_by(markers, spacing_Mb) %>%
  summarise(across(c(het_est, het_d, het_z, alt_d, ref_d),
                   list(mean = mean, se = ~ sd(.x) / sqrt(reps))),
            .groups = "drop") %>% arrange(markers)
write_csv(summ, file.path(outdir, "marker_sweep_summary.csv"))

cat(sprintf("\nN=%d, reps=%d, genome=%.0f Mb, expected Het=%.5f\n\n", N, reps, genome_mb, p0["Het"]))
show <- summ %>% transmute(markers, spacing_Mb = round(spacing_Mb, 2),
  het_est = round(het_est_mean, 4),
  het_d = sprintf("%+.3f±%.3f", het_d_mean, het_d_se),
  het_z = sprintf("%+.2f", het_z_mean),
  alt_d = sprintf("%+.3f", alt_d_mean), ref_d = sprintf("%+.3f", ref_d_mean))
print(as.data.frame(show), row.names = FALSE)
cat(sprintf("\nWrote: %s/marker_sweep_{per_rep,summary}.csv\n", outdir))
