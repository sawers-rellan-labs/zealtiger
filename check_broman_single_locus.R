#!/usr/bin/env Rscript
# Validate the Broman/simcross BC2S3 simulation against the single-locus
# breeding-scheme expectation, using Mb-based (bp-weighted) genome fractions.
# Thin driver: all reusable logic is in R/07_single_locus.R. See
# docs/13-single-locus-validation.md for the full rationale (why this is a
# circularity/correctness check, why Hotelling not Pearson, why Mb on real
# marker coordinates).
#
# Usage:
#   Rscript check_broman_single_locus.R [--n=1500] [--m=10] [--p=0]
#       [--markers=20000] [--markers-file=PATH] [--map=results/maize_map_v5_clean.rds]
#       [--seed=1] [--out=results/single_locus_check]
#
# --markers-file: a real panel (TSV) to evaluate on its own genomic coordinates;
#   auto-detects chr/pos columns (chr|chr_v5|CONTIG, bp|pos|pos_v5|POSITION).
#   When omitted, a synthetic bp-even grid of --markers points is used.

suppressMessages({ library(dplyr); library(tibble); library(readr) })

getopt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), commandArgs(TRUE), value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
N        <- as.integer(getopt("--n", "1500"))
m_int    <- as.integer(getopt("--m", "10"))
p_esc    <- as.numeric(getopt("--p", "0"))
n_mark   <- as.integer(getopt("--markers", "20000"))
mfile    <- getopt("--markers-file", "")
map_rds  <- getopt("--map", "results/maize_map_v5_clean.rds")
seed     <- as.integer(getopt("--seed", "1"))
outdir   <- getopt("--out", "results/single_locus_check")
set.seed(seed)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source("R/02_simulate.R")          # bc2s3_pedigree
source("R/05_make_rtiger_input.R") # build_marker_grid, make_bp_to_cm, simulate_sample_dosage, truth_segments_from_dosage
source("R/07_single_locus.R")      # single_locus_p0, state_fractions_mb, hotelling_fractions

stopifnot(file.exists(map_rds))
anchors <- readRDS(map_rds)
stopifnot(all(c("chr", "cm", "bp") %in% names(anchors)))
ped <- bc2s3_pedigree()
chr_cm_lengths <- anchors %>% group_by(chr) %>%
  summarise(L = max(cm), .groups = "drop") %>% arrange(chr)

# Load a real panel on its own coordinates, or build a synthetic bp-even grid.
load_panel <- function(path, anchors) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  cc <- intersect(c("chr", "chr_v5", "CONTIG", "contig"), names(df))[1]
  pc <- intersect(c("bp", "pos", "pos_v5", "POSITION", "position"), names(df))[1]
  stopifnot(!is.na(cc), !is.na(pc))
  m <- tibble(chr = as.integer(gsub("[^0-9]", "", as.character(df[[cc]]))),
              bp  = as.integer(df[[pc]])) %>%
    filter(!is.na(chr), !is.na(bp), chr %in% anchors$chr) %>% arrange(chr, bp)
  ab <- split(anchors, anchors$chr); m$cm <- NA_real_
  for (ch in unique(m$chr)) {
    r <- which(m$chr == ch)
    m$cm[r] <- make_bp_to_cm(ab[[as.character(ch)]])(m$bp[r])
  }
  filter(m, !is.na(cm))
}

if (nzchar(mfile)) {
  markers <- load_panel(mfile, anchors); panel <- basename(mfile)
} else {
  markers <- build_marker_grid(anchors, n_markers = n_mark); panel <- sprintf("synthetic(%d)", n_mark)
}
cat(sprintf("Map: %s | panel: %s | %d markers | m=%d p=%g | N=%d\n",
            basename(map_rds), panel, nrow(markers), m_int, p_esc, N))

p0 <- single_locus_p0(n_bc = 2L, n_self = 3L)
cat(sprintf("Single-locus p0: REF=%.6f Het=%.6f ALT=%.6f (dosage=%.4f)\n",
            p0["REF"], p0["Het"], p0["ALT"], p0["ALT"] + 0.5 * p0["Het"]))

cat(sprintf("Simulating %d NILs ...\n", N))
F <- t(vapply(seq_len(N), function(i)
  state_fractions_mb(markers, simulate_sample_dosage(ped, chr_cm_lengths, markers,
                                                     m = m_int, p = p_esc)),
  numeric(3)))
colnames(F) <- c("REF", "Het", "ALT")
write_csv(as_tibble(F) %>% mutate(nil_id = row_number(), .before = 1),
          file.path(outdir, "fractions.csv"))

res <- hotelling_fractions(F, p0)
cat("\n== Per-state Mb fractions (mean over NILs) ==\n")
print(res$table, row.names = FALSE, digits = 5)
cat(sprintf("\n== Hotelling/Wald (Het,ALT), df=2 ==\n  T^2=%.3f chi2 p=%.4g [F=%.3f p=%.4g]\n",
            res$T2, res$p_chisq, res$F_stat, res$p_F))
cat(sprintf("  N_eff(ALT) = p(1-p)/sd_fi^2 = %.1f  (vs %d markers)\n",
            with(res$table[res$table$state == "ALT", ], expected * (1 - expected) / sd_fi^2),
            nrow(markers)))
cat(sprintf("\n%s\n", if (all(res$table$in_ci))
  "PASS: all expectations inside Wald CIs (sampler is Mendelian)." else
  "CHECK: an expectation outside its CI -- inspect pedigree/donor/union logic."))
cat(" (Large N can give small joint p even when correct; judge by CIs + effect size.)\n")
write_csv(res$table, file.path(outdir, "summary.csv"))
cat(sprintf("\nWrote: %s/{fractions,summary}.csv\n", outdir))
