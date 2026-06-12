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
# No hard package deps beyond tidyverse/data.table: Viterbi is inlined; rebmix is
# used for the GMM if installed, else stats::kmeans (base R) is the fallback.

suppressMessages({library(tidyverse); library(data.table)})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

# Exact 3-state Viterbi (== HMM::viterbi for the same params; avoids the HMM dep).
# obs_idx in 1:3 (observed REF/HET/ALT symbol); start/trans/emiss over REF/HET/ALT.
viterbi3 <- function(obs_idx, start, trans, emiss) {
  n <- length(obs_idx); K <- 3L
  if (n == 1L) return(which.max(log(start) + log(emiss[, obs_idx[1]])))
  ls <- log(start); lt <- log(trans); le <- log(emiss)
  delta <- matrix(-Inf, n, K); psi <- matrix(0L, n, K)
  delta[1, ] <- ls + le[, obs_idx[1]]
  for (t in 2:n) for (j in 1:K) {
    v <- delta[t - 1, ] + lt[, j]
    psi[t, j] <- which.max(v); delta[t, j] <- max(v) + le[j, obs_idx[t]]
  }
  path <- integer(n); path[n] <- which.max(delta[n, ])
  for (t in (n - 1):1) path[t] <- psi[t + 1, path[t + 1]]
  path
}

# Evaluate expr while swallowing both stdout and stderr chatter (e.g. REBMIX banners).
# Outer capture = stdout (also swallows the inner call's visible return); inner = stderr.
quiet <- function(expr) {
  val <- NULL
  utils::capture.output(
    utils::capture.output(val <- expr, type = "message"))
  val
}

args     <- commandArgs(trailingOnly = TRUE)
bins_dir <- if ("--bins" %in% args) args[which(args == "--bins") + 1] else
            "results/wideseq_benchmark/bins"
out_root <- "results/wideseq_benchmark"
out_dir  <- file.path(out_root, "scoring")
fs::dir_create(out_dir)
has_rebmix <- requireNamespace("rebmix", quietly = TRUE)

# ---- BC2S3 Mendelian priors (REF/HET/ALT), per call_ancestry.Rmd -----------
calculate_nil_frequencies <- function(bc = 2, s = 3, crossing_pop = c(AA = 0, Aa = 1, aa = 0)) {
  backcross_AA <- matrix(c(1, 1/2, 0,  0, 1/2, 1,  0, 0, 0), nrow = 3, byrow = TRUE)
  selfing      <- matrix(c(1, 1/4, 0,  0, 1/2, 0,  0, 1/4, 1), nrow = 3, byrow = TRUE)
  pop <- matrix(crossing_pop, ncol = 1)
  for (i in seq_len(bc)) pop <- backcross_AA %*% pop
  for (i in seq_len(s))  pop <- selfing %*% pop
  f <- as.vector(pop); names(f) <- c("AA", "Aa", "aa")
  # donor = aa (ALT); recurrent = AA (REF)
  c(REF = f[["AA"]], HET = f[["Aa"]], ALT = f[["aa"]])
}
bc2s3 <- calculate_nil_frequencies(bc = 2, s = 3)

# ---- cluster non-zero bins into REF/HET/ALT by ascending mean ALT_FREQ ------
relabel_clusters <- function(clusters, alt_freq) {
  cluster_means <- tapply(alt_freq, clusters, mean)
  ordered <- order(cluster_means)
  map <- rep(NA_character_, length(unique(clusters)))
  lev <- c("REF", "HET", "ALT")
  for (k in seq_along(ordered)) map[ordered[k]] <- lev[min(k, 3)]
  factor(map[as.integer(factor(clusters))], levels = lev)
}

gmm_clusters <- function(x, K = 3) {
  if (length(unique(round(x, 4))) < K) return(as.integer(factor(rank(x, ties.method = "min"))))
  if (has_rebmix) {
    cl <- tryCatch(quiet({          # REBMIX/RCLRMIX print version banners; swallow them
      est <- rebmix::REBMIX(Dataset = list(data.frame(Value = x)),
                            Preprocessing = "histogram", cmin = K, cmax = K,
                            Criterion = "BIC", pdf = "normal")
      as.integer(rebmix::RCLRMIX(x = est)@Zp)
    }), error = function(e) NULL)
    if (!is.null(cl) && length(unique(cl)) >= 2) return(cl)
  }
  stats::kmeans(x, centers = min(K, length(unique(x))))$cluster   # fallback
}

# ---- 3-state HMM smoothing per chromosome (call_ancestry smooth_ancestry) ---
smooth_hmm <- function(genotypes, transitions = c(0.995, 0.005)) {
  geno_num <- as.numeric(factor(genotypes, levels = c("REF", "HET", "ALT"))) - 1
  trans <- matrix(c(transitions[1], transitions[2]/2, transitions[2]/2,
                    transitions[2]/2, transitions[1], transitions[2]/2,
                    transitions[2]/2, transitions[2]/2, transitions[1]),
                  nrow = 3, byrow = TRUE)
  emiss <- matrix(c(0.9, 0.08, 0.02,
                    0.1, 0.8,  0.1,
                    0.02, 0.08, 0.9), nrow = 3, byrow = TRUE)
  path <- viterbi3(geno_num + 1L, start = bc2s3, trans = trans, emiss = emiss)
  factor(c("REF", "HET", "ALT")[path], levels = c("REF", "HET", "ALT"))
}

# ---- call one sample's bin table -> per-bin REF/HET/ALT ---------------------
call_sample <- function(df) {
  df <- df %>% mutate(CONTIG = factor(CONTIG, levels = paste0("chr", 1:10))) %>%
    arrange(CONTIG, BIN_POS)
  df$call <- factor("REF", levels = c("REF", "HET", "ALT"))
  nz <- which(df$ALT_FREQ > 0)
  if (length(nz) >= 3) {
    cl <- gmm_clusters(df$ALT_FREQ[nz], K = 3)
    df$call[nz] <- relabel_clusters(cl, df$ALT_FREQ[nz])
  } else if (length(nz) > 0) {
    df$call[nz] <- ifelse(df$ALT_FREQ[nz] >= 0.75, "ALT", "HET")
  }
  df %>% group_by(CONTIG) %>%
    mutate(call_hmm = smooth_hmm(call)) %>% ungroup()
}

# ---- assemble call sets across all samples ---------------------------------
truth_bins <- readr::read_csv(file.path(out_root, "truth_bins.csv"), show_col_types = FALSE)
truth_seg  <- readr::read_csv(file.path(out_root, "truth_segments.csv"), show_col_types = FALSE)
BIN <- 1e6
geno2num <- c(REF = 0L, HET = 1L, ALT = 2L)

bin_files <- fs::dir_ls(bins_dir, glob = "*_bin_genotypes.tsv")
stopifnot(length(bin_files) > 0)
cat(sprintf("Scoring %d samples (caller=%s)\n", length(bin_files),
            if (has_rebmix) "Kgmm_HMM (rebmix)" else "Kkmeans_HMM (rebmix missing)"))

called_list <- vector("list", length(bin_files))
for (i in seq_along(bin_files)) {
  obs  <- readr::read_tsv(bin_files[i], show_col_types = FALSE)
  nm   <- obs$SAMPLE[1]
  cs   <- call_sample(obs)
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
