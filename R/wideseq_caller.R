#!/usr/bin/env Rscript
# Wideseq (method B) 3-state Kgmm_HMM ancestry caller, factored from BzeaSeq
# docs/call_ancestry.Rmd so score_wideseq.R and the joint comparison share one
# implementation. Operates on a per-sample 1Mb bin table (CONTIG, BIN_POS, ALT_FREQ):
# ALT_FREQ==0 -> REF; non-zero bins clustered k=3 by a Gaussian mixture (rebmix; kmeans
# fallback) and relabelled REF/HET/ALT by cluster mean; 3-state HMM smoothing per
# chromosome with BC2S3 Mendelian start priors, fixed emissions, high self-transition.

.ws_has_rebmix <- function() requireNamespace("rebmix", quietly = TRUE)

# Evaluate expr while swallowing stdout + stderr chatter (e.g. REBMIX banners).
ws_quiet <- function(expr) {
  val <- NULL
  utils::capture.output(utils::capture.output(val <- expr, type = "message"))
  val
}

# Exact 3-state Viterbi (== HMM::viterbi for the same params); obs_idx in 1:3.
ws_viterbi3 <- function(obs_idx, start, trans, emiss) {
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

# BC2S3 Mendelian start priors -> c(REF, HET, ALT); donor = aa, recurrent = AA.
calculate_nil_frequencies <- function(bc = 2, s = 3, crossing_pop = c(AA = 0, Aa = 1, aa = 0)) {
  backcross_AA <- matrix(c(1, 1/2, 0,  0, 1/2, 1,  0, 0, 0), nrow = 3, byrow = TRUE)
  selfing      <- matrix(c(1, 1/4, 0,  0, 1/2, 0,  0, 1/4, 1), nrow = 3, byrow = TRUE)
  pop <- matrix(crossing_pop, ncol = 1)
  for (i in seq_len(bc)) pop <- backcross_AA %*% pop
  for (i in seq_len(s))  pop <- selfing %*% pop
  f <- as.vector(pop); names(f) <- c("AA", "Aa", "aa")
  c(REF = f[["AA"]], HET = f[["Aa"]], ALT = f[["aa"]])
}

# Map non-zero-bin clusters to REF/HET/ALT by ascending mean ALT_FREQ.
relabel_clusters <- function(clusters, alt_freq) {
  cluster_means <- tapply(alt_freq, clusters, mean)
  ordered <- order(cluster_means)
  lev <- c("REF", "HET", "ALT"); map <- rep(NA_character_, length(unique(clusters)))
  for (k in seq_along(ordered)) map[ordered[k]] <- lev[min(k, 3)]
  factor(map[as.integer(factor(clusters))], levels = lev)
}

gmm_clusters <- function(x, K = 3) {
  if (length(unique(round(x, 4))) < K) return(as.integer(factor(rank(x, ties.method = "min"))))
  if (.ws_has_rebmix()) {
    cl <- tryCatch(ws_quiet({          # REBMIX/RCLRMIX print version banners; swallow
      est <- rebmix::REBMIX(Dataset = list(data.frame(Value = x)),
                            Preprocessing = "histogram", cmin = K, cmax = K,
                            Criterion = "BIC", pdf = "normal")
      as.integer(rebmix::RCLRMIX(x = est)@Zp)
    }), error = function(e) NULL)
    if (!is.null(cl) && length(unique(cl)) >= 2) return(cl)
  }
  stats::kmeans(x, centers = min(K, length(unique(x))))$cluster   # fallback
}

# 3-state HMM smoothing of a REF/HET/ALT call vector; `start` = BC2S3 priors.
smooth_hmm <- function(genotypes, start, transitions = c(0.995, 0.005)) {
  geno_num <- as.numeric(factor(genotypes, levels = c("REF", "HET", "ALT"))) - 1
  trans <- matrix(c(transitions[1], transitions[2]/2, transitions[2]/2,
                    transitions[2]/2, transitions[1], transitions[2]/2,
                    transitions[2]/2, transitions[2]/2, transitions[1]), nrow = 3, byrow = TRUE)
  emiss <- matrix(c(0.9, 0.08, 0.02,  0.1, 0.8, 0.1,  0.02, 0.08, 0.9), nrow = 3, byrow = TRUE)
  path <- ws_viterbi3(geno_num + 1L, start = start, trans = trans, emiss = emiss)
  factor(c("REF", "HET", "ALT")[path], levels = c("REF", "HET", "ALT"))
}

# Call one sample's bin table -> adds $call (pre-HMM) and $call_hmm. `start` = BC2S3 priors.
ws_call_sample <- function(df, start) {
  df <- df |>
    dplyr::mutate(CONTIG = factor(CONTIG, levels = paste0("chr", 1:10))) |>
    dplyr::arrange(CONTIG, BIN_POS)
  df$call <- factor("REF", levels = c("REF", "HET", "ALT"))
  nz <- which(df$ALT_FREQ > 0)
  if (length(nz) >= 3) {
    df$call[nz] <- relabel_clusters(gmm_clusters(df$ALT_FREQ[nz], 3), df$ALT_FREQ[nz])
  } else if (length(nz) > 0) {
    df$call[nz] <- ifelse(df$ALT_FREQ[nz] >= 0.75, "ALT", "HET")
  }
  df |> dplyr::group_by(CONTIG) |>
    dplyr::mutate(call_hmm = smooth_hmm(call, start = start)) |> dplyr::ungroup()
}
