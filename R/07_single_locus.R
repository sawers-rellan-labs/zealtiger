# Stage 7 — single-locus validation of the simulated population
#
# The BC2S3 simulation's single-locus MARGINAL genotype distribution equals the
# breeding-scheme transition-matrix expectation, exactly and independent of the
# map and of interference m (a marginal quantity; see proof in
# docs/13-single-locus-validation.md). These functions compute that expectation,
# tabulate per-NIL Mb-based REF/Het/ALT fractions on a marker panel, and test the
# population means against it with the design-correct Hotelling chi-square
# (sample = unit; EMPIRICAL between-sample covariance, NOT per-marker binomial —
# linkage overdispersion makes a naive Pearson/binomial reject a correct sampler).

#' Single-locus genotype expectation for a BC(n_bc) S(n_self) scheme
#'
#' Iterates the backcross-to-recurrent and selfing transition matrices on
#' genotype states ordered (AA = donor hom, Aa = het, aa = recurrent hom),
#' starting from the F1 (Aa). Map- and interference-independent.
#'
#' @param n_bc Number of backcrosses to the recurrent parent (BC2 -> 2).
#' @param n_self Number of selfing generations (S3 -> 3).
#'
#' @return Named numeric `c(REF, Het, ALT)` summing to 1.
#' @export
single_locus_p0 <- function(n_bc = 2L, n_self = 3L) {
  B <- rbind(c(0, 1, 0), c(0, 0.5, 0.5), c(0, 0, 1))      # backcross to aa
  S <- rbind(c(1, 0, 0), c(0.25, 0.5, 0.25), c(0, 0, 1))  # selfing
  v <- c(0, 1, 0)                                          # F1 = Aa; order (AA,Aa,aa)
  for (i in seq_len(n_bc))   v <- as.numeric(v %*% B)
  for (i in seq_len(n_self)) v <- as.numeric(v %*% S)
  c(REF = v[3], Het = v[2], ALT = v[1])
}

#' Per-NIL Mb-based REF/Het/ALT fractions from a dosage vector
#'
#' Mb is read off the marker bp positions (run-length truth segments), so the
#' fraction is bp-weighted on the panel's own coordinates. Requires
#' `truth_segments_from_dosage()` (R/05).
#'
#' @param markers Marker grid (chr, bp).
#' @param dosage Integer 0/1/2 aligned to `markers`.
#'
#' @return Named numeric `c(REF, Het, ALT)` Mb fractions summing to 1.
#' @export
state_fractions_mb <- function(markers, dosage) {
  seg <- truth_segments_from_dosage(markers, dosage)
  seg$mb <- (seg$end_bp - seg$start_bp) / 1e6
  tot <- tapply(seg$mb, factor(seg$state, levels = 0:2), sum)
  tot[is.na(tot)] <- 0
  f <- as.numeric(tot) / sum(tot)
  c(REF = f[1], Het = f[2], ALT = f[3])
}

#' Design-correct goodness-of-fit of population fractions to an expectation
#'
#' Sample = unit. Hotelling/Wald chi-square on the (Het, ALT) mean vector using
#' the EMPIRICAL between-sample covariance (df = 2); per-state Wald 95% CI with
#' `SE = sd(f_i)/sqrt(N)`. NOT a per-marker Pearson/binomial test — linkage makes
#' those use a too-small variance and an enormous fake N, rejecting a correct
#' sampler every time. `sd_fi` is the ergodicity-limited between-sample spread
#' (N_eff = p(1-p)/sd_fi^2 ~ #recombination blocks, not #markers).
#'
#' @param F N x 3 matrix of per-NIL fractions, columns REF/Het/ALT.
#' @param p0 Named expectation `c(REF, Het, ALT)`.
#' @param comp Two components for the joint test (REF is redundant: sums to 1).
#'
#' @return list: T2, df, p_chisq, F_stat, p_F, table (per-state CIs), N.
#' @export
hotelling_fractions <- function(F, p0, comp = c("Het", "ALT")) {
  stopifnot(ncol(F) == 3, all(c("REF", "Het", "ALT") %in% colnames(F)))
  N <- nrow(F); phat <- colMeans(F)
  d  <- phat[comp] - p0[comp]
  S2 <- stats::cov(F[, comp])
  T2 <- as.numeric(N * t(d) %*% solve(S2) %*% d)
  Ff <- (N - 2) / (2 * (N - 1)) * T2
  cv <- stats::cov(F); se <- sqrt(diag(cv) / N); z <- stats::qnorm(0.975)
  st <- c("REF", "Het", "ALT")
  tab <- data.frame(
    state = st, expected = as.numeric(p0[st]), estimate = as.numeric(phat[st]),
    se = as.numeric(se[st]),
    ci_lo = as.numeric(phat[st] - z * se[st]),
    ci_hi = as.numeric(phat[st] + z * se[st]),
    sd_fi = as.numeric(sqrt(diag(cv))[st]))
  tab$in_ci <- tab$expected >= tab$ci_lo & tab$expected <= tab$ci_hi
  list(T2 = T2, df = 2L, p_chisq = stats::pchisq(T2, 2, lower.tail = FALSE),
       F_stat = Ff, p_F = stats::pf(Ff, 2, N - 2, lower.tail = FALSE),
       table = tab, N = N)
}

#' Two-panel forest plot of the fraction validation
#'
#' Fractions span ~2 orders of magnitude (REF ~0.86 vs Het ~0.03), so a single
#' linear forest plot can show magnitude or intervals but not both. This splits
#' them: (A) absolute fractions on a log10 axis (magnitude + expected marker),
#' (B) standardized deviation (estimate - expected)/SE with a shaded +/-1.96 band
#' (the scale-free pass criterion). Shared y (states).
#'
#' @param tab The `table` element from `hotelling_fractions()`.
#' @param title Optional plot title.
#'
#' @return A patchwork object (A | B).
#' @export
plot_fraction_forest <- function(tab, title = NULL) {
  ord <- c("REF", "ALT", "Het")                     # high -> low, top -> bottom
  tab$state <- factor(tab$state, levels = rev(ord))
  tab$z <- (tab$estimate - tab$expected) / tab$se
  tab$pass <- ifelse(tab$in_ci, "in CI", "out")
  pal <- c("in CI" = "#2c7fb8", "out" = "#d7301f")

  pA <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$state)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$ci_lo, xend = .data$ci_hi,
                                       yend = .data$state, colour = .data$pass),
                          linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$estimate, colour = .data$pass), size = 2.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$expected), shape = 124, size = 6,
                        colour = "black") +   # expected = vertical tick
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::labs(x = "Mb fraction (log scale) — point=estimate ±95% CI, tick=expected",
                  y = NULL, subtitle = "A. Magnitude") +
    ggplot2::theme_bw(base_size = 11)

  pB <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$state)) +
    ggplot2::annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -Inf, ymax = Inf,
                      fill = "grey85", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$z, yend = .data$state,
                                       colour = .data$pass), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$z, colour = .data$pass), size = 2.6) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::labs(x = "standardized deviation (estimate − expected)/SE",
                  y = NULL, subtitle = "B. Inside 95% CI?  (band = ±1.96)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())

  p <- patchwork::wrap_plots(pA, pB, widths = c(1, 1))
  if (!is.null(title)) p <- p + patchwork::plot_annotation(title = title)
  p
}
