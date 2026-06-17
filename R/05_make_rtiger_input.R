# Stage 5 (benchmark) — generate RTIGER input from simulated NILs
#
# Produces ground-truth-labelled, RTIGER-format allele-count files so RTIGER can
# be benchmarked against a known introgression structure. Models the SNP50K
# regime: donor = teosinte = ALT, recurrent = B73 = REF; coverage-driven
# missingness via the exp-floor model  missing(lambda) = pi + (1-pi) e^{-k lambda}.
#
# RTIGER input file format (6 cols, tab-separated, NO header):
#   Chr  Position  RefBase  RefCount  AltBase  AltCount
# Missing/uncovered sites are kept as rows with RefCount = AltCount = 0.

#' Monotone bp -> cM interpolator for one chromosome
#'
#' Inverse of `make_cm_to_bp()`: maps a physical position to genetic position,
#' so a marker's bp can be placed on the simulated (cM) mosaic.
#'
#' Hyman-filtered monotone cubic spline, matching `make_cm_to_bp()` so one
#' Marey-map representation serves both directions. Anchors must already be
#' monotone (`enforce_monotone()`); `method = "hyman"` errors otherwise (a useful
#' guard). Strictly monotone (flat cM segments stay flat), unlike `monoH.FC`.
#' Queries outside the anchor range are clamped to the end values (the
#' `rule = 2` behaviour of the old `approxfun`). In cold regions the cM range per
#' anchor gap is tiny, so this is ~indistinguishable from linear here; the spline
#' matters in the inverse (cM -> bp) direction.
#'
#' @param chr_anchors Clean anchors for ONE chromosome (columns cm, bp).
#'
#' @return A function bp -> cM (monotone spline, clamped at the ends).
#' @export
make_bp_to_cm <- function(chr_anchors) {
  stopifnot(all(c("cm", "bp") %in% names(chr_anchors)))
  d <- chr_anchors %>%
    dplyr::group_by(.data$bp) %>%
    dplyr::summarise(cm = mean(.data$cm), .groups = "drop") %>%
    dplyr::arrange(.data$bp)
  stopifnot(nrow(d) >= 2)
  f   <- stats::splinefun(d$bp, d$cm, method = "hyman")
  rng <- range(d$bp)
  function(bp) f(pmin(pmax(bp, rng[1]), rng[2]))  # clamp at ends (rule = 2 equiv)
}

#' Build a physical marker grid with per-marker cM and ref/alt bases
#'
#' Allocates `n_markers` across chromosomes in proportion to physical span and
#' spaces them evenly in bp within each chromosome's anchored range. Each marker
#' gets a cM position (for mosaic lookup) and arbitrary REF/ALT bases (cosmetic;
#' RTIGER uses counts, not bases).
#'
#' @param anchors_clean Clean anchor set (columns chr, cm, bp).
#' @param n_markers Target marker count (e.g. 50000).
#' @param chr_prefix Chromosome label prefix to match RTIGER seqlengths names.
#'
#' @return Tibble: chr (int), chr_label, bp (int), cm, ref_base, alt_base,
#'   sorted by chr, bp.
#' @export
build_marker_grid <- function(anchors_clean, n_markers, chr_prefix = "chr") {
  stopifnot(all(c("chr", "cm", "bp") %in% names(anchors_clean)), n_markers > 0)

  spans <- anchors_clean %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(min_bp = min(.data$bp), max_bp = max(.data$bp),
                     span = max(.data$bp) - min(.data$bp), .groups = "drop") %>%
    dplyr::mutate(quota = pmax(2L, as.integer(round(n_markers * .data$span /
                                                       sum(.data$span)))))

  bases <- c("A", "C", "G", "T")
  anchors_by_chr <- split(anchors_clean, anchors_clean$chr)

  purrr::map(seq_len(nrow(spans)), function(i) {
    ch <- spans$chr[i]
    q  <- spans$quota[i]
    bp <- round(seq(spans$min_bp[i], spans$max_bp[i], length.out = q))
    bp2cm <- make_bp_to_cm(anchors_by_chr[[as.character(ch)]])
    ref <- sample(bases, q, replace = TRUE)
    alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1L), character(1))
    tibble::tibble(
      chr       = as.integer(ch),
      chr_label = paste0(chr_prefix, ch),
      bp        = as.integer(bp),
      cm        = bp2cm(bp),
      ref_base  = ref,
      alt_base  = unname(alt)
    )
  }) %>%
    purrr::list_rbind() %>%
    dplyr::arrange(.data$chr, .data$bp)
}

#' Founder allele on one haplotype at given cM positions (vectorised)
#'
#' @param hap Haplotype list with `$alleles` and `$locations` (right endpoints).
#' @param cm Numeric vector of cM query positions.
#'
#' @return Integer vector of founder labels at each cM.
#' @export
hap_allele_at <- function(hap, cm) {
  idx <- findInterval(cm, hap$locations, left.open = TRUE) + 1L
  idx[idx > length(hap$alleles)] <- length(hap$alleles)
  idx[idx < 1L] <- 1L
  hap$alleles[idx]
}

#' Donor dosage (0/1/2) at each marker for one simulated NIL
#'
#' @param ped Pedigree from `bc2s3_pedigree()`.
#' @param chr_cm_lengths Tibble chr, L (arranged by chr).
#' @param markers Marker grid from `build_marker_grid()` (needs chr, cm).
#' @param m,p Stahl interference parameters.
#' @param donor_allele Founder label for the donor (default 2).
#' @param nil_id Final-NIL individual id as character (default "8").
#'
#' @return Integer vector (length nrow(markers)) of donor dosage 0/1/2.
#' @export
simulate_sample_dosage <- function(ped, chr_cm_lengths, markers, m = 10, p = 0,
                                   donor_allele = 2L, nil_id = "8") {
  L <- chr_cm_lengths$L
  sim <- simcross::sim_from_pedigree(ped, L = L, m = m, p = p)
  dosage <- integer(nrow(markers))
  for (ci in seq_along(L)) {
    ch   <- chr_cm_lengths$chr[ci]
    rows <- which(markers$chr == ch)
    if (!length(rows)) next
    nil <- sim[[ci]][[nil_id]]
    cm  <- markers$cm[rows]
    dosage[rows] <- (hap_allele_at(nil$mat, cm) == donor_allele) +
                    (hap_allele_at(nil$pat, cm) == donor_allele)
  }
  dosage
}

#' Draw REF/ALT read counts under the SNP50K coverage + missingness model
#'
#' Matches BOTH moments of the target regime exactly, by construction:
#'   * mean depth over all markers = `lambda`;
#'   * missingness = pi + (1-pi) e^{-k lambda}  (the exp-floor fit).
#' Since lambda is the mean over ALL markers (including the uncovered ones),
#' present markers must carry a higher conditional mean depth
#' `lambda / present_prob` (~1.4 reads at lambda = 0.43, so a het site usually
#' shows a single read and looks homozygous — the regime that stresses RTIGER).
#' Present-site depth is drawn as 1 + Poisson(cond_mean - 1) so every present
#' marker has >= 1 read. ALT is the donor allele: p(ALT) is 0 / 0.5 / 1 for
#' dosage 0 / 1 / 2, perturbed by sequencing `error`.
#'
#' @param dosage Integer vector of donor dosage 0/1/2.
#' @param lambda Mean per-marker coverage for this sample.
#' @param pi_floor Structural missingness floor.
#' @param k_decay Coverage decay constant from the exp-floor fit.
#' @param error Per-base sequencing error rate.
#'
#' @return List with integer vectors `ref` and `alt` (read counts).
#' @export
draw_allele_counts <- function(dosage, lambda, pi_floor, k_decay, error) {
  n <- length(dosage)
  present_prob <- (1 - pi_floor) * (1 - exp(-k_decay * lambda))
  cond_mean    <- lambda / present_prob          # mean depth | present (>= 1)
  present <- stats::runif(n) < present_prob
  depth <- integer(n)
  depth[present] <- 1L + stats::rpois(sum(present), max(cond_mean - 1, 0))
  p_alt <- c(0, 0.5, 1)[dosage + 1L]             # donor allele = ALT
  p_eff <- p_alt * (1 - error) + (1 - p_alt) * error
  alt <- stats::rbinom(n, depth, p_eff)
  list(ref = as.integer(depth - alt), alt = as.integer(alt))
}

#' Write one sample's RTIGER allele-count file (6 cols, no header)
#'
#' @param markers Marker grid (chr_label, bp, ref_base, alt_base).
#' @param ref,alt Integer read-count vectors.
#' @param path Output file path.
#'
#' @return Invisibly, `path`.
#' @export
write_rtiger_sample <- function(markers, ref, alt, path) {
  out <- data.frame(
    chr = markers$chr_label, pos = markers$bp,
    ref_base = markers$ref_base, ref_count = ref,
    alt_base = markers$alt_base, alt_count = alt
  )
  readr::write_tsv(out, path, col_names = FALSE)
  invisible(path)
}

#' Ground-truth segments (runs of constant dosage state) for one sample
#'
#' Run-length-encodes the dosage along markers per chromosome, giving segments
#' directly comparable to RTIGER's state calls. Boundaries are reported at the
#' flanking marker bp (so resolution is limited by marker spacing, as for the
#' real data).
#'
#' @param markers Marker grid (chr, bp).
#' @param dosage Integer dosage vector aligned to `markers`.
#'
#' @return Tibble: chr, start_bp, end_bp, state (0/1/2), n_markers.
#' @export
truth_segments_from_dosage <- function(markers, dosage) {
  df <- markers %>%
    dplyr::transmute(.data$chr, .data$bp, state = dosage) %>%
    dplyr::arrange(.data$chr, .data$bp)
  df %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::mutate(run = cumsum(c(TRUE, diff(.data$state) != 0))) %>%
    dplyr::group_by(.data$chr, .data$run, .data$state) %>%
    dplyr::summarise(start_bp = min(.data$bp), end_bp = max(.data$bp),
                     n_markers = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(.data$chr, .data$start_bp) %>%
    dplyr::select("chr", "start_bp", "end_bp", "state", "n_markers")
}
