# Stage 3 — cM -> Mb conversion
#
# Functions only. See spec section 6.

#' Build a monotone cM -> bp interpolator for one chromosome
#'
#' Collapses cM ties to mean bp (so the relation is a function), then returns a
#' Hyman-filtered monotone cubic spline. Anchors must already be monotone (see
#' `enforce_monotone()`) -- `method = "hyman"` errors on non-monotone input,
#' which is a useful guard. Unlike a linear ramp, the spline follows the local
#' curvature of the Marey map, which matters in recombination-cold regions where
#' a tiny cM range inverts to a huge bp range (dx/dy can exceed ~40 Mb/cM).
#' Hyman is strictly monotone (flat cM segments stay flat), unlike `monoH.FC`,
#' which overshoots by ~1e-3 at flat->rising transitions. Queries outside the
#' anchor range are clamped to the end values (the `rule = 2` behaviour of the
#' old `approxfun`).
#'
#' @param chr_anchors Clean anchors for ONE chromosome (columns cm, bp).
#'
#' @return A function cM -> bp.
#' @export
make_cm_to_bp <- function(chr_anchors) {
  stopifnot(all(c("cm", "bp") %in% names(chr_anchors)))
  d <- chr_anchors %>%
    dplyr::group_by(.data$cm) %>%
    dplyr::summarise(bp = mean(.data$bp), .groups = "drop") %>%
    dplyr::arrange(.data$cm)
  stopifnot(nrow(d) >= 2)
  f   <- stats::splinefun(d$cm, d$bp, method = "hyman")
  rng <- range(d$cm)
  function(cm) f(pmin(pmax(cm, rng[1]), rng[2]))  # clamp at ends (rule = 2 equiv)
}

#' Convert donor cM segments to physical Mb
#'
#' @param segments_cm Tibble with columns nil_id, chr, start_cm, end_cm.
#' @param anchors_clean Full clean anchor set (columns chr, cm, bp).
#'
#' @return `segments_cm` plus start_bp, end_bp, mb.
#' @export
segments_to_mb <- function(segments_cm, anchors_clean) {
  stopifnot(all(c("chr", "start_cm", "end_cm") %in% names(segments_cm)),
            all(c("chr", "cm", "bp") %in% names(anchors_clean)))

  cm_to_bp <- anchors_clean %>%
    split(.$chr) %>%
    purrr::map(make_cm_to_bp)

  segments_cm %>%
    dplyr::mutate(
      start_bp = purrr::map2_dbl(.data$chr, .data$start_cm,
                                 ~ cm_to_bp[[as.character(.x)]](.y)),
      end_bp   = purrr::map2_dbl(.data$chr, .data$end_cm,
                                 ~ cm_to_bp[[as.character(.x)]](.y)),
      mb       = (.data$end_bp - .data$start_bp) / 1e6
    )
}
