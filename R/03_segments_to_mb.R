# Stage 3 — cM -> Mb conversion
#
# Functions only. See spec section 6.

#' Build a monotone cM -> bp interpolator for one chromosome
#'
#' Collapses cM ties to mean bp (so the relation is a function), then returns a
#' linear interpolator clamped at the ends (`rule = 2`).
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
  stats::approxfun(d$cm, d$bp, rule = 2)  # rule = 2 clamps at ends
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
