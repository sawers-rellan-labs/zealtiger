# Stage 1 — Compile the maize map (B73 NAM v5)
#
# Functions only; the orchestration lives in nil_introgression_1400.qmd / run_all.R.

#' Parse one MaizeGDB consensus-map chromosome file to v5 anchors
#'
#' Reads a single `consensus_chrNN.txt`, locates the B73 NAM v5 columns by
#' name (robust to per-file column differences and dirty whitespace), and
#' returns tidy (locus, cM, chr, bp) anchors with physical coordinates.
#'
#' @param path Path to a consensus_chrNN.txt file.
#'
#' @return A tibble with columns locus, cm, chr (int), bp (numeric).
#' @export
parse_consensus_v5 <- function(path) {
  stopifnot(fs::file_exists(path))

  raw <- readr::read_tsv(path, skip = 1, show_col_types = FALSE)
  names(raw) <- stringr::str_squish(names(raw))

  v5 <- "Zm-B73-REFERENCE-NAM-5.0"
  col_chr   <- paste0(v5, "_chr")
  col_start <- paste0(v5, "_start")
  col_end   <- paste0(v5, "_end")

  stopifnot(all(c("Locus", "Coordinate", col_chr, col_start, col_end)
                %in% names(raw)))

  raw %>%
    dplyr::transmute(
      locus = .data$Locus,
      cm    = suppressWarnings(as.numeric(.data$Coordinate)),
      chr   = suppressWarnings(
        as.integer(stringr::str_remove(.data[[col_chr]], "(?i)^chr"))),
      start = suppressWarnings(as.numeric(.data[[col_start]])),
      end   = suppressWarnings(as.numeric(.data[[col_end]]))
    ) %>%
    dplyr::filter(
      !is.na(.data$cm), !is.na(.data$chr),
      !is.na(.data$start), !is.na(.data$end)
    ) %>%
    dplyr::mutate(bp = (.data$start + .data$end) / 2) %>%
    dplyr::select(locus, cm, chr, bp)
}

#' Enforce a monotone Marey map within each chromosome
#'
#' Sorts markers by physical position and retains only those whose genetic
#' position is non-decreasing (running-max filter), removing map inversions
#' from mismapped or paralogous loci.
#'
#' @param anchors Tibble from `parse_consensus_v5()` (possibly many chrs).
#'
#' @return Tibble of the same columns, monotone in (bp, cm) per chr.
#' @export
enforce_monotone <- function(anchors) {
  stopifnot(is.data.frame(anchors),
            all(c("chr", "cm", "bp") %in% names(anchors)))

  anchors %>%
    dplyr::arrange(.data$chr, .data$bp) %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::mutate(cm_runmax = cummax(.data$cm)) %>%
    dplyr::filter(.data$cm >= .data$cm_runmax) %>%
    dplyr::select(-cm_runmax) %>%
    dplyr::ungroup()
}

#' Thin clean anchors to a ~n_target marker grid
#'
#' Allocates a per-chromosome quota proportional to each chromosome's cM span,
#' then selects markers as evenly spaced in cM as possible (nearest available
#' marker to evenly-spaced cM targets). If a chromosome has fewer clean markers
#' than its quota, all are kept (2500 is a cap/target, not a guarantee).
#'
#' @param anchors Monotone anchor tibble (locus, cm, chr, bp).
#' @param n_target Target total marker count (default 2500).
#'
#' @return A thinned tibble, same columns as `anchors`.
#' @export
thin_markers <- function(anchors, n_target = 2500) {
  stopifnot(is.data.frame(anchors),
            all(c("locus", "cm", "chr", "bp") %in% names(anchors)),
            n_target > 0)

  spans <- anchors %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(span = max(.data$cm) - min(.data$cm),
                     n_avail = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      quota = pmax(2L, as.integer(round(n_target * .data$span / sum(.data$span))))
    )

  anchors %>%
    dplyr::arrange(.data$chr, .data$cm) %>%
    dplyr::group_split(.data$chr) %>%
    purrr::map(function(df) {
      q <- spans$quota[spans$chr == df$chr[1]]
      if (nrow(df) <= q) return(df)
      targets <- seq(min(df$cm), max(df$cm), length.out = q)
      idx <- vapply(targets, function(t) which.min(abs(df$cm - t)), integer(1))
      df[sort(unique(idx)), ]
    }) %>%
    purrr::list_rbind()
}

#' Per-chromosome map QC report
#'
#' @param raw Raw anchors (pre-monotonicity) from `parse_consensus_v5()`.
#' @param clean Monotone anchors from `enforce_monotone()`.
#'
#' @return A tibble: chr, n_raw, n_clean, dropped, pct_dropped, cm_span,
#'   mb_span, cm_per_mb.
#' @export
report_map_qc <- function(raw, clean) {
  stopifnot(is.data.frame(raw), is.data.frame(clean))

  clean_stats <- clean %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(
      n_clean = dplyr::n(),
      cm_span = max(.data$cm) - min(.data$cm),
      mb_span = (max(.data$bp) - min(.data$bp)) / 1e6,
      .groups = "drop"
    )

  raw %>%
    dplyr::count(.data$chr, name = "n_raw") %>%
    dplyr::left_join(clean_stats, by = "chr") %>%
    dplyr::mutate(
      dropped     = .data$n_raw - .data$n_clean,
      pct_dropped = round(100 * .data$dropped / .data$n_raw, 1),
      cm_per_mb   = round(.data$cm_span / .data$mb_span, 3)
    ) %>%
    dplyr::arrange(.data$chr)
}
