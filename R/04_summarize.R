# Stage 4 — Distribution + summaries
#
# Functions only. See spec section 7 and the validation asserts in section 9.

#' Per-chromosome physical lengths (Mb) from clean anchors
#'
#' @param anchors_clean Clean anchor set (columns chr, bp).
#'
#' @return Tibble: chr, chr_mb (max bp span observed in the map, in Mb).
#' @export
chr_phys_lengths <- function(anchors_clean) {
  stopifnot(all(c("chr", "bp") %in% names(anchors_clean)))
  anchors_clean %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(chr_mb = max(.data$bp) / 1e6, .groups = "drop")
}

#' Assert no donor segment exceeds its chromosome's physical length
#'
#' @param segments_mb Tibble with columns chr, mb.
#' @param anchors_clean Clean anchor set (for per-chr Mb spans).
#'
#' @return Invisibly TRUE if the assertion holds; errors otherwise.
#' @export
assert_within_chr_length <- function(segments_mb, anchors_clean) {
  lens <- chr_phys_lengths(anchors_clean)
  chk <- segments_mb %>%
    dplyr::left_join(lens, by = "chr") %>%
    dplyr::filter(.data$mb > .data$chr_mb + 1e-6)
  if (nrow(chk) > 0) {
    stop("Segment(s) exceed chromosome physical length: ",
         paste0("chr", chk$chr[1], " mb=", round(chk$mb[1], 2),
                " > ", round(chk$chr_mb[1], 2)), call. = FALSE)
  }
  invisible(TRUE)
}

#' Crossover summary for RTIGER autotune
#'
#' Aggregates detectable crossover (donor<->recurrent transition) counts into
#' the statistics RTIGER's autotune needs: COs per diploid genome and COs per
#' megabase. RTIGER asks for the *highest* CO/Mb across samples as the
#' conservative `crossovers_per_megabase` input, so both mean and max are
#' returned.
#'
#' @param co_per_nil Tibble: nil_id, n_co (total detectable COs per genome).
#' @param genome_mb Total physical genome length in Mb.
#'
#' @return A one-row tibble with mean/sd/median/min/max COs per genome and the
#'   mean and (conservative) max COs per Mb.
#' @export
crossover_summary <- function(co_per_nil, genome_mb) {
  stopifnot(all(c("nil_id", "n_co") %in% names(co_per_nil)), genome_mb > 0)
  co_per_nil <- dplyr::mutate(co_per_nil, co_per_mb = .data$n_co / genome_mb)
  tibble::tibble(
    genome_mb        = genome_mb,
    mean_co_genome   = mean(co_per_nil$n_co),
    sd_co_genome     = stats::sd(co_per_nil$n_co),
    median_co_genome = stats::median(co_per_nil$n_co),
    min_co_genome    = min(co_per_nil$n_co),
    max_co_genome    = max(co_per_nil$n_co),
    mean_co_per_mb   = mean(co_per_nil$co_per_mb),
    max_co_per_mb    = max(co_per_nil$co_per_mb)
  )
}

#' Per-NIL summary (one row per line)
#'
#' Includes lines with zero donor genome (left join against the full id set),
#' since the unconditioned BC2S3 distribution legitimately contains NILs with
#' little or no introgression (spec section 8.1).
#'
#' @param segments_mb Tibble with columns nil_id, mb.
#' @param n_nil Total number of simulated NILs.
#'
#' @return Tibble: nil_id, n_segments, total_mb.
#' @export
summarise_nils <- function(segments_mb, n_nil) {
  stopifnot(all(c("nil_id", "mb") %in% names(segments_mb)))
  per_line <- segments_mb %>%
    dplyr::group_by(.data$nil_id) %>%
    dplyr::summarise(n_segments = dplyr::n(),
                     total_mb = sum(.data$mb), .groups = "drop")
  tibble::tibble(nil_id = seq_len(n_nil)) %>%
    dplyr::left_join(per_line, by = "nil_id") %>%
    dplyr::mutate(
      n_segments = tidyr::replace_na(.data$n_segments, 0L),
      total_mb   = tidyr::replace_na(.data$total_mb, 0)
    )
}
