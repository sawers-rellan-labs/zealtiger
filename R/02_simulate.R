# Stage 2 — Pedigree + simulation (BC2S3 NILs, no selection)
#
# Functions only. See spec sections 5 and 9.
#
# API notes confirmed against simcross 0.8 in this environment:
#   * check_pedigree() must be called with ignore_sex = TRUE: selfing is
#     encoded as mom == dad, which makes a single parent both "mom" and "dad",
#     so the default male/female consistency check fails. Sex is immaterial for
#     autosomal maize (spec section 5).
#   * sim_from_pedigree(ped, L = <vector>, m, p) with a length-K L returns a
#     list indexed by CHROMOSOME (1..K); each element is a list indexed by
#     individual id ("1".."8"); each of those has $mat / $pat with $alleles and
#     $locations. So the NIL's chr-c haplotypes are x[[c]][["8"]].

#' Build the BC2S3 NIL pedigree (single replicate)
#'
#' Two inbred founders (1 = recurrent, 2 = donor), F1, two backcrosses to the
#' recurrent founder, then three generations of selfing. Final NIL is id 8.
#'
#' @return A base data frame with columns id, mom, dad, sex, gen. A base
#'   data.frame (not a tibble) is required: simcross's check_pedigree() uses
#'   `pedigree$col` in a way that returns a 1-column tibble for tibble input,
#'   breaking its internal id matching.
#' @export
bc2s3_pedigree <- function() {
  data.frame(
    id  = 1:8,
    mom = c(0, 0, 1, 1, 1, 5, 6, 7),
    dad = c(0, 0, 2, 3, 4, 5, 6, 7),
    sex = c(0, 1, 0, 1, 0, 1, 0, 1),
    gen = c(0, 0, 1, 2, 3, 4, 5, 6)
  )
}

#' Extract donor (founder-2) segments for one chromosome of one NIL
#'
#' Takes the simcross result for the final individual on one chromosome and
#' returns donor intervals (per haplotype) in cM.
#'
#' @param nil_chr One chromosome's result for the NIL: list with `$mat` and
#'   `$pat`, each having `$alleles` and `$locations`.
#' @param donor_allele Founder label for the donor (default 2).
#'
#' @return Tibble: hap, start_cm, end_cm for donor-carrying segments.
#' @export
donor_segments_cm <- function(nil_chr, donor_allele = 2L) {
  stopifnot(all(c("mat", "pat") %in% names(nil_chr)))

  one_hap <- function(hap, label) {
    locs <- hap$locations
    alle <- hap$alleles
    left <- c(0, utils::head(locs, -1))
    tibble::tibble(
      hap      = label,
      start_cm = left,
      end_cm   = locs,
      allele   = alle
    ) %>%
      dplyr::filter(.data$allele == donor_allele) %>%
      dplyr::select(-allele)
  }

  dplyr::bind_rows(one_hap(nil_chr$mat, "mat"),
                   one_hap(nil_chr$pat, "pat"))
}

#' Union of (possibly overlapping) intervals
#'
#' Collapses a set of [start, end] intervals into non-overlapping intervals.
#' Used to merge donor segments across the mat and pat haplotypes: a region is
#' "introgressed" if either haplotype carries the donor allele.
#'
#' @param start Numeric vector of interval starts.
#' @param end Numeric vector of interval ends (same length as `start`).
#'
#' @return Tibble: start, end (sorted, non-overlapping).
#' @export
union_intervals <- function(start, end) {
  stopifnot(length(start) == length(end))
  if (length(start) == 0L) {
    return(tibble::tibble(start = numeric(0), end = numeric(0)))
  }
  o <- order(start)
  s <- start[o]
  e <- end[o]

  out_s <- numeric(0)
  out_e <- numeric(0)
  cur_s <- s[1]
  cur_e <- e[1]
  for (i in seq_along(s)[-1]) {
    if (s[i] <= cur_e) {
      cur_e <- max(cur_e, e[i])
    } else {
      out_s <- c(out_s, cur_s)
      out_e <- c(out_e, cur_e)
      cur_s <- s[i]
      cur_e <- e[i]
    }
  }
  out_s <- c(out_s, cur_s)
  out_e <- c(out_e, cur_e)
  tibble::tibble(start = out_s, end = out_e)
}

#' Count crossover (allele-transition) events on one haplotype
#'
#' With two founders, every boundary between adjacent mosaic intervals of
#' different founder origin is a detectable donor<->recurrent transition, i.e.
#' a crossover as RTIGER would call it from resequencing. Same-allele adjacency
#' (if simcross ever records it) is not counted.
#'
#' @param hap A haplotype list with `$alleles`.
#'
#' @return Integer count of allele transitions.
#' @export
n_crossovers_hap <- function(hap) {
  a <- hap$alleles
  if (length(a) < 2L) return(0L)
  sum(a[-1] != a[-length(a)])
}

#' Simulate one BC2S3 NIL and return its donor introgressions + crossovers
#'
#' Runs `sim_from_pedigree` once (all chromosomes via a length-K cM vector),
#' pulls the final individual (id 8), and for each chromosome returns the
#' donor intervals as the UNION across the mat and pat haplotypes (the
#' physically-detected introgression model; spec section 5 note).
#'
#' Also reports:
#'   * `donor_cm_hap`: total donor cM over BOTH haplotypes (pre-union) — the
#'     dosage that validates against the allele frequency 0.125 (spec §9).
#'   * `co_chr`: per-chromosome crossover count = mat + pat allele transitions.
#'     This is the detectable-CO count per diploid genome that feeds RTIGER's
#'     `crossovers_per_megabase` autotune parameter.
#'
#' @param ped Pedigree from `bc2s3_pedigree()`.
#' @param chr_cm_lengths Tibble with columns chr, L (max cM per chromosome),
#'   arranged by chr.
#' @param m,p Stahl interference parameters (simcross default m = 10, p = 0).
#' @param donor_allele Founder label for the donor (default 2).
#' @param nil_id Individual id (as character) of the final NIL (default "8").
#'
#' @return A list with `$segments` (tibble: chr, start_cm, end_cm — unioned
#'   donor introgressions), `$donor_cm_hap` (scalar), and `$co_chr` (tibble:
#'   chr, n_co).
#' @export
simulate_nil <- function(ped, chr_cm_lengths, m = 10, p = 0,
                         donor_allele = 2L, nil_id = "8") {
  stopifnot(all(c("chr", "L") %in% names(chr_cm_lengths)))
  L <- chr_cm_lengths$L
  sim <- simcross::sim_from_pedigree(ped, L = L, m = m, p = p)

  per_chr <- purrr::map(seq_along(L), function(ci) {
    nil_chr <- sim[[ci]][[nil_id]]
    hap_segs <- donor_segments_cm(nil_chr, donor_allele = donor_allele)
    donor_cm_hap <- sum(hap_segs$end_cm - hap_segs$start_cm)
    u <- union_intervals(hap_segs$start_cm, hap_segs$end_cm)
    n_co <- n_crossovers_hap(nil_chr$mat) + n_crossovers_hap(nil_chr$pat)
    list(
      segments = tibble::tibble(
        chr      = chr_cm_lengths$chr[ci],
        start_cm = u$start,
        end_cm   = u$end
      ),
      donor_cm_hap = donor_cm_hap,
      co = tibble::tibble(chr = chr_cm_lengths$chr[ci], n_co = n_co)
    )
  })

  list(
    segments     = purrr::list_rbind(purrr::map(per_chr, "segments")),
    donor_cm_hap = sum(purrr::map_dbl(per_chr, "donor_cm_hap")),
    co_chr       = purrr::list_rbind(purrr::map(per_chr, "co"))
  )
}
