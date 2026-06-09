# Stage 6 (benchmark) — score RTIGER calls against simulated ground truth
#
# State encoding (internal, 0/1/2):
#   0 = homozygous recurrent (B73 / REF / RTIGER "AA" / "pat")
#   1 = heterozygous          (RTIGER "AB" / "het")
#   2 = homozygous donor       (teosinte / ALT / RTIGER "BB" / "mat")

# RTIGER label -> internal state ------------------------------------------------
.score_map <- c(AA = 0L, AB = 1L, BB = 2L)          # CompleteBlock .bed scores
.vit_map   <- c(pat = 0L, het = 1L, mat = 2L)       # @Viterbi metadata states

#' Read RTIGER called segments into a tidy table
#'
#' Accepts either (a) a directory containing per-sample
#' `*/CompleteBlock-state-*.bed` files (written by RTIGER save.results), or
#' (b) a fitted RTIGER object (or path to an .rds holding one).
#'
#' @param x Directory path, RTIGER object, or .rds path.
#'
#' @return Tibble: name, chr, start_bp, end_bp, state (0/1/2), sorted.
#' @export
read_rtiger_segments <- function(x) {
  if (is.character(x) && length(x) == 1 && fs::file_exists(x) &&
      grepl("\\.rds$", x, ignore.case = TRUE)) {
    x <- readRDS(x)
  }
  if (methods::is(x, "RTIGER")) return(segments_from_rtiger_object(x))
  stopifnot(fs::dir_exists(x))

  beds <- fs::dir_ls(x, recurse = TRUE, glob = "*CompleteBlock-state-*.bed")
  stopifnot(length(beds) > 0)
  purrr::map(beds, function(f) {
    nm <- stringr::str_remove(stringr::str_remove(fs::path_file(f),
            "^CompleteBlock-state-"), "\\.bed$")
    d <- readr::read_tsv(f, col_names = c("chr", "start0", "end", "score"),
                         show_col_types = FALSE, progress = FALSE)
    tibble::tibble(name = nm, chr = d$chr, start_bp = d$start0 + 1L,
                   end_bp = d$end, state = unname(.score_map[d$score]))
  }) %>%
    purrr::list_rbind() %>%
    dplyr::arrange(.data$name, .data$chr, .data$start_bp)
}

#' Extract called segments from a fitted RTIGER object
#'
#' @param res A fitted RTIGER S4 object (slot `Viterbi`).
#'
#' @return Tibble: name, chr, start_bp, end_bp, state (0/1/2).
#' @export
segments_from_rtiger_object <- function(res) {
  samples <- names(res@Viterbi)
  purrr::map(samples, function(s) {
    gr <- res@Viterbi[[s]]
    d <- tibble::tibble(
      chr   = as.character(GenomicRanges::seqnames(gr)),
      start = GenomicRanges::start(gr),
      end   = GenomicRanges::end(gr),
      state = unname(.vit_map[as.character(
        S4Vectors::elementMetadata(gr)$Viterbi)])
    ) %>% dplyr::arrange(.data$chr, .data$start)
    # collapse consecutive equal-state markers into segments per chr
    d %>%
      dplyr::group_by(.data$chr) %>%
      dplyr::mutate(run = cumsum(c(TRUE, diff(.data$state) != 0))) %>%
      dplyr::group_by(.data$chr, .data$run, .data$state) %>%
      dplyr::summarise(start_bp = min(.data$start), end_bp = max(.data$end),
                       .groups = "drop") %>%
      dplyr::transmute(name = s, .data$chr, .data$start_bp, .data$end_bp,
                       .data$state)
  }) %>%
    purrr::list_rbind() %>%
    dplyr::arrange(.data$name, .data$chr, .data$start_bp)
}

#' Assign each truth marker the RTIGER-called state (overlap join)
#'
#' @param called Tidy called segments (name, chr, start_bp, end_bp, state).
#' @param truth_markers Per-marker truth (name, chr, bp, dosage).
#'
#' @return `truth_markers` plus `called_state` (NA where uncovered).
#' @export
assign_called_state <- function(called, truth_markers) {
  ct <- data.table::as.data.table(called)
  ct[, chr := as.character(chr)]
  data.table::setkey(ct, name, chr, start_bp, end_bp)
  mk <- data.table::as.data.table(truth_markers)
  mk[, `:=`(chr = as.character(chr), start_bp = bp, end_bp = bp)]
  ov <- data.table::foverlaps(mk, ct, by.x = c("name", "chr", "start_bp",
        "end_bp"), by.y = c("name", "chr", "start_bp", "end_bp"),
        type = "within", mult = "first", nomatch = NA)
  tibble::as_tibble(ov) %>%
    dplyr::transmute(.data$name, .data$chr, .data$bp, true_state = .data$dosage,
                     called_state = .data$state)
}

#' Confusion matrix + per-state metrics at marker level
#'
#' @param mk Output of `assign_called_state()`.
#'
#' @return List: `confusion` (table), `accuracy`, `per_state` tibble
#'   (state, recall, precision), and `uncovered` fraction.
#' @export
score_markers <- function(mk) {
  cov <- dplyr::filter(mk, !is.na(.data$called_state))
  cm <- table(true = factor(cov$true_state, 0:2),
              called = factor(cov$called_state, 0:2))
  per_state <- purrr::map_dfr(0:2, function(s) {
    tp <- cm[as.character(s), as.character(s)]
    tibble::tibble(state = s,
                   recall    = tp / sum(cm[as.character(s), ]),
                   precision = tp / sum(cm[, as.character(s)]))
  })
  list(confusion = cm, accuracy = sum(diag(cm)) / sum(cm),
       per_state = per_state, uncovered = mean(is.na(mk$called_state)))
}

#' Segment-level precision/recall/FDR for an event class via reciprocal overlap
#'
#' An event is a maximal run whose state is in `target`. A truth event is a TP
#' if some called event reciprocally overlaps it by >= `min_overlap` (both
#' directions). Boundary errors are reported for matched pairs.
#'
#' @param called,truth Tidy segment tables (name, chr, start_bp, end_bp, state).
#' @param target Integer states defining the event class (e.g. 2L, or c(1L,2L)).
#' @param min_overlap Reciprocal-overlap fraction to count a match (default 0.5).
#' @param label A name for the class (for the output row).
#'
#' @return List: `summary` one-row tibble (TP/FP/FN/precision/recall/F1/FDR,
#'   boundary error) and `matched` pairs tibble.
#' @export
score_segments <- function(called, truth, target, min_overlap = 0.5,
                           label = "event") {
  events <- function(df) {
    df %>% dplyr::filter(.data$state %in% target) %>%
      dplyr::mutate(len = .data$end_bp - .data$start_bp + 1)
  }
  te <- events(truth); ce <- events(called)
  if (nrow(te) == 0 && nrow(ce) == 0)
    return(list(summary = .empty_seg_summary(label), matched = tibble::tibble()))

  ov <- .pair_overlaps(te, ce)            # all overlapping (truth,called) pairs
  matched <- ov %>%
    dplyr::filter(.data$ov_len / .data$len_t >= min_overlap,
                  .data$ov_len / .data$len_c >= min_overlap)
  tp_truth <- dplyr::distinct(matched, .data$ti)
  tp_call  <- dplyr::distinct(matched, .data$ci)
  TP <- nrow(tp_truth)
  FN <- nrow(te) - TP
  FP <- nrow(ce) - nrow(tp_call)
  prec <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  rec  <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(rec) && prec + rec > 0)
            2 * prec * rec / (prec + rec) else NA_real_
  best <- matched %>% dplyr::group_by(.data$ti) %>%
    dplyr::slice_max(.data$ov_len, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
  list(
    summary = tibble::tibble(
      class = label, n_truth = nrow(te), n_called = nrow(ce),
      TP = TP, FP = FP, FN = FN,
      precision = prec, recall = rec, F1 = f1,
      FDR = if (TP + FP > 0) FP / (TP + FP) else NA_real_,
      med_start_err_kb = stats::median(abs(best$start_t - best$start_c)) / 1e3,
      med_end_err_kb   = stats::median(abs(best$end_t - best$end_c)) / 1e3
    ),
    matched = best)
}

.empty_seg_summary <- function(label) tibble::tibble(
  class = label, n_truth = 0L, n_called = 0L, TP = 0L, FP = 0L, FN = 0L,
  precision = NA_real_, recall = NA_real_, F1 = NA_real_, FDR = NA_real_,
  med_start_err_kb = NA_real_, med_end_err_kb = NA_real_)

# all overlapping (truth, called) event pairs within the same name+chr
.pair_overlaps <- function(te, ce) {
  t_dt <- data.table::as.data.table(te)[, .(ti = .I, name, chr = as.character(chr),
            start_t = start_bp, end_t = end_bp, len_t = len)]
  c_dt <- data.table::as.data.table(ce)[, .(ci = .I, name, chr = as.character(chr),
            start_c = start_bp, end_c = end_bp, len_c = len)]
  data.table::setkey(c_dt, name, chr, start_c, end_c)
  ov <- data.table::foverlaps(t_dt, c_dt,
          by.x = c("name", "chr", "start_t", "end_t"),
          by.y = c("name", "chr", "start_c", "end_c"),
          type = "any", nomatch = NULL)
  ov[, ov_len := pmin(end_t, end_c) - pmax(start_t, start_c) + 1]
  tibble::as_tibble(ov)
}

#' Crossover-count comparison (called vs true detectable COs)
#'
#' @param called,truth Tidy segment tables.
#'
#' @return Tibble: name, true_co, called_co (state transitions = n_seg - 1
#'   summed over chromosomes).
#' @export
score_cos <- function(called, truth) {
  co <- function(df, col) df %>% dplyr::group_by(.data$name, .data$chr) %>%
    dplyr::summarise(n = dplyr::n() - 1, .groups = "drop") %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarise(!!col := sum(.data$n), .groups = "drop")
  dplyr::full_join(co(truth, "true_co"), co(called, "called_co"), by = "name") %>%
    tidyr::replace_na(list(true_co = 0, called_co = 0))
}

#' Degrade truth into a mock RTIGER call set (for self-testing the scorer)
#'
#' Merges segments shorter than `min_markers` into the preceding segment
#' (mimicking RTIGER's rigidity floor) and jitters breakpoints, so scoring on
#' the result yields recall < 1 and nonzero boundary error.
#'
#' @param truth_seg Truth segments (name, chr, start_bp, end_bp, state, n_markers).
#' @param min_markers Segments below this are absorbed into the previous one.
#' @param jitter_bp Uniform breakpoint jitter (+/-).
#'
#' @return A called-segment tibble (name, chr, start_bp, end_bp, state).
#' @export
degrade_truth <- function(truth_seg, min_markers = 30L, jitter_bp = 50000L) {
  truth_seg %>%
    dplyr::group_by(.data$name, .data$chr) %>%
    dplyr::group_modify(function(g, ...) {
      st <- g$state
      keep <- g$n_markers >= min_markers
      keep[1] <- TRUE
      for (i in seq_along(st)) if (!keep[i]) st[i] <- st[i - 1]
      g$state <- st
      # re-collapse runs
      g$run <- cumsum(c(TRUE, diff(g$state) != 0))
      g %>% dplyr::group_by(.data$run) %>%
        dplyr::summarise(start_bp = min(.data$start_bp),
                         end_bp = max(.data$end_bp),
                         state = dplyr::first(.data$state), .groups = "drop") %>%
        dplyr::select(-run)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$name, .data$chr) %>%
    dplyr::mutate(
      start_bp = ifelse(dplyr::row_number() == 1, .data$start_bp,
                        .data$start_bp + round(stats::runif(dplyr::n(),
                          -jitter_bp, jitter_bp))),
      start_bp = pmax(1, .data$start_bp)) %>%
    dplyr::ungroup() %>%
    dplyr::select("name", "chr", "start_bp", "end_bp", "state")
}
