#!/usr/bin/env Rscript
# Driver for the BC2S3 NIL introgression-size simulation.
#
# Usage:
#   Rscript run_all.R            # Stage 1 only (map compile + Marey QC gate)
#   Rscript run_all.R --sim      # Stages 1-4 (full 1400-NIL run)
#   Rscript run_all.R --pilot    # Stage 1 + n=20 pilot (dosage sanity check)
#
# Run from the nil_introgression/ project directory.

suppressMessages({
  library(tidyverse)
  library(simcross)
  library(furrr)
  library(patchwork)
})
purrr::walk(fs::dir_ls("R", glob = "*.R"), source)

# ---- parameters (mirror nil_introgression.qmd) ----------------------------
params <- list(
  map_dir        = "data/maizegdb_consensus_map",
  n_nil          = 1500L,
  n_markers      = 2500L,
  m_interference = 10,    # simcross default; positive interference
  p_escape       = 0,
  donor_allele   = 2L,
  seed           = 20260609L
)

args   <- commandArgs(trailingOnly = TRUE)
do_sim   <- "--sim"   %in% args
do_pilot <- "--pilot" %in% args || do_sim

set.seed(params$seed)
fs::dir_create("results/figures")

map_dir <- params$map_dir
stopifnot(fs::dir_exists(map_dir))

# ===========================================================================
# Stage 1 — Compile the map
# ===========================================================================
map_files <- fs::dir_ls(map_dir, glob = "*chr*.txt")

anchors_raw <- map_files %>%
  purrr::map(parse_consensus_v5) %>%
  purrr::list_rbind() %>%
  dplyr::arrange(chr, bp)

anchors_clean <- enforce_monotone(anchors_raw)
qc <- report_map_qc(anchors_raw, anchors_clean)

marker_grid <- thin_markers(anchors_clean, n_target = params$n_markers)

chr_cm_lengths <- anchors_clean %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(L = max(cm), .groups = "drop") %>%
  dplyr::arrange(chr)

saveRDS(anchors_clean, "results/maize_map_v5_clean.rds")
saveRDS(marker_grid,   "results/maize_map_v5_grid2500.rds")

marey <- anchors_clean %>%
  ggplot(aes(bp / 1e6, cm)) +
  geom_point(size = 0.4, alpha = 0.5, colour = "grey30") +
  facet_wrap(~ chr, scales = "free") +
  labs(x = "Physical position (Mb, B73 NAM v5)", y = "Genetic position (cM)",
       title = "Marey maps (monotone clean anchors)") +
  theme_bw()
ggsave("results/marey_maps.pdf", marey, width = 11, height = 8)

cat("\n================ STAGE 1: MAP QC ================\n")
cat("Files parsed:", length(map_files),
    "| raw anchors:", nrow(anchors_raw),
    "| clean anchors:", nrow(anchors_clean), "\n")
cat("Marker grid:", nrow(marker_grid), "(target", params$n_markers, ")\n")
cat("Total map length (cM):", round(sum(chr_cm_lengths$L), 1),
    "(expected maize range ~1400-1600)\n\n")
print(qc)
cat("\nPer-chr cM length (L for simcross):\n")
print(chr_cm_lengths)
cat("\nWrote: results/maize_map_v5_clean.rds, ",
    "results/maize_map_v5_grid2500.rds, results/marey_maps.pdf\n", sep = "")

if (!do_pilot) {
  cat("\n>>> Stage 1 only. Inspect results/marey_maps.pdf, then re-run with",
      "--sim. <<<\n")
  quit(save = "no")
}

# ===========================================================================
# Stage 2 — Pedigree + simulation
# ===========================================================================
ped <- bc2s3_pedigree()
stopifnot(check_pedigree(ped, ignore_sex = TRUE))  # selfing => ignore_sex

plan(multisession)
total_cm  <- sum(chr_cm_lengths$L)
genome_mb <- sum(chr_phys_lengths(anchors_clean)$chr_mb)  # map-spanned Mb

sim_one <- function(i) {
  res <- simulate_nil(ped, chr_cm_lengths,
                      m = params$m_interference, p = params$p_escape,
                      donor_allele = params$donor_allele)
  res$segments$nil_id <- i
  list(segments = res$segments, donor_cm_hap = res$donor_cm_hap)
}

# ---- pilot (n = 20): dosage must approach 0.125 ----------------------------
pilot <- furrr::future_map(seq_len(20), sim_one,
                           .options = furrr_options(seed = TRUE))
# dosage = donor cM over BOTH haplotypes / (2 * genome cM) = allele frequency
pilot_dosage <- purrr::map_dbl(pilot, "donor_cm_hap") / (2 * total_cm)
cat("\n================ STAGE 2: PILOT (n=20) ================\n")
cat(sprintf("Mean donor dosage: %.4f  (analytical expectation 0.1250)\n",
            mean(pilot_dosage)))
stopifnot(abs(mean(pilot_dosage) - 0.125) < 0.03)

if (!do_sim) {
  cat("\n>>> Pilot passed. Re-run with --sim for the full 1400-NIL run. <<<\n")
  quit(save = "no")
}

# ---- full run, reusable over interference settings ------------------------
# Runs n_nil NILs at (m, p), then Stages 3-4. Returns segments_mb, nil_summary,
# dosage. Used for the primary (m = 10) run and the m = 0 sensitivity pass.
run_full <- function(m, p, label) {
  sims <- furrr::future_map(
    seq_len(params$n_nil),
    function(i) {
      res <- simulate_nil(ped, chr_cm_lengths, m = m, p = p,
                          donor_allele = params$donor_allele)
      res$segments$nil_id <- i
      res
    },
    .options = furrr_options(seed = TRUE)
  )
  segments_cm <- purrr::list_rbind(purrr::map(sims, "segments"))
  dosage <- purrr::map_dbl(sims, "donor_cm_hap") / (2 * total_cm)

  segments_mb <- segments_to_mb(
    dplyr::select(segments_cm, nil_id, chr, start_cm, end_cm), anchors_clean)
  assert_within_chr_length(segments_mb, anchors_clean)
  nil_summary <- summarise_nils(segments_mb, params$n_nil)

  # detectable crossovers per genome (for RTIGER crossovers_per_megabase)
  co_per_nil <- purrr::imap(sims, ~ dplyr::mutate(.x$co_chr, nil_id = .y)) %>%
    purrr::list_rbind()
  co_genome <- co_per_nil %>%
    dplyr::group_by(nil_id) %>%
    dplyr::summarise(n_co = sum(n_co), .groups = "drop")
  co_stats <- crossover_summary(co_genome, genome_mb)

  cat(sprintf(
    "[%s | m=%g p=%g] dosage %.4f (exp 0.1250) | segs %d | seg Mb med %.2f / 90th %.2f / max %.2f | tot Mb/line %.1f | CO/genome mean %.1f max %d | CO/Mb mean %.4f max %.4f\n",
    label, m, p, mean(dosage), nrow(segments_mb),
    median(segments_mb$mb), quantile(segments_mb$mb, 0.9),
    max(segments_mb$mb), mean(nil_summary$total_mb),
    co_stats$mean_co_genome, co_stats$max_co_genome,
    co_stats$mean_co_per_mb, co_stats$max_co_per_mb))

  list(label = label, m = m, p = p, segments_mb = segments_mb,
       nil_summary = nil_summary, dosage = dosage,
       co_per_chr = co_per_nil, co_genome = co_genome, co_stats = co_stats)
}

cat("\n================ STAGE 2-4: FULL RUNS ================\n")
primary <- run_full(params$m_interference, params$p_escape, "default")
noint   <- run_full(0, 0, "no_interference")

# ---- canonical outputs come from the primary (default-interference) run ----
segments_mb <- primary$segments_mb
nil_summary <- primary$nil_summary
q90 <- quantile(segments_mb$mb, 0.9)

saveRDS(segments_mb, "results/nil_segments.rds")
readr::write_csv(segments_mb, "results/nil_segments.csv")
saveRDS(nil_summary, "results/nil_summary.rds")
readr::write_csv(nil_summary, "results/nil_summary.csv")

# ---- crossover counts per genome (for RTIGER crossovers_per_megabase) ------
co_genome   <- primary$co_genome                       # one row per NIL
co_per_chr  <- primary$co_per_chr                       # nil_id x chr counts
co_stats    <- primary$co_stats
co_chr_stats <- co_per_chr %>%                          # per-chromosome CO dist
  dplyr::group_by(chr) %>%
  dplyr::summarise(mean_co = mean(n_co), max_co = max(n_co), .groups = "drop") %>%
  dplyr::left_join(chr_phys_lengths(anchors_clean), by = "chr") %>%
  dplyr::mutate(mean_co_per_mb = mean_co / chr_mb,
                max_co_per_mb  = max_co  / chr_mb)
readr::write_csv(dplyr::mutate(co_genome, co_per_mb = n_co / genome_mb),
                 "results/nil_crossovers.csv")
readr::write_csv(co_chr_stats, "results/crossovers_per_chr.csv")

p_total <- ggplot(nil_summary, aes(total_mb)) +
  geom_histogram(bins = 40, fill = "steelblue") +
  labs(x = "Total donor genome per NIL (Mb)", y = "NILs") + theme_bw()
p_seg <- ggplot(segments_mb, aes(mb)) +
  geom_histogram(bins = 50, fill = "darkorange") +
  labs(x = "Donor segment size (Mb)", y = "Segments") + theme_bw()
p_ecdf <- ggplot(segments_mb, aes(mb)) +
  stat_ecdf() +
  geom_vline(xintercept = median(segments_mb$mb), linetype = 2) +
  geom_vline(xintercept = q90, linetype = 3) +
  labs(x = "Donor segment size (Mb)", y = "ECDF",
       caption = "dashed = median, dotted = 90th pct") + theme_bw()
dist_panel <- p_total / p_seg / p_ecdf
ggsave("results/figures/introgression_size.pdf", dist_panel, width = 7, height = 10)
ggsave("results/figures/introgression_size.png", dist_panel, width = 7, height = 10,
       dpi = 150)

# ---- interference sensitivity comparison ----------------------------------
cmp <- dplyr::bind_rows(
  dplyr::mutate(segments_mb, model = "default (m=10)"),
  dplyr::mutate(noint$segments_mb, model = "no interference (m=0)")
)
cmp_summary <- cmp %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    n_seg     = dplyr::n(),
    median_mb = median(mb),
    p90_mb    = quantile(mb, 0.9),
    p99_mb    = quantile(mb, 0.99),
    max_mb    = max(mb),
    .groups = "drop"
  )
readr::write_csv(cmp_summary, "results/interference_comparison.csv")

p_cmp <- ggplot(cmp, aes(mb, colour = model)) +
  stat_ecdf() +
  labs(x = "Donor segment size (Mb)", y = "ECDF",
       title = "Interference sensitivity: segment-size tail",
       colour = NULL) +
  theme_bw() + theme(legend.position = "bottom")
ggsave("results/figures/interference_comparison.pdf", p_cmp, width = 7, height = 5)
ggsave("results/figures/interference_comparison.png", p_cmp, width = 7, height = 5,
       dpi = 150)

jsonlite::write_json(
  list(
    seed              = params$seed,
    n_nil             = params$n_nil,
    n_markers         = nrow(marker_grid),
    total_cM          = total_cm,
    genome_Mb         = genome_mb,
    simcross          = as.character(packageVersion("simcross")),
    primary           = list(m = primary$m, p = primary$p,
                             mean_donor_dosage = mean(primary$dosage),
                             mean_total_mb = mean(primary$nil_summary$total_mb),
                             mean_co_per_genome = co_stats$mean_co_genome,
                             max_co_per_genome  = co_stats$max_co_genome,
                             mean_co_per_mb     = co_stats$mean_co_per_mb,
                             max_co_per_mb      = co_stats$max_co_per_mb),
    sensitivity       = list(m = noint$m, p = noint$p,
                             mean_donor_dosage = mean(noint$dosage),
                             mean_total_mb = mean(noint$nil_summary$total_mb),
                             mean_co_per_genome = noint$co_stats$mean_co_genome,
                             max_co_per_genome  = noint$co_stats$max_co_genome,
                             max_co_per_mb      = noint$co_stats$max_co_per_mb),
    rtiger_crossovers_per_megabase = co_stats$max_co_per_mb,
    date              = as.character(Sys.time())
  ),
  "results/run_metadata.json", pretty = TRUE, auto_unbox = TRUE
)

cat("\n================ INTERFERENCE COMPARISON ================\n")
print(cmp_summary)
cat("\n================ CROSSOVERS (primary, m=10) ================\n")
print(co_stats)
cat("\nPer-chromosome CO (mean / max counts, mean / max per Mb across NILs):\n")
print(co_chr_stats)
cat("\nConservative per-chromosome CO/Mb vector (max over samples):\n")
cat(paste(sprintf("chr%d=%.4f", co_chr_stats$chr, co_chr_stats$max_co_per_mb),
          collapse = ", "), "\n")
cat(sprintf(
  "\n>>> RTIGER crossovers_per_megabase (genome-wide conservative): %.4f <<<\n",
  co_stats$max_co_per_mb))
cat("\nWrote nil_segments.{csv,rds}, nil_summary.{csv,rds},",
    "nil_crossovers.csv, crossovers_per_chr.csv,\n",
    "interference_comparison.csv, run_metadata.json,",
    "figures/*.{pdf,png}\n")
