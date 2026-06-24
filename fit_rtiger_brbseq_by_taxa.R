#!/usr/bin/env Rscript
# Fit BrB-RTI PER SPECIES (taxon) on the cached BRB counts, B73/Purple as separate
# check groups — mirrors the SNP50K per-taxon scheme (fit_rtiger_by_taxa.R +
# fit_rtiger_checks.R). One joint EM WITHIN each taxon (no cross-species pooling, which
# mis-specifies emissions and destabilises the B73 control). Downstream code SUBSETS each
# sample's call from the combined output by `name` — never re-fits an ad-hoc cohort.
#
# Input : results/sim_calibration/brbseq_ks/counts/<name>.tsv  (full 379, cached)
#         agent/brb_species_map.csv  (name, sid, accession, taxon, is_check)
#         data/rtiger_50K/seqlengths.csv
# Output: results/sim_calibration/brbseq_taxa_r109/calls_common_schema.csv
#           source = BRBseq_r109_MLE ; donor = taxon (Zx/Zv/Zd/Zl/Zh) or B73/Purple
# Usage : Rscript fit_rtiger_brbseq_by_taxa.R [--r=109]
suppressMessages({library(RTIGER); library(readr); library(tidyverse); library(GenomicRanges); library(S4Vectors)})
source("R/06_score_rtiger.R"); set.seed(42)
JULIA_HOME <- Sys.getenv("JULIA_HOME", "")
if (!nzchar(JULIA_HOME) || !file.exists(file.path(JULIA_HOME, "julia"))) {
  l <- path.expand("~/.juliaup/bin/julia")
  if (file.exists(l)) JULIA_HOME <- system(paste(shQuote(l), "--startup-file=no -e",
                                                  shQuote("print(Sys.BINDIR)")), intern = TRUE)
}
setupJulia(JULIA_HOME = JULIA_HOME); sourceJulia()

argv <- commandArgs(trailingOnly = TRUE)
go <- function(f, d) { h <- grepl(paste0("^", f, "="), argv); if (any(h)) sub(paste0("^", f, "="), "", argv[h][1]) else d }
r       <- as.integer(go("--r", "109"))
cnt     <- "results/sim_calibration/brbseq_ks/counts"
outdir  <- go("--outdir", "results/sim_calibration/brbseq_taxa_r109")
sl   <- read_csv("data/rtiger_50K/seqlengths.csv", show_col_types = FALSE); seqv <- setNames(sl$len, sl$chr_label)

map <- read_csv("agent/brb_species_map.csv", show_col_types = FALSE) %>%
  mutate(group = ifelse(is_check, ifelse(grepl("[Pp]urple", accession), "Purple", "B73"), taxon))

# Coverage QC: RTIGER aborts if any sample has a chromosome with < 2*r covered markers.
covered_ok <- function(file) {
  d <- suppressMessages(read_tsv(file, col_names = FALSE, show_col_types = FALSE))
  cov <- d[d$X4 + d$X6 > 0, ]; per_chr <- table(cov$X1)
  length(per_chr) == 10 && min(per_chr) >= 2L * r
}

groups <- map %>% count(group, name = "n") %>% arrange(desc(n))
cat(sprintf("r=%d | %d groups: %s\n", r, nrow(groups),
            paste(sprintf("%s(%d)", groups$group, groups$n), collapse = " ")))
all_seg <- list(); summ <- list()
for (grp in groups$group) {
  nm <- map$name[map$group == grp]
  files <- normalizePath(file.path(cnt, paste0(nm, ".tsv")), mustWork = FALSE)
  ok <- file.exists(files); nm <- nm[ok]; files <- files[ok]
  keep <- vapply(files, covered_ok, logical(1))
  if (any(!keep)) cat(sprintf("  [%s] coverage QC drops %d/%d: %s\n", grp, sum(!keep), length(keep),
                              paste(nm[!keep], collapse = ", ")))
  nm <- nm[keep]; files <- files[keep]
  if (!length(nm)) { cat(sprintf("  [%s] EMPTY after QC — skipped\n", grp)); next }
  ed <- data.frame(files = files, name = nm)
  fit_dir <- file.path(outdir, "fit", grp); dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  res <- RTIGER(expDesign = ed, outputdir = fit_dir, seqlengths = seqv, rigidity = r,
                autotune = FALSE, save.results = FALSE, threads = 1L, verbose = FALSE)
  seg <- segments_from_rtiger_object(res) %>%
    mutate(source = sprintf("BRBseq_r%d_MLE", r), donor = grp) %>%
    select(source, donor, name, chr, start_bp, end_bp, state)
  all_seg[[grp]] <- seg
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  summ[[grp]] <- tibble(group = grp, n = length(nm), segments = nrow(seg), seconds = round(secs, 1))
  cat(sprintf("  [%s] fit %d samples in %.1fs -> %d segments\n", grp, length(nm), secs, nrow(seg)))
}
out <- bind_rows(all_seg) %>% arrange(donor, name, chr, start_bp)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_csv(out, file.path(outdir, "calls_common_schema.csv"))
write_csv(bind_rows(summ), file.path(outdir, "run_summary.csv"))
cat(sprintf("\nWROTE %s/calls_common_schema.csv: %d samples, %d segments across %d groups\n",
            outdir, n_distinct(out$name), nrow(out), length(all_seg)))
