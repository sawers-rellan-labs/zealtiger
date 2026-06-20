#!/usr/bin/env Rscript
# Fit RTIGER on the BRB-seq (3' bulk RNA-seq) data using REAL allele counts.
#
# BRB-seq is low-coverage RNA sequencing: GATK CollectAllelicCounts was run over
# the brbseq joint-called SNP sites, giving per-sample REF/ALT read counts with
# REF = B73 NAM v5 reference allele and ALT = the observed alternate. For a
# teosinte-introgression NIL on a B73 background, ALT reads at a segregating site
# are donor (teosinte) reads, so the counts are already in the RTIGER polarity
# (REF = recurrent/B73, ALT = donor/teosinte) like the SNP50K skim (Source 2) and
# the MolBreeding GBTS (Source 5) — no wideseq panel filter is applied here, the
# brbseq sites are used as-is.
#
# Mirrors fit_rtiger_molbreeding.R (real-count pipeline, MLE E-step, ~20x cap) but
# reads the per-sample *.allelicCounts.tsv files instead of one stacked GATK table.
#
# Sample naming: the brbseq files (PN4_SID290 ...) follow the wideseq/MolBreeding
# plate scheme, so donor accessions come from data/molbreeding_45k_sample_map.tsv.
# By default every sample_map entry that has an allelicCounts file is fitted (the
# samples that overlap the MolBreeding-vs-skim comparison); pass names to override.
#
# DEFAULTS:  --r=8  --estimator=mle  --cap=20   (cap no-ops; brbseq is ~2-3x)
#
# Input : <AC dir>/<name>.allelicCounts.tsv   (CONTIG POSITION REF_COUNT ALT_COUNT
#                                               REF_NUCLEOTIDE ALT_NUCLEOTIDE; SAM @ header)
#         data/molbreeding_45k_sample_map.tsv, data/rtiger_50K/seqlengths.csv
# Output: data/brbseq/{counts/,expDesign.csv,fit/,calls_common_schema.csv}
#           source = BRBseq_r<r>_<METHOD>[_cap<N>]
#
# Usage:  Rscript fit_rtiger_brbseq.R [name ...] [--ac=DIR] [--estimator=mle|mom] [--cap=N] [--r=N]

suppressMessages({
  library(data.table); library(RTIGER); library(readr)
  library(tidyverse); library(GenomicRanges); library(S4Vectors)
})
source("R/06_score_rtiger.R")   # segments_from_rtiger_object()
set.seed(42)                    # reproducible binomial thinning

## ---- args ------------------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  hit <- grepl(paste0("^", flag, "="), argv)
  if (any(hit)) sub(paste0("^", flag, "="), "", argv[hit][1]) else default
}
r         <- as.integer(getopt("--r", "8"))
estimator <- tolower(getopt("--estimator", "mle"))
cap_arg   <- as.integer(getopt("--cap", "20"))
user_cap  <- if (is.na(cap_arg) || cap_arg <= 0) NA_integer_ else cap_arg
ac_dir    <- getopt("--ac",
  "/Volumes/rsstu/users/r/rrellan/BZea/brbseq/run_full/results/allelic_counts")
out_root  <- getopt("--outdir", "data/brbseq")
stopifnot(estimator %in% c("mle", "mom"))
names_arg <- argv[!grepl("^--", argv) & nzchar(argv)]

## ---- Julia -----------------------------------------------------------------
JULIA_HOME <- Sys.getenv("JULIA_HOME", "")
if (!nzchar(JULIA_HOME) || !file.exists(file.path(JULIA_HOME, "julia"))) {
  launcher <- path.expand("~/.juliaup/bin/julia")
  if (file.exists(launcher))
    JULIA_HOME <- system(paste(shQuote(launcher), "--startup-file=no -e",
                               shQuote("print(Sys.BINDIR)")), intern = TRUE)
}
setupJulia(JULIA_HOME = JULIA_HOME); sourceJulia()

.mstep <- "mle"
set_mstep <- function(method) {
  if (method == "mom" && .mstep != "mom") {
    if (!file.exists("mom_mstep.jl")) stop("mom_mstep.jl not found (required for --estimator=mom)")
    JuliaCall::julia_source(normalizePath("mom_mstep.jl")); .mstep <<- "mom"
  } else if (method == "mle" && .mstep != "mle") {
    sourceJulia(); .mstep <<- "mle"
  }
}

## ---- cap helper ------------------------------------------------------------
cap_counts <- function(ref, alt, cap) {                # binomial down-thinning to ~cap x
  tot <- ref + alt
  p   <- pmin(1, cap / pmax(tot, 1L))
  list(ref = rbinom(length(ref), ref, p), alt = rbinom(length(alt), alt, p))
}

## ---- shared inputs ---------------------------------------------------------
smap <- fread("data/molbreeding_45k_sample_map.tsv")
sl   <- read_csv("data/rtiger_50K/seqlengths.csv", show_col_types = FALSE)
seqv <- setNames(sl$len, sl$chr_label)
chr_levels <- sl$chr_label

donor_of <- function(nm) {                              # name -> donor accession (Zx.0340 etc.)
  g <- smap$genotype[match(nm, smap$sample_id)]
  ifelse(is.na(g), NA_character_, sub("_.*$", "", g))
}

# samples to fit: explicit names, else every sample_map entry with an AC file
if (length(names_arg)) {
  samples <- names_arg
} else {
  samples <- smap$sample_id[file.exists(file.path(ac_dir, paste0(smap$sample_id,
                                                                  ".allelicCounts.tsv")))]
}
stopifnot(length(samples) > 0)
cat(sprintf("BRB-seq samples to fit: %d  [%s]\n", length(samples),
            paste(samples, collapse = ", ")))

cnt_dir <- file.path(out_root, "counts"); fit_dir <- file.path(out_root, "fit")
dir.create(cnt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

## ---- per-sample RTIGER count files -----------------------------------------
# CollectAllelicCounts tsv: SAM-style @ header, then
#   CONTIG POSITION REF_COUNT ALT_COUNT REF_NUCLEOTIDE ALT_NUCLEOTIDE
read_ac <- function(path) {
  dt <- fread(path, skip = "CONTIG", sep = "\t")        # skip @HD/@SQ header block
  dt <- dt[CONTIG %in% chr_levels]                      # chr1..chr10, drop scaffolds/organelles
  setorderv(dt, c("CONTIG", "POSITION"))
  dt[, CONTIG := factor(CONTIG, levels = chr_levels)]
  setorder(dt, CONTIG, POSITION)
  dt[, CONTIG := as.character(CONTIG)]
  dt[]
}

ed <- data.table(files = character(), name = character(), donor = character())
depth_tab <- data.table()
for (nm in samples) {
  ac <- read_ac(file.path(ac_dir, paste0(nm, ".allelicCounts.tsv")))
  ref <- ac$REF_COUNT; alt <- ac$ALT_COUNT
  md  <- mean(ref + alt)
  if (!is.na(user_cap) && md > user_cap) {              # cap only bites above target depth
    th <- cap_counts(ref, alt, user_cap); ref <- th$ref; alt <- th$alt
  }
  six <- data.table(chr = ac$CONTIG, pos = ac$POSITION,
                    RefBase = ac$REF_NUCLEOTIDE, RefCount = ref,
                    AltBase = ac$ALT_NUCLEOTIDE, AltCount = alt)
  f <- normalizePath(file.path(cnt_dir, paste0(nm, ".tsv")), mustWork = FALSE)
  fwrite(six, f, sep = "\t", col.names = FALSE)
  ed <- rbind(ed, data.table(files = f, name = nm, donor = donor_of(nm)))
  depth_tab <- rbind(depth_tab, data.table(
    name = nm, sites = nrow(ac), covered = sum(ref + alt > 0), mean_depth = md))
}
fwrite(ed, file.path(out_root, "expDesign.csv"))
cat("Per-sample coverage:\n"); print(depth_tab)
cat(sprintf("estimator %s | cap %s | rigidity %d\n", toupper(estimator),
            if (is.na(user_cap)) "off" else paste0(user_cap, "x (no-op when below)"), r))
set_mstep(estimator)

## ---- fit (one joint EM over the samples) -----------------------------------
prog <- file.path(fit_dir, "fit_progress.log"); if (file.exists(prog)) file.remove(prog)
t0 <- Sys.time()
res <- RTIGER(expDesign = as.data.frame(ed[, .(files, name)]),
              outputdir = fit_dir, seqlengths = seqv, rigidity = r,
              autotune = FALSE, progress_log = prog, save.results = FALSE,
              threads = 1L, verbose = TRUE)
saveRDS(res, file.path(fit_dir, "rtiger_result.rds"))
cat(sprintf("Fit done in %.2f min (%d samples, M-step=%s).\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            length(res@Viterbi), toupper(estimator)))

## ---- common schema ---------------------------------------------------------
tag <- sprintf("BRBseq_r%d_%s%s", r, toupper(estimator),
               if (!is.na(user_cap)) "" else "")
seg <- segments_from_rtiger_object(res)
seg$donor  <- ed$donor[match(seg$name, ed$name)]
seg$source <- tag
seg <- seg |> select(source, donor, name, chr, start_bp, end_bp, state) |>
  arrange(name, chr, start_bp)
out_csv <- file.path(out_root, "calls_common_schema.csv"); write_csv(seg, out_csv)

l <- seg$end_bp - seg$start_bp + 1
dose <- seg |> mutate(len = l) |> group_by(name) |>
  summarise(d = sum(len * state / 2) / sum(len), .groups = "drop")
cat(sprintf("Wrote %s [%s]: %d samples, %d segments, mean dosage %.3f\n",
            out_csv, tag, n_distinct(seg$name), nrow(seg), mean(dose$d)))
print(dose)
cat("\nALL DONE.\n")