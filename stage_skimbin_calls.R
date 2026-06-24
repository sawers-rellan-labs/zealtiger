#!/usr/bin/env Rscript
# Stage Skim-BIN ancestry calls (BzeaSeq 1 Mb-bin GMM-HMM, SAME skim reads as Skim-RTI)
# into the common call schema used across skim_vs_brbseq.qmd.
#
# Source : data/bzeaseq_ancestry/all_samples_ancestry_segments.bed  (full 1,434-sample set,
#          staged from /Volumes/BZea/bzeaseq/ancestry/all_samples_ancestry_segments.bed)
#          cols: sample, chrom(chr1..chr10), start, end, genotype{REF/HET/ALT}, score, alt_freq
# Output : results/sim_calibration/skimbin_calls/calls_common_schema.csv
#          schema: source, donor, name, chr(int), start_bp, end_bp, state(REF=0/HET=1/ALT=2)
suppressMessages(library(tidyverse))
src <- "data/bzeaseq_ancestry/all_samples_ancestry_segments.bed"
if (!file.exists(src))
  stop("missing ", src, " — stage from /Volumes/BZea/bzeaseq/ancestry/all_samples_ancestry_segments.bed")
bed <- read_tsv(src, show_col_types = FALSE)
out <- bed %>%
  transmute(source = "SkimBIN_1Mb_GMMHMM", donor = NA_character_, name = sample,
            chr = as.integer(sub("chr", "", chrom)), start_bp = start, end_bp = end,
            state = match(genotype, c("REF", "HET", "ALT")) - 1L) %>%
  filter(!is.na(chr), !is.na(state)) %>% arrange(name, chr, start_bp)
dir.create("results/sim_calibration/skimbin_calls", showWarnings = FALSE, recursive = TRUE)
write_csv(out, "results/sim_calibration/skimbin_calls/calls_common_schema.csv")
cat(sprintf("wrote skimbin_calls: %d samples, %d segments | states %s\n",
            n_distinct(out$name), nrow(out), paste(sort(unique(out$state)), collapse = "/")))
