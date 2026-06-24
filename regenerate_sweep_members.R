#!/usr/bin/env Rscript
# Regenerate the coverage-sweep members + their ancestry call sets for
# skim_vs_brbseq.qmd, using the FIXED grid_series (picks track the evenly-spaced
# grid targets, so a dense coverage region no longer piles several steps onto the
# same value). Replaces the stale, externally-precomputed
# results/sim_calibration/{coverage_sweep_members.csv, skim_sweep_calls/, brbseq_sweep_calls_r109/}.
#
# Both platforms' calls are SUBSET from per-(teosinte-)species fits — never an ad-hoc cohort:
#   skim -> calls_taxa_r5 + the B73 control from the check fit.
#   BrB  -> brbseq_taxa_r109 (per species, B73/Purple as checks; produced once by
#           fit_rtiger_brbseq_by_taxa.R over the cached brbseq_ks counts). This script
#           subsets the sweep members AND the 8-sample MolB-overlap pair set from it.
# The sweep BRB sample list is recorded in agent/brb_sweep_fit.txt (survives reboots, /tmp doesn't).
suppressMessages({library(readr); library(dplyr)})
chop <- function(x) sub("\\.B$", "", x)
SKIM_B73 <- "PN10_SID893"; BRB_B73 <- "PN1_SID15"

# ---- cov: pedigree-matched NIL coverage on both platforms (as in the doc) ----
brbmeta <- read_csv("data/BRBseq_metadata.csv", show_col_types = FALSE) %>%
  transmute(sid = sample_id, nil = chop(genotype))
skimped <- read_csv("data/skim_sample_pedigree.csv", show_col_types = FALSE)
brb_cov <- read_csv("results/sim_calibration/brbseq_coverage.csv", show_col_types = FALSE) %>%
  mutate(sid = as.integer(sub(".*_SID", "", name)),
         nil = brbmeta$nil[match(sid, brbmeta$sid)]) %>%
  filter(!is.na(nil)) %>% transmute(nil, brb_name = name, brb_cov = lambda)
skim_cov <- read_csv("results/sim_calibration/real_skim_lambda.csv", show_col_types = FALSE) %>%
  mutate(nil = skimped$pedigree[match(name, skimped$sample)]) %>% filter(!is.na(nil)) %>%
  group_by(nil) %>% summarise(skim_name = first(name), skim_cov = first(lambda), .groups = "drop")
cov <- inner_join(brb_cov, skim_cov, by = "nil")
# restrict to NILs that actually have a Skim-RTI call (so the skim lane is never blank)
skim_calls <- read_csv("data/rtiger_50K/calls_taxa_r5.csv", show_col_types = FALSE)
cov <- cov %>% filter(skim_name %in% skim_calls$name)
ms <- mean(cov$skim_cov); mb <- mean(cov$brb_cov)

# ---- grid_series: nearest NIL to (grid_target, hold_mean) in screen space ----
# Must stay identical to the grid-sample chunk in skim_vs_brbseq.qmd.
W <- c(skim_cov = 1, brb_cov = 0.1)   # per-unit screen weight (skim:BRB ≈ 10:1)
grid_series <- function(df, vary, hold, hold_mean, n = 10, q = c(0.025, 0.975)) {
  v <- df[[vary]]; h <- df[[hold]]; wv <- W[[vary]]; wh <- W[[hold]]
  gr <- seq(quantile(v, q[1]), quantile(v, q[2]), length.out = n)
  chosen <- integer(0)
  for (g in gr) {
    free <- setdiff(seq_along(v), chosen)
    D <- (wv * (v[free] - g))^2 + (wh * (h[free] - hold_mean))^2
    chosen <- c(chosen, free[which.min(D)])
  }
  df[chosen, , drop = FALSE] %>% mutate(grid_target = gr)
}
vary_brb  <- grid_series(cov, "brb_cov",  "skim_cov", ms)
vary_skim <- grid_series(cov, "skim_cov", "brb_cov", mb)

# ---- members file ----
members <- bind_rows(
  vary_brb  %>% mutate(sweep = "vary_brb",  vary_cov = brb_cov),
  vary_skim %>% mutate(sweep = "vary_skim", vary_cov = skim_cov)) %>%
  select(sweep, nil, skim_name, brb_name, skim_cov, brb_cov, vary_cov, grid_target)
write_csv(members, "results/sim_calibration/coverage_sweep_members.csv")

# ---- skim sweep calls: subset calls_taxa_r5 for members + B73 control ----
skim_ids <- unique(c(members$skim_name, SKIM_B73))
chk <- read_csv("data/rtiger_50K/calls_checks_common_schema.csv", show_col_types = FALSE)
skim_sweep <- bind_rows(
  skim_calls %>% filter(name %in% members$skim_name),
  chk %>% filter(name == SKIM_B73) %>% select(any_of(names(skim_calls))))
dir.create("results/sim_calibration/skim_sweep_calls", showWarnings = FALSE, recursive = TRUE)
write_csv(skim_sweep, "results/sim_calibration/skim_sweep_calls/calls_common_schema.csv")

# ---- BRB sweep calls: subset the per-species fit for members + B73 control ----
TAXA_CALLS <- "results/sim_calibration/brbseq_taxa_r109/calls_common_schema.csv"
if (!file.exists(TAXA_CALLS))
  stop("missing ", TAXA_CALLS, " — run fit_rtiger_brbseq_by_taxa.R first.")
brb_taxa  <- read_csv(TAXA_CALLS, show_col_types = FALSE)
brb_ids   <- unique(c(members$brb_name, BRB_B73))
brb_sweep <- brb_taxa %>% filter(name %in% brb_ids)
miss_sw <- setdiff(brb_ids, unique(brb_sweep$name))
if (length(miss_sw)) warning("sweep BRB names absent from per-species fit: ", paste(miss_sw, collapse = ", "))
dir.create("results/sim_calibration/brbseq_sweep_calls_r109", showWarnings = FALSE, recursive = TRUE)
write_csv(brb_sweep, "results/sim_calibration/brbseq_sweep_calls_r109/calls_common_schema.csv")
writeLines(brb_ids, "agent/brb_sweep_fit.txt")   # record of the sweep BRB samples (was /tmp)

# ---- BRB pair calls: re-source the doc's 8-sample MolB-overlap pairs from the same fit ----
PAIR_SAMPLES <- c("PN4_SID290", "PN4_SID291", "PN4_SID293", "PN4_SID309",
                  "PN4_SID320", "PN4_SID326", "PN4_SID330", "PN4_SID359")
brb_pairs <- brb_taxa %>% filter(name %in% PAIR_SAMPLES)
write_csv(brb_pairs, "results/sim_calibration/brbseq_calls_r109/calls_common_schema.csv")

cat(sprintf("matched NILs with skim calls: %d | skim mean %.3f× | brb mean %.2f×\n", nrow(cov), ms, mb))
cat("vary_skim picks (skim_cov):\n"); print(round(sort(vary_skim$skim_cov), 3))
cat("vary_brb picks (brb_cov):\n");   print(round(sort(vary_brb$brb_cov), 2))
cat(sprintf("skim_sweep_calls samples: %d (incl B73 %s)\n", n_distinct(skim_sweep$name), SKIM_B73))
cat(sprintf("brbseq_sweep_calls_r109 samples: %d (incl B73 %s); list -> agent/brb_sweep_fit.txt\n",
            n_distinct(brb_sweep$name), BRB_B73))
cat(sprintf("brbseq_calls_r109 (pairs) re-sourced: %d samples\n", n_distinct(brb_pairs$name)))
