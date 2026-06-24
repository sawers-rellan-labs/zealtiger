#!/usr/bin/env Rscript
# Regenerate the coverage-sweep members + their ancestry call sets for
# skim_vs_brbseq.qmd, using the FIXED grid_series (picks track the evenly-spaced
# grid targets, so a dense coverage region no longer piles several steps onto the
# same value). Replaces the stale, externally-precomputed
# results/sim_calibration/{coverage_sweep_members.csv, skim_sweep_calls/, brbseq_sweep_calls_r109/}.
#
# Skim sweep calls are free (subset calls_taxa_r5 + the B73 control from the check fit).
# BRB sweep calls must be fit at r=109 from allelicCounts — this script writes the BRB
# partner list to /tmp/brb_sweep_fit.txt; fit it with fit_rtiger_brbseq.R, then the
# painting reads the regenerated brbseq_sweep_calls_r109.
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

# ---- FIXED grid_series: pick among the k nearest-to-target unused NILs ----
grid_series <- function(df, vary, hold, hold_mean, n = 10, q = c(0.025, 0.975), k = 3) {
  v <- df[[vary]]; h <- df[[hold]]
  gr <- seq(quantile(v, q[1]), quantile(v, q[2]), length.out = n)
  chosen <- integer(0)
  for (g in gr) {
    free <- setdiff(seq_along(v), chosen)
    near <- free[order(abs(v[free] - g))][seq_len(min(k, length(free)))]
    chosen <- c(chosen, near[which.min(abs(h[near] - hold_mean))])
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

# ---- BRB partner list to fit at r=109 (members + B73 control) ----
brb_ids <- unique(c(members$brb_name, BRB_B73))
writeLines(brb_ids, "/tmp/brb_sweep_fit.txt")

cat(sprintf("matched NILs with skim calls: %d | skim mean %.3f× | brb mean %.2f×\n", nrow(cov), ms, mb))
cat("vary_skim picks (skim_cov):\n"); print(round(sort(vary_skim$skim_cov), 3))
cat("vary_brb picks (brb_cov):\n");   print(round(sort(vary_brb$brb_cov), 2))
cat(sprintf("skim_sweep_calls samples: %d (incl B73 %s)\n", n_distinct(skim_sweep$name), SKIM_B73))
cat(sprintf("BRB to fit at r=109: %d samples -> /tmp/brb_sweep_fit.txt\n", length(brb_ids)))
