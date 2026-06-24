#!/usr/bin/env Rscript
# Provisional cross-experiment CHECK correspondence by ancestry-call similarity.
#
# Checks carry no pedigree, so they can't be matched across experiments by id.
# Match them by the SHAPE of their ancestry calls instead — the PN4_SID330
# mislabel idea (agent/find_pn4_sid330_match.R) generalised to a similarity matrix
# over every check sample (Skim-RTI r=5 + BrB-RTI r=109).
#
# Skim (50K array sites) and BrB (RNA SNP sites) are on DIFFERENT marker grids, so
# similarity is on a shared 1 Mb BIN grid (segment overlap per bin), not markers.
#
# Per pair, over bins where BOTH define a state:
#   conc3  = 3-state agreement (REF/HET/ALT)
#   jaccH  = Jaccard of HET bins (Purple/B73 are HET-vs-REF, never hom-donor)
# Plus a per-sample genome-wide het fraction (the class discriminant).
suppressMessages({library(tidyverse); library(GenomicRanges)})

sl <- read_csv("data/rtiger_50K/seqlengths.csv", show_col_types = FALSE)
grid <- sl %>% mutate(chr = as.integer(sub("chr", "", chr_label))) %>%
  rowwise() %>% mutate(mid = list(seq(5e5, len, by = 1e6))) %>%
  unnest(mid) %>% transmute(chr, bp = as.integer(mid)) %>% arrange(chr, bp)
gg <- GRanges(grid$chr, IRanges(grid$bp, width = 1))
state_on_bins <- function(seg) seg$state[
  findOverlaps(gg, GRanges(seg$chr, IRanges(seg$start_bp, seg$end_bp)), select = "first")]

read_calls <- function(path, tag) if (!file.exists(path)) NULL else
  read_csv(path, show_col_types = FALSE) %>%
    transmute(name, chr = as.integer(chr), start_bp, end_bp, state, exp = tag)

all_calls <- bind_rows(
  read_calls("data/rtiger_50K/calls_checks_common_schema.csv", "SkimRTI"),
  read_calls("data/brbseq_checks/Purple/calls_common_schema.csv", "BrBRTI"),
  read_calls("data/brbseq_checks/B73/calls_common_schema.csv",    "BrBRTI"))

grp_of <- read_csv("results/sample_correspondence/checks_to_resolve.csv", show_col_types = FALSE) %>%
  transmute(name = prefix, exp = recode(experiment, "Skim"="SkimRTI", "BRB-seq"="BrBRTI"), group)

keys <- all_calls %>% distinct(exp, name)
mat  <- map(seq_len(nrow(keys)), ~ state_on_bins(filter(all_calls, exp==keys$exp[.x], name==keys$name[.x])))
names(mat) <- paste(keys$exp, keys$name, sep = "::")

het_frac <- map_dbl(mat, ~ { v <- .x[!is.na(.x)]; if (!length(v)) NA_real_ else mean(v == 1) })
feat <- tibble(key = names(mat), het_frac = round(het_frac, 3)) %>%
  separate(key, c("exp","name"), sep="::", remove=FALSE) %>%
  left_join(grp_of, by = c("exp","name"))
write_csv(feat, "results/sample_correspondence/check_het_fraction.csv")

conc3 <- function(a,b){ok<-!is.na(a)&!is.na(b);if(!any(ok))return(NA_real_);mean(a[ok]==b[ok])}
jaccH <- function(a,b){ok<-!is.na(a)&!is.na(b);A<-a[ok]==1;B<-b[ok]==1;u<-sum(A|B);if(u==0)return(NA_real_);sum(A&B)/u}

sk <- names(mat)[startsWith(names(mat),"SkimRTI::")]
br <- names(mat)[startsWith(names(mat),"BrBRTI::")]
pairs <- expand_grid(brb=br, skim=sk) %>%
  mutate(brb_name=sub("BrBRTI::","",brb), skim_name=sub("SkimRTI::","",skim),
         conc3 =map2_dbl(brb,skim,~conc3(mat[[.x]],mat[[.y]])),
         jaccH =map2_dbl(brb,skim,~jaccH(mat[[.x]],mat[[.y]]))) %>%
  left_join(grp_of%>%filter(exp=="BrBRTI") %>%transmute(brb_name=name, brb_group=group), by="brb_name") %>%
  left_join(grp_of%>%filter(exp=="SkimRTI")%>%transmute(skim_name=name,skim_group=group),by="skim_name")
write_csv(pairs, "results/sample_correspondence/check_similarity_skimRTI_vs_brbRTI.csv")

cat("=== genome-wide het fraction by experiment x group (class discriminant) ===\n")
feat %>% group_by(exp, group) %>%
  summarise(n=n(), het_mean=round(mean(het_frac),3),
            het_min=round(min(het_frac),3), het_max=round(max(het_frac),3), .groups="drop") %>%
  as.data.frame() %>% print()

cat("\n=== label-vs-ancestry anomalies (het_frac off its group) ===\n")
feat %>% filter((group=="B73" & het_frac>0.15) | (group=="Purple" & het_frac<0.30)) %>%
  arrange(group, desc(het_frac)) %>% as.data.frame() %>% print()

cat("\n=== mean het-Jaccard, BrB x Skim, by group pairing ===\n")
pairs %>% group_by(brb_group, skim_group) %>%
  summarise(jaccH=round(mean(jaccH,na.rm=TRUE),3), conc3=round(mean(conc3,na.rm=TRUE),3), .groups="drop") %>%
  as.data.frame() %>% print()

p <- pairs %>% mutate(brb_lab=paste0(brb_group,":",brb_name), skim_lab=paste0(skim_group,":",skim_name)) %>%
  ggplot(aes(skim_lab, brb_lab, fill=jaccH)) + geom_tile() +
  scale_fill_viridis_c(limits=c(0,1), name="HET\nJaccard") +
  labs(title="Check correspondence — BrB-RTI (r=109) vs Skim-RTI (r=5), 1 Mb bins",
       subtitle="bright = shared HET footprint; Purple is het-rich, B73 recurrent",
       x="Skim-RTI check", y="BrB-RTI check") +
  theme_minimal(base_size=8) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave("results/sample_correspondence/check_similarity_heatmap.png", p, width=11, height=5, dpi=150)
cat("\nwrote check_het_fraction.csv, check_similarity_skimRTI_vs_brbRTI.csv, check_similarity_heatmap.png\n")