#!/usr/bin/env Rscript
# Convert MolBreeding 45K array genotypes → RTIGER 6-col input at B73 v5 sites.
#
# Array calls are near-certain genotypes (no reads/pileup), so we encode each
# call as a 0/1/2 dosage of PSEUDO allele-counts at a fixed depth D:
#   RefRef (state 0) → (RefCount=D, AltCount=0)
#   het    (state 1) → (D/2, D/2)
#   AltAlt (state 2) → (0, D)
#   missing/invalid  → (0, 0)
# Polarity from snp.annot.xls: Ref = B73 (recurrent), Alt = donor (teosinte) —
# so state matches the rest of the project (0=recurrent, 1=het, 2=donor).
# Markers are placed at their lifted v5 positions (molbreeding_liftover.R output).
#
# Inputs : data/GSER2026030032P01/SNP/{All.GT.ped, All.GT.map, snp.annot.xls}
#          data/molbreeding_45k_v5_sites.tsv      (run molbreeding_liftover.R first)
#          data/molbreeding_45k_sample_map.tsv    (well_id → PN#_SID#/genotype)
# Outputs: data/molbreeding_45k/counts/<sample>.tsv  (RTIGER 6-col, no header)
#          data/molbreeding_45k/expDesign.csv  (files,name,donor=accession)
#          data/molbreeding_45k/seqlengths.csv (copied v5 lengths)
#
# Usage:  Rscript molbreeding_to_rtiger.R [pseudo_depth]    # default depth 20
#
# To fit:  RTIGER(expDesign[,c("files","name")], outputdir=..., seqlengths=<v5>,
#                  rigidity=3, autotune=FALSE, save.results=FALSE)   # 16 samples, one quick fit

suppressMessages({library(data.table)})

snp_dir   <- "data/GSER2026030032P01/SNP"
lift_path <- "data/molbreeding_45k_v5_sites.tsv"
smap_path <- "data/molbreeding_45k_sample_map.tsv"
seq_src   <- "data/rtiger_50K/seqlengths.csv"      # B73 v5 chr lengths (reuse)
out_root  <- "data/molbreeding_45k"

args <- commandArgs(trailingOnly = TRUE)
D  <- if (length(args) >= 1) as.integer(args[1]) else 20L
Dh <- as.integer(round(D / 2))

dir.create(file.path(out_root, "counts"), recursive = TRUE, showWarnings = FALSE)

# --- marker order (v3) + Ref/Alt polarity + lifted v5 position ---------------
map <- fread(snp_dir |> file.path("All.GT.map"), header = FALSE)
setnames(map, c("chr", "marker", "cm", "bp"))
nmark <- nrow(map)

annot <- fread(file.path(snp_dir, "snp.annot.xls"))          # Chrom Start End Ref Alt ...
amkey <- paste0(annot$Chrom, "_", annot$Start)
mkey  <- paste0(map$chr, "_", map$bp)
ai    <- match(mkey, amkey)
map[, ref := annot$Ref[ai]]
map[, alt := annot$Alt[ai]]

lift <- fread(lift_path)                                     # marker chr_v5 pos_v5 ...
li   <- match(map$marker, lift$marker)
map[, v5chr := lift$chr_v5[li]]
map[, v5pos := lift$pos_v5[li]]

# usable markers: lifted, polarisable (ref/alt known and distinct)
ok <- !is.na(map$v5pos) & !is.na(map$ref) & !is.na(map$alt) & map$ref != map$alt
cat(sprintf("markers: %d total | %d lifted | %d usable (lifted & ref!=alt)\n",
            nmark, sum(!is.na(map$v5pos)), sum(ok)))

# stable v5 sort order among usable markers
ord <- order(as.integer(map$v5chr[ok]), map$v5pos[ok])
mk  <- map[ok][ord]                                          # usable markers, v5-sorted
ref_v <- mk$ref; alt_v <- mk$alt
keep_idx <- which(ok)[ord]                                   # original marker indices

# --- ped: each marker is ONE tab-delimited field "allele1 allele2" -----------
# (this MolBreeding .ped is tab-separated with a space between the two alleles
#  inside each cell, i.e. 6 + nmark columns — not the 2-col-per-marker layout)
ped <- fread(file.path(snp_dir, "All.GT.ped"), header = FALSE, sep = "\t")
iids <- as.character(ped[[2]])
gm   <- as.matrix(ped[, 7:(6 + nmark)])                     # 16 × nmark, cells like "A G"

# --- sample naming: well_id → PN#_SID#, accession from genotype prefix -------
smap <- fread(smap_path)
si   <- match(iids, smap$well_id)
sample_name <- ifelse(is.na(si), iids, smap$sample_id[si])
genotype    <- ifelse(is.na(si), iids, smap$genotype[si])
accession   <- sub("_.*$", "", genotype)                    # "Zx.0340_…" → "Zx.0340"; "B73-bulk" stays

# --- write one RTIGER 6-col file per sample ----------------------------------
ed <- data.table(files = character(), name = character(), donor = character())
dosage <- numeric(length(iids))
for (i in seq_along(iids)) {
  sp <- tstrsplit(gm[i, keep_idx], " ", fixed = TRUE)       # split "A G" → alleles
  a1 <- sp[[1]]; a2 <- sp[[2]]
  cref <- (a1 == ref_v) + (a2 == ref_v)
  calt <- (a1 == alt_v) + (a2 == alt_v)
  valid <- (cref + calt) == 2                                # both alleles ∈ {ref,alt}
  state <- ifelse(!valid, NA_integer_, ifelse(cref == 2, 0L, ifelse(calt == 2, 2L, 1L)))
  RefCount <- ifelse(is.na(state), 0L, ifelse(state == 0L, D, ifelse(state == 1L, Dh, 0L)))
  AltCount <- ifelse(is.na(state), 0L, ifelse(state == 2L, D, ifelse(state == 1L, Dh, 0L)))

  out <- data.table(chr = paste0("chr", mk$v5chr), pos = mk$v5pos,   # match seqlengths chr1..chr10
                    RefBase = ref_v, RefCount = RefCount,
                    AltBase = alt_v, AltCount = AltCount)
  f <- normalizePath(file.path(out_root, "counts", paste0(sample_name[i], ".tsv")),
                     mustWork = FALSE)
  fwrite(out, f, sep = "\t", col.names = FALSE)
  ed <- rbind(ed, data.table(files = f, name = sample_name[i], donor = accession[i]))

  cov <- !is.na(state)
  dosage[i] <- sum(state[cov]) / (2 * sum(cov))              # donor dosage
  cat(sprintf("  %-13s %-9s n_called=%d  dosage=%.3f\n",
              sample_name[i], accession[i], sum(cov), dosage[i]))
}

fwrite(ed, file.path(out_root, "expDesign.csv"))
file.copy(seq_src, file.path(out_root, "seqlengths.csv"), overwrite = TRUE)

cat(sprintf("\nWrote %d sample files + expDesign.csv (depth D=%d) → %s\n",
            nrow(ed), D, out_root))
cat(sprintf("mean donor dosage across samples = %.3f\n", mean(dosage)))
cat("Fit:  ed <- read.csv('data/molbreeding_45k/expDesign.csv'); sl <- read.csv('data/molbreeding_45k/seqlengths.csv')\n")
cat("      RTIGER(ed[,c('files','name')], outputdir='data/molbreeding_45k/fit', \n")
cat("             seqlengths=setNames(sl$len, sl$chr_label), rigidity=3, autotune=FALSE, save.results=FALSE)\n")
