# External data sources & resources

The zealtiger repo tracks code + rendered notebooks only; the actual inputs/outputs
(`/data/`, `/results/`) and the working area (`/agent/`) are git-ignored. This page
catalogs the **external resources** (outside the repo, or under git-ignored `data/`) that
the analyses depend on, so they're retrievable from the documentation.

Most live on the NCSU **`rsstu`** mount: `/Volumes/rsstu/users/r/rrellan/` (paths below
abbreviate this as `…/rrellan/`). Others are in Google Drive or local `~/ref`.

## Sequencing data & panels

| Resource | What it is | Used by |
|---|---|---|
| `…/rrellan/BZea/bzeaseq/50K/allelic_counts50K.tsv` | SNP50K skim allele counts, 1,439 samples | Source 2 (RTIGER skim) |
| `…/rrellan/BZea/bzeaseq/wideseq_ref/` (`wideseq_all.pos`, `wideseq_chr*.vcf.gz`) | ~27.6M teosinte-vs-B73 variant panel (Schnable 2023), B73 v5, ~12,950 sites/Mb | Source 3 (wideseq); Source 4 site list |
| `…/rrellan/BZea/bzeaseq/ancestry/all_samples_bin_genotypes.tsv` (309 MB) | real wideseq per-1Mb-bin calls + counts, 1,434 samples (`Kgmm_HMM` etc.) | Source 3 real calls; λ fit; ergodicity check |
| `…/rrellan/BZea/bzeaseq/ancestry/PN#_SID#_bin_genotypes.tsv` (1,441 files) | per-sample wideseq bin tables | per-sample inspection |
| `…/rrellan/BZea/bzeaseq/WGSmetrics_summary.tsv` | Picard per-sample genome-wide coverage, 1,437 samples (λ≈0.39) | coverage-distribution fit |
| `…/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/` | BRBseq: `metadata.txt`, `clean_reads/` (379 SE fastqs), `alignments/` (STAR BAMs) | Source 4 (BRBseq) |
| `data/GSER2026030032P01/` (repo `data/`, git-ignored) | MolBreeding GBTS: `mSNP/` (264,553 sites) + `SNP/` (44,935), real `ref_depth`/`alt_depth` ~200×, `report.html` | Source 5 |

## References & reference genomes

| Resource | What it is | Used by |
|---|---|---|
| `…/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa` (+ `.fai`, `.dict`), `NAM5_CHR` STAR index, `…56.gtf`, `zea_mays.vcf.gz` | B73 NAM v5 reference set | RTIGER seqlengths; rnagt/STAR |
| `~/ref/zea/chain_files/` (`B73_RefGen_v4_to_…NAM-5.0.chain`, `AGPv3_to_B73_RefGen_v4.chain`) + `~/ref/zea/GCA_000005005.6_B73_RefGen_v4_genomic.fna` | liftover chains (v4→v5, v3→v4) + B73 v4 genome | MolBreeding marker liftover |

## Method code (repos)

| Resource | What it is | Used by |
|---|---|---|
| `~/Library/CloudStorage/GoogleDrive-…/My Drive/repos/BzeaSeq/` | wideseq method: `docs/call_ancestry.Rmd` (3-state GMM+HMM), `docs/missing_data.Rmd` (coverage model), `R/ancestry_processing.R` | Source 3 method-of-record |
| `~/Library/CloudStorage/GoogleDrive-…/My Drive/repos/rnagt/` (GitHub `sawers-rellan-labs/rnagt`) | Nextflow GATK RNA-seq genotyping (LSF/Hazel) | Source 4 pipeline |
| `…/rrellan/tlaloc/nilhifi/` | chromosome-painting / 2-state Gaussian-HMM (`nextflow/bin/chromosome_painting.R`, `chrpainting_report_arpae.Rmd`) | painting reference (mean-0 REF emission) |
| `…/rrellan/DOE_CAREER/inv4m/nilhifimi21/` | HiFi genome **assembly** project (hifiasm/ragtag/anchorwave) — not the painter | context |
| RTIGER fork `faustovrz/RTIGER@optimize-julia-core` (installed R pkg); Julia via `~/.juliaup/bin/julia` | optimized RTIGER (streaming M-step, autotune, progress log) | all RTIGER fits |

## Compute & public refs

| Resource | What it is |
|---|---|
| **Hazel** (NC State HPC, LSF); conda env `/share/maize/frodrig4/conda/env/rnagt` | cluster for rnagt (Source 4) |
| rpubs.com/faustovrz/1306822 ("BZea WideSeq") | rendered wideseq method |
| rpubs.com/faustovrz/1337797 | rendered ancestry analysis; consumes the **GATK table** (GATK `CollectAllelicCounts` tidy/long readcounts: `SAMPLE CONTIG POSITION REF_NUCLEOTIDE ALT_NUCLEOTIDE REF_COUNT ALT_COUNT TOTAL_COUNT`) as input |
| GitHub `sawers-rellan-labs/BzeaSeq` `docs/getting_teosinte_variants_from_Schnable2023.md` | teosinte-variant site filtering |
| bioRxiv 2026.02.20.707111v1 | wideseq-style method paper |

> Working handover (task status + plan of action) is in the git-ignored
> `agent/HANDOVER.md`.
