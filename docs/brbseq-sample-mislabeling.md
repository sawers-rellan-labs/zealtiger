# BRB-seq sample mislabeling — error report & pipeline fix

**Status:** confirmed by genotype. **Scope:** the BRB-seq `run_full` allelic-count
files (`…/brbseq/run_full/results/allelic_counts/PN<plate>_SID<n>.allelicCounts.tsv`).
**Analysis:** `brbseq_sample_identity.qmd` → `docs/brbseq_sample_identity.html`.
**Corrected map:** `brbseq_corrected_labels.csv`.

## Summary of the error

The BRB-seq outputs are named with the **MolBreeding 45K plate scheme**
(`PN<plate>_SID<n>`, resolved through `data/molbreeding_45k_sample_map.tsv`). That
scheme does **not** identify the sample that was physically in each BRB-seq well.
The authoritative identity for a BRB-seq well is the **pedigree string** in
`data/BRBseq_metadata.csv` (`sample_id → genotype/accession`, with `plate`,
`plate_pos`, and the i7/i5F/i5R barcodes). The two disagree for **every** sample.

The `PN#_SID#`/plate index is experiment-local: the same `SID` number maps to a
different NIL in the MolBreeding 45K map than in the BRB-seq metadata. Only the
**pedigree string** is a stable cross-experiment identifier.

## Evidence (genotype-based)

For each BRB sample we matched its RTIGER donor pattern (donor Jaccard on the 9,157-
site wideseq grid) against the full SNP50K skim pool (1,393 samples) and the
MolBreeding set, and compared the recovered identity to (a) the `PN#_SID#` plate
label and (b) the `BRBseq_metadata` accession:

| BRB file | plate-label acc (wrong) | BRBseq_metadata pedigree | genotype match (acc, J) | metadata = genotype |
|---|---|---|---|---|
| PN4_SID290 | Zx.0340 | Zx.0560_P4_P3_P3.1.1.1 | Zl.0040 (0.59) | ✗ (low-J, RNA noise) |
| PN4_SID291 | Zx.0340 | Zx.0570_P2_P1_P2.1.1.1 | Zx.0570 (0.73) | ✓ |
| PN4_SID293 | Zx.0340 | Zd.0030_P2_P3_P2.1.1.1 | Zd.0030 (0.93) | ✓ |
| PN4_SID309 | Zx.0150 | Zx.0370_P2_P5_P1.1.1.1 | Zx.0040 (0.73) | ✗ (low-J, RNA noise) |
| PN4_SID320 | Zx.0030 | Zx.0160_P2_P3_P2.2.1.1 | Zx.0160 (0.98) | ✓ |
| PN4_SID326 | **B73-bulk** | Zx.0020_P3_P3_P1.2.1.1 | Zx.0020 (0.91) | ✓ |
| PN4_SID330 | Zx.0030 | Zx.0580_P2_P3_P3.1.1.1 | Zx.0580 (0.59) | ✓ |
| PN4_SID359 | Zx.0560 | Zd.0010_P3_P2_P3.1.1.1 | Zd.0010 (0.81) | ✓ |

- **plate-label accession = genotype-recovered identity: 0 / 8.**
- **`BRBseq_metadata` accession = genotype-recovered identity: 6 / 8** (the two misses
  are the lowest-Jaccard matches — BRB-seq RNA sparsity, not a metadata error).
- `PN4_SID326`, labelled **"B73-bulk"**, is genotypically **Zx.0020** (teosinte) — it
  is not a recurrent-parent control. This is why that "check" never read ≈ 0 donor.

Conclusion: the sequencing and genotype calls are fine; the **sample→name mapping**
is wrong. The pedigree string in `BRBseq_metadata.csv` is the correct identifier.

## Quick fix (local, interim)

Re-key the existing call files with `brbseq_corrected_labels.csv`
(`brb_file → metadata_pedigree / metadata_acc`). This corrects the labels without
re-running anything. It does **not** make the BRB samples comparable per-NIL to the
MolBreeding/skim PN4 set — under correct labels they are *different* NILs, so any
"same-NIL on 3 platforms" comparison built on the old labels is invalid.

## True fix — for the agent running the BRB-seq pipeline on Hazel (HPC)

Root cause to fix in the pipeline (the `rnagt` Nextflow GATK RNA-seq genotyping
workflow, LSF/Hazel; `~/…/repos/rnagt`):

1. **Use `BRBseq_metadata.csv` as the single source of truth for sample identity.**
   Build the run sample sheet by joining demultiplexed reads → well by the
   **i7/i5F/i5R barcodes + `plate`/`plate_pos`** in `BRBseq_metadata.csv`, and name
   every output by the **pedigree string** (`genotype`), not by any `PN#_SID#`
   plate index borrowed from another experiment. Do **not** route names through
   `molbreeding_45k_sample_map.tsv` — that map belongs to the MolBreeding 45K run
   and its `SID` numbering is unrelated to the BRB-seq plate.
2. **Carry the pedigree string end-to-end.** Emit `…/allelic_counts/<pedigree>.allelicCounts.tsv`
   (or keep `SID` only as a within-run well tag, e.g. `<pedigree>__well<plate_pos>`),
   so downstream joins key on the stable identifier.
3. **Add a genotyping-identity QC gate** to the pipeline (cheap, catches this class
   of bug): after calling, for each sample compute donor-marker concordance against
   the matching DNA genotype (SNP50K skim / MolBreeding) for the **same pedigree**;
   flag any sample whose own-label concordance is not the top match. A true B73/
   recurrent control must read ≈ 0 donor — assert it.
4. **Re-run `run_full`** with the corrected sample sheet and regenerate
   `results/vcf` + `results/allelic_counts`. After the rerun, the cross-platform
   comparison (skim / MolBreeding / BRB-seq) can be rebuilt on genuinely shared NILs.

### Suggested prompt for the pipeline agent

> The BRB-seq `run_full` outputs were named with the MolBreeding 45K plate scheme
> (`PN#_SID#` via `molbreeding_45k_sample_map.tsv`), which does not identify the
> samples actually sequenced. Rebuild the sample sheet from `BRBseq_metadata.csv`
> only — map demultiplexed reads to wells by the i7/i5F/i5R barcodes + plate/plate_pos
> and name every output by the `genotype` pedigree string. Add a post-calling
> identity-QC step that checks each sample's donor pattern against the same-pedigree
> DNA genotype and asserts the B73 control reads ≈ 0 donor. Re-run `run_full` and
> regenerate the VCF and allelic counts. Reference: `docs/brbseq-sample-mislabeling.md`,
> `brbseq_corrected_labels.csv`, `docs/brbseq_sample_identity.html`.
