# SNP50K cohort — provenance & genotype-content analysis

What the "SNP50K" dataset actually is, how it was made, and why its hard
genotypes are unusable at face value (so the ancestry callers run on allele
**counts**, not `GT`). Settles the recurring question of whether the matrix is
imputed. **It is not.**

Source file analysed: `…/rrellan/BZea/bzeaseq/50K/results/joint/cohort.vcf.gz`
(the raw joint-genotyping output; the `…_cohort`/`…_cohort_ref` files are
reshapes of it — see [File inventory](#file-inventory)).

## TL;DR

- **No imputation, at any stage.** `cohort.vcf.gz` is raw `bcftools mpileup`
  + `call -mv`. The `_cohort_ref` file *merges in* the Grzybowski (2023)
  diversity panel as **extra samples**; it never fills genotypes.
- **The hard `GT` matrix is 70% missing** and, among called cells,
  **97% REF / 3% HET / 0.00% ALT-hom** — because ~0.4× coverage means 75% of
  called genotypes rest on a **single read**, and homozygous-ALT is essentially
  uncallable.
- Therefore the ancestry callers consume the per-sample **allele-count** tables
  (`allelic_counts/*_allele_counts.tsv`), not `GT`.
- The "imputed REF positions cause over-fragmentation / only-HET" hypothesis is
  **refuted**: there are no imputed positions; REF cells are genuine single-read
  hom-ref pileup calls; over-fragmentation is count-emission noise, controlled by
  the `err` knob (settled skim config `r=4e-4, err=0.01`), not the transition
  rate. See [15-rigidity-and-scoring.md](15-rigidity-and-scoring.md).

## Provenance chain (from the VCF header)

```
bcftools mpileup -f Zm-B73-REFERENCE-NAM-5.0.fa \
                 -R HQ_BZEA.vcf.gz  -Ou  -b bam.list      # bcftools 1.21
  → bcftools call -mv -Oz -o cohort.vcf.gz                 # Thu Aug 21 2025
```

- **Reference:** B73 NAM v5 (`Zm-B73-REFERENCE-NAM-5.0.fa`).
- **Target panel:** `…/rrellan/BZea/bzeaseq/nilhmm/vcf/HQ_BZEA.vcf.gz` =
  **51,991** high-quality teosinte-vs-B73 SNP sites. `-R` restricts the pileup
  to these — this is "the 50K".
- **`call -mv`** emits variant sites only → **50,995** SNP records
  (1,993 multiallelic) across **1,439** cohort samples.
- **FORMAT:** `GT`, `AD` (allelic depths), `PL`. **INFO:** `DP`, `DP4`, `MQ`,
  `AC`, `AN`, plus mapping/base-quality bias metrics.
- No `beagle`/`impute`/phasing step appears in any header, and no phased file
  exists in the directory. Confirmed against the README and
  `get_cohort_and_reference_vcf.sh`.

## Genotype content (all 1,439 × 50,995 = 73,381,805 cells)

| class            | count       | % of all | % of *called* |
|------------------|------------:|---------:|--------------:|
| REF-hom `0/0`    | 21,168,920  |   28.85  |         97.02 |
| HET              |    648,539  |    0.88  |          2.97 |
| ALT-hom          |      1,001  |    0.00  |          0.00 |
| missing `./.`    | 51,563,345  |   70.27  |            —  |

There is only a "before" — **no imputation to compare against.**

## Why "only HET, no ALT" — it is coverage, not imputed REF

Depth = sum of `AD` per genotype cell:

- mean depth **0.393×** across all cells;
- mean depth **1.33×** among covered cells;
- only **29.7%** of cells have ≥1 read;
- of *called* genotypes, **74.8% rest on a single read**, 19.3% on two, 5.7% on ≥3.

Per-genotype depth × class (cell counts):

| depth | REF        | HET     | ALT | missing     |
|-------|-----------:|--------:|----:|------------:|
| 0     |     40,524 |       1 |   0 | 51,563,345  |
| 1     | 15,868,677 | 448,805 | 355 |           0 |
| 2     |  4,057,642 | 154,596 | 138 |           0 |
| 3–5   |  1,184,345 |  44,361 | 343 |           0 |
| 6–10  |     17,559 |     727 | 163 |           0 |
| >10   |        173 |      49 |   2 |           0 |

A single read can only support REF (or a marginal HET); calling a homozygous-ALT
genotype needs several concordant ALT reads, which almost never occur here —
hence **1,001 ALT-hom cells in 73 million (0.00%)**. The REF cells are not
spurious fills; they are true single-read hom-ref calls, and most of a BC2S3
NIL genome genuinely *is* B73.

## Why the callers use counts, not GT

Because the hard `GT` above is dead (70% missing, ALT unresolvable), the ancestry
callers consume per-sample allele counts:

- **`…/50K/results/allelic_counts/`** — **1,439** tables
  `SAMPLE_allele_counts.tsv`, columns `chr  pos  ref  n_ref  alt  n_alt`,
  built by `get_allelic_counts.sh` (per-sample `bcftools query` of
  `FORMAT/AD{0}` and `AD{1}` over the cohort positions).
- These feed the count / BetaBinomial emission (nilHMM `nnil`, RTIGER), which
  recover the donor footprint that the hard genotypes cannot.

## File inventory (`…/50K/results/joint/`)

| File | What it is |
|---|---|
| `cohort.vcf.gz` | **Raw** joint call (mpileup+call), full INFO/FORMAT (GT+AD+PL). 1,439 × 50,995. The source of everything below. |
| `cohort.vcf.stats.txt` | `bcftools stats` summary (no per-sample counts). |
| `bzea_50K_cohort.vcf.gz` | Simplified: biallelic SNPs, **GT-only**, minimal header. |
| `bzea_50K_ref_panel.vcf.gz` | Grzybowski (2023) diversity panel, subset to cohort positions. |
| `bzea_50K_cohort_ref.vcf.gz` | **Merge** of cohort + ref panel (extra samples at the same sites). *Not imputed.* |
| `bzea_50K_cohort_ref_metadata.csv` / `.md` | Per-sample metadata (cohort vs reference, taxa, checks, donor, pedigree). |
| `allelic_counts/*_allele_counts.tsv` | 1,439 per-sample count tables — the count-path input. |
| `get_cohort_and_reference_vcf.sh` | The reshape/merge script (Steps 1–4). |

## Key facts to carry forward

- **"SNP50K" = 51,991 HQ teosinte-vs-B73 panel sites** pileup-genotyped at
  ~0.4× skim coverage; **not an array, not imputed.**
- **Hard GT is unusable** (70% missing, 0% ALT-hom); the callers run on
  `allelic_counts/` reads.
- **Only-HET** = depth-1 zygosity degeneracy; **over-fragmentation** = emission
  noise (fix with `err`, not `r`). Corroborates the settled skim calibration
  `r=4e-4, err=0.01`.

## Reproduce

```bash
D=/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq/50K/results/joint

# provenance + dimensions
bcftools view -h "$D/cohort.vcf.gz" | grep -iE '^##(bcftools|reference|FORMAT|INFO)'
bcftools query -l "$D/cohort.vcf.gz" | wc -l          # 1439 samples
bcftools index -n "$D/cohort.vcf.gz"                  # 50995 sites

# genotype-class frequencies
bcftools query -f '[%GT\n]' "$D/cohort.vcf.gz" | awk -F'[/|]' '
  {n++; if($1=="."||$2=="."){m++} else if($1==0&&$2==0){r++}
        else if($1==$2){a++} else {h++}}
  END{printf "REF %d  HET %d  ALT %d  miss %d  (n=%d)\n",r,h,a,m,n}'

# depth × class (sum of AD per cell)
bcftools query -f '[%GT %AD\n]' "$D/cohort.vcf.gz" | awk '
  {gt=$1; n=split($2,A,","); dp=0; for(i=1;i<=n;i++) if(A[i]~/^[0-9]+$/) dp+=A[i];
   print gt, dp}'   # bin/tabulate as needed
```

---
*Analysis run 2026-07-02 against `cohort.vcf.gz` (mounted read-only). See also
[data-sources.md](data-sources.md) and
[06-rtiger-benchmark-dataset.md](06-rtiger-benchmark-dataset.md) (the SNP50K
coverage model).*