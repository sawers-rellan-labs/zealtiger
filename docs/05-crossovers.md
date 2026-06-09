# Crossovers per genome (for RTIGER autotune)

**Code:** `n_crossovers_hap` + `crossover_summary` in `R/02_simulate.R` /
`R/04_summarize.R`, driven by `run_all.R --sim`.
**Output:** `results/nil_crossovers.csv`, `results/crossovers_per_chr.csv`,
crossover fields in `results/run_metadata.json`.

## What is counted

RTIGER detects crossovers as **donor↔recurrent allele transitions** along the
resequenced genome. With two founders, every mosaic segment boundary is such a
transition, so per haplotype a CO = `sum(diff(alleles) != 0)`; per diploid NIL we
sum both haplotypes over all chromosomes. This is the *detectable* CO count
RTIGER's `crossovers_per_megabase` autotune input expects — **not** the total
meiotic CO count (crossovers inside the homozygous-recurrent background, ~87% of
the genome, leave no transition and are invisible).

## Results (continuous mosaic, m = 10)

| statistic | value |
|---|---|
| Mean detectable COs / diploid NIL | **33.3** (sd 9.0, median 33, range 8–62) |
| Mean CO / Mb | 0.0156 |
| **Max CO / Mb across NILs** | **0.0291** |
| Genome length used | 2,129.6 Mb (map-spanned) |

Per chromosome (mean / max CO, and the conservative max CO/Mb):

| chr | mean CO | max CO/Mb | chr | mean CO | max CO/Mb |
|---|---|---|---|---|---|
| 1 | 5.25 | 0.0714 | 6 | 2.64 | 0.0665 |
| 2 | 3.45 | 0.0657 | 7 | 2.93 | 0.0754 |
| 3 | 3.83 | 0.0589 | 8 | 2.96 | 0.0768 |
| 4 | 3.49 | 0.0560 | 9 | 3.04 | 0.0737 |
| 5 | 3.22 | 0.0530 | 10 | 2.48 | 0.0788 |

No-interference (m=0) is essentially the same mean (33.4) with a fatter tail
(max CO/Mb 0.0362) — interference respaces COs, it doesn't change the mean count.

## Using it in RTIGER

RTIGER wants the **highest** CO/Mb across samples (higher = more conservative FP
estimate → higher chosen rigidity). Pass the per-chromosome vector if giving a
per-chr rate, or the genome max (0.029); the m=0 value (0.036) is the
conservative ceiling.

> **Scale note.** Per-chromosome max rates (0.053–0.079) exceed the genome-wide
> max (0.029) because the genome figure averages busy and quiet chromosomes
> within a line. Use per-chr rates with per-chr `seqlengths`; don't mix scales.

These are *continuous-mosaic* COs. At the 50K-marker, 70%-missing benchmark
resolution the detectable count drops to ~21.5/line — see
[06](06-rtiger-benchmark-dataset.md) and [09](09-rigidity-selection.md).
