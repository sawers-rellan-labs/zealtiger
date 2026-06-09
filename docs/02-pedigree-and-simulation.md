# Stage 2 — BC2S3 pedigree + simulation

**Repo:** [sawers-rellan-labs/zealtiger](https://github.com/sawers-rellan-labs/zealtiger)

**Code:** `R/02_simulate.R`, driven by `run_all.R --sim`.
**Output (via Stages 3–4):** `results/nil_segments.*`, `results/nil_summary.*`

Simulate **1,400 independent BC2S3 NILs** from a two-founder cross (founder 1 =
recurrent B73, founder 2 = donor teosinte), **no selection**. The final NIL is
individual id 8 (F1 → BC1 → BC2 → 3× self).

## Three gotchas fixed (vs the handover spec)

1. **`check_pedigree` needs `ignore_sex = TRUE`.** Selfing is encoded as
   `mom == dad`, so a single parent is simultaneously "mom" and "dad"; the
   default male/female consistency check fails (`Some moms are male`). Sex is
   immaterial for autosomal maize.
2. **The pedigree must be a base `data.frame`, not a tibble.** `check_pedigree`
   indexes `pedigree$col` in a way that returns a 1-column tibble for tibble
   input, breaking its id matching (`Some moms not found`).
3. **`sim_from_pedigree(ped, L = <vector>)` returns chromosome-outer, then
   individual-inner**: `sim[[chr]][["8"]]$mat/$pat`, each with `$alleles`
   (founder label per interval) and `$locations` (right-endpoint cM). The spec
   implied "a list per individual"; the extractor indexes `sim[[ci]][["8"]]`.

## Segment model: union vs dosage (important)

`donor_segments_cm` returns donor intervals **per haplotype**. An *introgression*
is the **union of donor intervals across both haplotypes** (a locus is
introgressed if either haplotype carries the donor allele) — `union_intervals`
collapses them per chromosome. This union span feeds the size distribution.

This is **distinct** from the donor *dosage* (allele frequency):

- **dosage** = donor cM summed over *both* haplotypes / (2 × total map cM) → **0.125**
- **union** fraction = P(at least one haplotype donor) ≈ f₂ + f₁ ≈ **0.14**

The skeleton notebook validated the *union* against the *dosage* expectation
(0.125), which is wrong (they differ). We validate **dosage** against 0.125 and
report **sizes** from the **union**. `simulate_nil` returns `$segments` (union),
`$donor_cm_hap` (for dosage), and `$co_chr` (crossovers, see
[05](05-crossovers.md)).

### Where 0.125 comes from (it's an allele frequency, not a genotype frequency)

**0.125 is the donor *allele* frequency (dosage)** — there is no single genotype,
nor any plain sum of two genotypes, equal to 0.125. It is the homozygous-donor
frequency *plus half* the heterozygous frequency. From the BC2S3 genotype
frequencies (`A` = recurrent/B73, `a` = donor/teosinte; breeding-matrix
recursion, [rpubs.com/faustovrz/1312473](https://rpubs.com/faustovrz/1312473)):

| generation | AA | Aa | aa | donor allele freq `aa + ½Aa` |
|---|---|---|---|---|
| BC2 | 0.7500 | 0.2500 | 0.0000 | 0.125 |
| BC2S1 | 0.8125 | 0.1250 | 0.0625 | 0.125 |
| BC2S2 | 0.8438 | 0.0625 | 0.0938 | 0.125 |
| **BC2S3** | **0.8594** | **0.0312** | **0.1094** | **0.125** |

So at BC2S3: `aa + ½·Aa = 0.1094 + ½(0.0312) = 0.125` — identical to the spec §9
form `f₂ + ½f₁ = 0.109375 + ½(0.03125)`.

**Why 0.125, and why it's invariant across the selfing generations:** each
backcross to the recurrent parent *halves* the donor allele frequency, and
selfing *preserves* allele frequency (it only redistributes `Aa → ½AA + ½aa`):

```
F1 q=0.5  →  BC1 q=0.25  →  BC2 q=0.125  (= (1/2)^3)
BC2S1 = BC2S2 = BC2S3 = 0.125            (selfing conserves q)
```

The same genotype frequencies give **both** quantities we use — the difference is
whether the heterozygote is weighted by ½ (dosage) or by 1 (any-donor-present):

| quantity | from genotypes | BC2S3 value | used for |
|---|---|---|---|
| donor **allele freq** (dosage) | `aa + ½·Aa` | **0.125** | §9 validation (simulation must hit this) |
| **donor-present** (union, any donor allele) | `aa + Aa` | **0.1406** | introgression-size distribution / total donor Mb |

This is why the mean **total donor genome ≈ 296 Mb/line** (`0.1406 × ~2,130 Mb`),
while the dosage the simulation validates is 0.125 — two different, both-correct
numbers from the same BC2S3 genotype distribution.

## Validation (must pass — spec §9)

| run | mean donor dosage | expected |
|---|---|---|
| default m=10 | **0.1248** | 0.125 |
| no interference m=0 | 0.1241 | 0.125 |

Pilot (n=20) checks dosage ≈ 0.125 before scaling. Also asserted: no donor
segment exceeds its chromosome's physical length. Replicates run in parallel
(`furrr::future_map`, `plan(multisession)`); 1,400 NILs in ~20 s.

## No selection (spec §8.1)

Confirmed with Fausto: report the genome-wide distribution from the neutral
pedigree as-is. Lines with little/zero donor genome are the correct
unconditioned result; we do **not** condition on a target locus.
