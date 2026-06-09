# Diagnosing a heterozygote excess in ancestry calls

**Repo:** [sawers-rellan-labs/zealtiger](https://github.com/sawers-rellan-labs/zealtiger)

An ancestry caller run on the real **SNP50K** VCF returned BC2S3 genotype
frequencies with **~3× more heterozygotes than expected**. This note works
through what can and cannot cause that, because the answer is diagnostic: it
separates *biology* (selection, outcrossing) from *artifacts* (mapping, the
caller) — which is exactly what the ground-truth simulation is built to resolve.

## The observation

`A` = recurrent/B73 allele, `a` = donor/teosinte allele. Expected BC2S3 from the
breeding-matrix recursion (see [02](02-pedigree-and-simulation.md)); observed =
the caller's genome-wide genotype frequencies on SNP50K:

| genotype | expected BC2S3 | observed (caller) | ratio |
|---|---|---|---|
| AA | 0.859 | 0.844 | 0.98 |
| **Aa (het)** | **0.031** | **0.094** | **3.0× ✗** |
| aa (homozygous donor) | 0.109 | 0.063 | 0.58 |
| donor allele freq `aa + ½Aa` | 0.125 | 0.109 | 0.87 |

Two features to explain together: **het is 3× too high**, and **homozygous-donor
is ~½ too low**, while the **donor dosage is only mildly reduced**.

## What it is *not*: directional selection against the donor

Intuitive first guess — teosinte alleles are selected against during NIL
development. But that predicts the **opposite** for het. Selection against the
donor disfavors the `a`-carrying genotypes (`Aa` *and* `aa`), pushing het **down**.
And after three selfings het is intrinsically tiny under *any* purifying scheme —
selfing destroys heterozygosity geometrically (0.25 → 0.125 → 0.0625 → 0.03125),
and selection against the donor cannot create more of it:

| BC2S3 scenario | het (Aa) | aa | donor dosage |
|---|---|---|---|
| neutral | 0.031 | 0.109 | 0.125 |
| complete selection against `aa` (recursion) | 0.035 | ~0 | ↓ |
| **observed** | **0.094** | **0.063** | 0.109 |

Even *complete* removal of homozygous-donor only nudges het to ~3.5% (while
driving `aa` → 0). So directional selection against the donor is excluded as the
cause of the **excess** — and the observed `aa = 0.063` is far too high for strong
selection against `aa` anyway.

## What *can* raise het: balancing selection (overdominant heterosis), and two artifacts

To exceed the ~3% selfing expectation you need a process that *makes* or
*retains* heterozygotes:

- **Overdominant / pseudo-overdominant heterosis = balancing selection.** If the
  heterozygote is genuinely fitter (true overdominance, or repulsion-linked
  deleterious recessives on the two haplotypes — plausible for teosinte
  introgressions), het lineages are preferentially propagated and het is
  *maintained* against the selfing decay. This is the one selection regime that
  elevates het.
  - **Caveat:** only this *overdominant* flavor maintains het. The more common
    **dominance-complementation** basis of maize heterosis is *directional* —
    deleterious recessives are purged as homozygosity rises, homozygotes fix, and
    het is **not** maintained. So "heterosis raises het" holds specifically for
    overdominance/pseudo-overdominance.
- **Paralog / CNV / structural-variant mismapping** (prime artifact suspect for
  teosinte × B73): reads from two divergent copies pile onto one B73 position
  carrying both alleles → a permanent spurious "het."
- **Outcrossing / pollen contamination:** re-introduces genuine het, genome-wide
  and at random.
- **Caller over-assignment:** an HMM/threshold that favors the het state where
  information is thin (we saw RTIGER manufacture het at high rigidity — het
  precision fell to 0.73, [08](08-scoring.md)/[09](09-rigidity-selection.md)).

## The arithmetic favors mis-calling

A simple mis-call model reproduces the whole pattern with no selection term:
move ~0.046 from true-`aa`→called-`het` and ~0.015 from true-`AA`→called-`het`.
That gives the observed `aa↓`, `AA↓`, `het↑`, **and** the dosage drop — each
`aa→het` loses ½ a dosage unit, each `AA→het` gains ½, netting
`−0.046·½ + 0.015·½ ≈ −0.0156`, i.e. 0.125 → 0.109, matching the observed dosage
to the digit. So a calling/mapping artifact is sufficient on its own; balancing
selection is *possible* but not *required* by the numbers.

## How to tell them apart — look at *where* the het sits, and its read signature

Balancing selection is locus-specific; artifacts and outcrossing are not.

| cause | spatial pattern | allele balance at het sites | total read depth |
|---|---|---|---|
| **overdominant heterosis** (balancing) | het tracts at *specific* loci, recombination-bounded, recurrent across lines | ~50/50 (real het) | normal |
| **paralog / CNV mismapping** | scattered at repetitive / SV sites | often *skewed* (e.g. 1:2) | *elevated* |
| **outcrossing / contamination** | uniform, random, genome-wide | ~50/50 | normal |
| **caller over-assignment** | wherever the model is fragile (low-info regions) | inherited from data | normal |

Operationally: real het (heterosis or outcrossing) shows **~50/50 allele balance
and normal depth in contiguous haplotype blocks**; paralog "het" shows **skewed
balance and excess depth at scattered repetitive positions**. Whether the het is
**clustered at consistent loci** (→ heterosis) or **uniform genome-wide** (→
outcrossing) then separates the two real-het causes.

## The simulation pins down the caller's contribution

The benchmark is **neutral** (no selection, no heterosis) and **paralog-free**,
with **truth het = 3%** (the simulation's realized het is 3.0%, matching the
matrix). Run the same ancestry caller on `results/rtiger_benchmark/counts/` and
score it ([08](08-scoring.md), [11](11-benchmarking-other-tools.md)):

- **caller returns ~3% het** → it is faithful, so the real-data excess is
  **biology/artifact** in the *data* — then use the table above (allele balance,
  depth, clustering) to decide heterosis vs paralogs vs outcrossing;
- **caller returns ~9% het on the clean simulation** → the caller **over-assigns
  het** regardless, and the het excess is (at least partly) the *method*.

That cleanly partitions the het excess across **method**, **messy data
(paralogs)**, **outcrossing**, and **balancing selection (overdominant
heterosis)** — which is the entire point of having a ground-truth simulation.

## Summary

- Directional **selection against the donor is excluded** — it lowers het; het is
  ~3% under any purifying scheme after three selfings.
- **Overdominant heterosis is balancing selection** and *can* elevate het — but
  only the overdominant/pseudo-overdominant flavor (not dominance complementation),
  and it predicts **locus-specific** het, not a uniform genome-wide excess.
- The numbers are **fully consistent with het mis-calling** (paralogs / caller),
  so that must be ruled out first — via allele balance, read depth, spatial
  pattern, and the **neutral-simulation control**.
