# GATK does not call genotypes on posteriors

A note settling a recurring question: **does GATK apply a genotype prior — i.e.
does it call genotypes on the posterior, the way bcftools does?** No. GATK assigns
each sample's genotype as the maximum-*likelihood* genotype; the population prior
only decides whether the *site* is variant, and never rewrites the per-sample call.

Companion to
[genotype_likelihoods_and_hmm.html](https://sawers-rellan-labs.github.io/zealtiger/genotype_likelihoods_and_hmm.html)
(§ "HWE prior at the cohort frequency → forces HET"), which covers the
posterior/HWE callers (bcftools, freebayes, ANGSD `-doPost`, Beagle/GLIMPSE) that
*do* call on the posterior.

## What GATK actually does

- **PL = genotype likelihoods** (normalized, Phred-scaled), not posteriors. Despite
  the GATK docs writing "P(Genotype | Data)", the emitted PL is a likelihood under
  an effectively flat genotype prior.
- **GT = argmax likelihood** (the smallest PL). From the docs, verbatim: *"The
  largest likelihood, which corresponds to the most likely genotype, is picked and
  assigned to the sample."*
- **GQ** = (2nd-smallest PL) − (smallest PL), capped at 99 — a likelihood-ratio
  confidence, not a posterior probability.
- The **heterozygosity / allele-frequency prior** (`--heterozygosity`, default
  0.001; infinite-sites Pr{AC = i} = Het / i) enters **only** the site-level variant
  probability: QUAL = −10 · log Pr{AC = 0 | D}. The docs are explicit: *"the
  heterozygosity parameter has nothing to do with the likelihood of any given sample
  having a heterozygous genotype … [it] doesn't actually change the output genotype
  likelihoods at all."*
- So the prior gates *"is this a variant site?"* — it does not touch the per-sample
  GT or the PLs. The only way to get prior-informed GATK genotypes is the opt-in
  `CalculateGenotypePosteriors` (adds a `PP` field and can re-assign GT using a
  population-AF or pedigree prior).

## Consequence for one ALT read

For a single ALT read the genotype likelihoods are ε/3 : ½ : (1 − ε) for
{hom-ref, het, hom-alt}. GATK picks the largest → **hom-ALT (1/1)**. This is the
*opposite* of a posterior caller like bcftools, whose rare-allele HWE prior converts
the same read into a **HET** call (see the HWE table in the companion doc). The
divergence is entirely the prior: same reads, same likelihoods, different call.

## The "GATK overcalls hets" impression — real, but a different mechanism

GATK het-excess is genuine and widely reported, but it is **data/likelihood-driven,
not prior-driven**, and it appears at **moderate-to-high coverage**, not ultra-low:

- Because GATK is max-likelihood, it calls HET only when the pileup **actually
  contains both alleles**. At ~1× (one read/site) it rarely does, so it hom-calls —
  it does *not* over-call HET at ultra-low coverage; if anything it over-calls
  hom-ALT on isolated alt reads.
- Where the het-excess comes from: **paralogs / segmental duplications / CNV**
  (reads from a second locus mismap in, injecting minority alt reads alongside ref
  reads → the likelihood legitimately favors HET) and **contamination** (foreign
  minority alleles). The classic collapsed-repeat "het-hell". These need a mixed
  pileup, i.e. depth.

## Contrast

| caller | genotype call | HET-inflation mechanism | regime | one ALT read → |
|---|---|---|---|---|
| bcftools `call -m` | MAP (**posterior**) | rare-allele **prior** pushes single alt reads to HET | ultra-low coverage | **HET** |
| GATK HaplotypeCaller | argmax **likelihood** (min PL) | max-likelihood picks up both alleles when the (often artifactual) pileup contains them | moderate+ coverage, paralogs/contamination | **hom-ALT** |

Both can inflate HET, for non-overlapping reasons. For ultra-low-coverage skim
(~0.4×) the bcftools/posterior mechanism is the one that bites; GATK at 1× would
instead scatter spurious hom-ALT calls on isolated alt reads. Neither is trustworthy
per-site at that depth — which is why the ancestry HMM (soft read counts pooled over
a segment) is the right unit, not the per-site genotype.

## Sources

- Assigning per-sample genotypes (HaplotypeCaller) —
  <https://gatk.broadinstitute.org/hc/en-us/articles/360035890511>
- Calculation of PL and GQ by HaplotypeCaller and GenotypeGVCFs —
  <https://gatk.broadinstitute.org/hc/en-us/articles/360035890451>
- GenotypeGVCFs (heterozygosity / infinite-sites prior) —
  <https://gatk.broadinstitute.org/hc/en-us/articles/360037594731>
