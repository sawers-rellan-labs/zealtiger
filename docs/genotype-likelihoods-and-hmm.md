# Genotype likelihoods, priors, and their use in an ancestry HMM

This note answers a specific methodological question: *can genotype likelihoods
(GLs), accumulated across a gene / haplotype stretch / fixed bin, rescue
genotype-call confidence at ultra-low coverage?* It is the reviewer question the
nilHMM paper must preempt.

**Short answer.** GLs are the correct currency, and accumulation across a stretch
*does* help — but only for a **shared latent state** (the ancestry of the segment),
which the count-based HMM already does because the read counts *are* the GL's
sufficient statistic. Accumulation does **not** manufacture per-site zygosity
information, and the dilution ceiling (`d ≈ 0.3`,
[snp50k-genotype-identifiability.md](snp50k-genotype-identifiability.md)) is
unchanged by moving from counts to GLs.

Notation: biallelic site, reference (B73) allele $A$, alternate (donor) allele $B$;
genotypes $\{AA, AB, BB\}=\{0/0, 0/1, 1/1\}$; per-base error
$\varepsilon = 10^{-Q/10}$.

---

## 1. The single-site genotype likelihood — identical core in GATK, ANGSD, samtools

All three tools compute the **likelihood** $P(\text{reads}\mid G)$ per genotype with
the same diploid-averaging model (ANGSD calls it "the GATK/dragon model"; samtools
adds a read-dependency term that vanishes at depth 1):

$$
P(\text{reads}\mid G=\{A_1,A_2\}) = \prod_{i}\Big[\tfrac12 P(b_i\mid A_1) + \tfrac12 P(b_i\mid A_2)\Big],
\qquad
P(b\mid A) = \begin{cases} 1-\varepsilon & b = A\\[2pt] \varepsilon/3 & b \neq A. \end{cases}
$$

Each allele of the diploid genotype contributes with weight $\tfrac12$ — this is
why the heterozygote likelihood of any single read is pinned near $\tfrac12$.

### One ALT read ($b=B$, $\varepsilon=0.01$)

| genotype | likelihood $P(\text{read}\mid G)$ | Phred $-10\log_{10}L$ | normalized PL |
|---|---|---:|---:|
| $AA$ | $\varepsilon/3 = 0.0033$ | 24.8 | **25** |
| $AB$ | $\tfrac12\!\cdot\!\tfrac{\varepsilon}{3} + \tfrac12(1-\varepsilon) \approx 0.498$ | 3.0 | **3** |
| $BB$ | $1-\varepsilon = 0.99$ | 0.04 | **0** |

- **PL** $= -10\log_{10} L$, normalized by subtracting the smallest so the best
  genotype is 0 (the GATK worked example: probabilities $100\times$ apart give PLs
  $20$ apart). PL is a *likelihood*, emitted under an effectively flat genotype
  prior — despite GATK's docs writing "$P(\text{Genotype}\mid\text{Data})$".
- **GQ** $=$ (2nd-smallest PL) $-$ (smallest PL), capped at 99. Here **GQ = 3**.

So *by likelihood alone*, one alt read points to **$BB$ (homozygous-alt), GQ 3** —
a $\approx 2{:}1$ call ($\tfrac{1-\varepsilon}{1/2}\approx 2 \Rightarrow 3$ PL).
**That 3-PL likelihood ratio is the entire information content of one read.** No
downstream step adds to it; it can only be combined with a prior or with *other
sites*.

---

## 2. Likelihood → posterior: where the prior lives, and what swapping it does

The genotype **call** and its meaning come from the posterior

$$
P(G \mid \text{reads}) \;\propto\; \underbrace{P(\text{reads}\mid G)}_{\text{GL — identical across tools}} \;\times\; \underbrace{\pi(G)}_{\text{prior — tool-specific, swappable}} .
$$

Where each tool applies $\pi$:

| tool | GL stage | prior stage |
|---|---|---|
| **GATK** | `HaplotypeCaller` PairHMM → PL (likelihood) | `GenotypeGVCFs` allele-freq/HWE prior (`--heterozygosity`, default $10^{-3}$); `CalculateGenotypePosteriors` for **custom** (pedigree / population-AF) priors |
| **ANGSD** | `-GL` models (GATK / samtools / SOAPsnp) → GLs, *no calling by default* | `-doPost 1` = HWE prior at estimated allele freq; `-doPost 2` = uniform |
| **bcftools** | `mpileup` (Li 2011) → GLs | `call -m` allele-freq prior with `-P` mutation-rate hyperprior |

### HWE prior at the cohort frequency → forces HET

With alt-allele frequency $p$ (measured cohort mean $\bar p = 0.015$), HWE gives
$\pi(AA,AB,BB)=((1-p)^2, 2p(1-p), p^2)$. For one alt read:

$$
\frac{P(AB\mid\text{read})}{P(BB\mid\text{read})}
= \frac{\tfrac12\cdot 2p(1-p)}{(1-\varepsilon)\,p^2}
= \frac{1-p}{(1-\varepsilon)\,p} \approx \frac1p \approx 66 .
$$

The prior overturns the $2{:}1$ likelihood and forces the call to **HET** — this is
the 97% REF / 0.00% ALT-hom pattern in the SNP50K matrix
([snp50k-cohort-provenance.md](snp50k-cohort-provenance.md)).

### Breeding prior — the single-locus NIL genotype frequencies

The NIL design is *not* HWE. Selfing after backcross drives heterozygotes to
homozygotes. Single-locus genotype frequencies (donor allele freq $0.125$; selfing
preserves it):

| design | $f_{AA}$ (REF-hom) | $f_{AB}$ (het) | $f_{BB}$ (donor-hom) | het : donor-hom \| introgressed |
|---|---:|---:|---:|---:|
| **BC2S2** (bulked skim) | 0.844 | 0.0625 | 0.0938 | 40 : 60 |
| **BC2S3** (nominal NIL) | 0.859 | 0.0312 | 0.1094 | 22 : 78 |

Note the sign: **$f_{BB} > f_{AB}$** — by design, an introgressed locus is
*more often donor-homozygous than heterozygous* (78% vs 22% for BC2S3), because
selfing removes heterozygosity. (This corrects a claim in an earlier draft that the
introgression is "mostly het" — the *design* expects mostly donor-hom; the callers
*report* mostly het only because donor-hom is undetectable at $0.4\times$, §5.)

### Swapping the prior flips the single-read call

Same one alt read, breeding (BC2S3) prior:

$$
\frac{P(AB\mid\text{read})}{P(BB\mid\text{read})}
= \frac{\tfrac12 f_{AB}}{(1-\varepsilon) f_{BB}}
= \frac{0.5 \times 0.0312}{0.99 \times 0.1094} \approx 0.14
\quad\Rightarrow\quad \textbf{DONOR-HOM} \ (\approx 7{:}1).
$$

So the *identical read and likelihood* is called **HET under the HWE prior** and
**DONOR-HOM under the breeding prior**. This is the crux for the "fix the prior"
argument: swapping HWE → breeding frequencies is correct and worthwhile (it
counteracts the coverage-induced het-bias), **but it only relabels the default —
it adds no zygosity information.** The data still speaks at only $\approx 3$ PL. One
read cannot resolve het from hom under *any* prior; the prior decides the tie.

---

## 3. From one site to a stretch: what accumulation can and cannot do

Accumulating across a gene / haplotype / bin is the right instinct. The subtlety is
**what latent variable the accumulation is evidence *for*.** Two targets:

**(a) The genotype at each individual site.** Borrowing across sites requires an
explicit haplotype/LD model — a reference panel (Beagle/GLIMPSE, Li–Stephens) or a
founder graph (PHG). Multiplying GLs across sites of *different* genotypes is
meaningless without such a model to link them. This route needs the donor
founders / phased panel we lack, and teosinte's short-range LD undercuts it.

**(b) A single latent state shared over the stretch** — the *ancestry*
(REF / HET / donor-hom segment, or "which founder"). If the stretch has one state
$S$, each site's reads are a conditionally-independent emission and the product
concentrates the posterior on $S$:

$$
P(S \mid \text{reads}_{1:n}) \;\propto\; \pi(S)\prod_{j=1}^{n}
\underbrace{\Big[\textstyle\sum_{G} P(\text{reads}_j\mid G)\,P(G\mid S)\Big]}_{\text{emission}_j \;=\; \sum (\text{GL}) \times (\text{state-conditional genotype dist})}.
$$

This is exactly RTIGER / binhmm / PHG's `findPaths`. It works **without a panel**
because the linkage is physical contiguity plus the breeding design (few
recombinations) — not population LD.

### The count emission *is* the GL emission

For a biallelic site under an independent per-read error model, the pair
$(n_{\text{ref}}, n_{\text{alt}})$ is the **sufficient statistic** of the per-site
genotype likelihood: the emission $\sum_G P(\text{reads}_j\mid G)P(G\mid S)$ depends
on the reads only through those counts. nilHMM's count / BetaBinomial emission is a
direct parameterization of this (the BetaBinomial simply adds overdispersion for
mapping/paralog noise). Therefore:

> **"Use genotype likelihoods in the HMM instead of counts" is, for a biallelic
> site, already implemented — phrased as counts.** The caller never hard-calls per
> site (so it never inherits the bcftools 0%-ALT-hom prior artifact); it feeds soft
> per-site evidence into the along-chromosome HMM and decodes the segment.

What *full* GLs add over raw counts is real but marginal at 1 read/site: **per-base
quality weighting** (down-weight low-$Q$ bases rather than counting equally) and
mapping-quality / multiallelic handling. The information axis is the same:
$(n_{\text{ref}}, n_{\text{alt}})$ and their qualities.

---

## 4. Two different HMMs — do not conflate GATK's PairHMM

GATK's HaplotypeCaller *does* use an HMM, but it is orthogonal to the one Ruben
wants:

| | GATK **PairHMM** | **ancestry HMM** (RTIGER / nilHMM / PHG) |
|---|---|---|
| hidden path over | **bases within one read** | **marker sites along the chromosome** |
| emission | read base vs candidate haplotype (indels, base-Q) | reads at a site vs latent ancestry state |
| purpose | local re-assembly / alignment likelihood $P(\text{read}\mid\text{hap})$ | segment ancestry / genotype imputation |

GATK implements the alignment PairHMM, **not** the along-genome ancestry HMM. "They
already do this" is true only for the alignment step; it is not a reusable
component for ancestry segmentation. Keep the two distinct in the paper.

---

## 5. What survives — the information ceiling

Accumulation via the ancestry HMM buys a **confident segment state**, and the
per-marker genotype is then that state **projected/imputed onto each marker** —
PHG-style, and valid for GWAS/QTL (the mapping unit is the segment/dosage
genotype). It is *not* an independent per-site zygosity measurement.

Three limits no GL/HMM refinement removes:

1. **The signal is compressed, not hidden by hard-calling.** Outbred-donor dilution
   puts the states at alt-fractions $\approx\{0,\,d/2,\,d\}=\{0,\,0.15,\,0.30\}$
   with $d\approx0.3$ (empirically confirmed,
   [snp50k-genotype-identifiability.md](snp50k-genotype-identifiability.md) §2.5).
   GLs cannot create information the reads do not carry.
2. **Het vs donor-hom stays weak** — $0.15$ vs $0.30$ separation. A breeding prior
   ($f_{BB}>f_{AB}$) helps counteract the coverage het-bias, but does not resolve
   the zygosity.
3. **The real levers are upstream**, not in the emission: (i) **site selection** to
   donor-informative markers (raise $d$ toward 1, restoring $\{0,\tfrac12,1\}$
   separation), (ii) **fit the emission means per taxon**, (iii) use the
   **BC2S2/BC2S3 prior** in place of HWE.

**Bottom line for the "GLs will rescue it" claim:** GLs are the right soft evidence
and the ancestry HMM is the right accumulator — and both are essentially what the
count-emission caller already does. They rescue the **segment/haplotype** call (and
the imputed marker genotypes), not the isolated site, and they do not lift the
dilution ceiling.

---

## References

- ANGSD genotype-likelihood models — <https://www.popgen.dk/angsd/index.php/Genotype_Likelihoods>
- GATK, *Calculation of PL and GQ by HaplotypeCaller and GenotypeGVCFs* —
  <https://gatk.broadinstitute.org/hc/en-us/articles/360035890451>
- Li H. (2011) *A statistical framework for SNP calling…* Bioinformatics 27:2987 (samtools/bcftools GL model).
- Jensen et al. (2020) sorghum PHG (`data/jensen2020.txt`); Hui et al. (2020) aDNA imputation (`data/hui2020.txt`).

---
*Companion to [snp50k-genotype-identifiability.md](snp50k-genotype-identifiability.md)
and [snp50k-cohort-provenance.md](snp50k-cohort-provenance.md). 2026-07-02.*