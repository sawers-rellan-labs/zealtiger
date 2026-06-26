# Local-ancestry admixture mapping for BZea NILs — ELAI + the "Two teosintes" approach, and a proposed methodology

**Date:** 2026-06-26
**Context:** The objective is to use **local ancestry / chromosome painting** (teosinte
introgression in a B73 background) to **associate with traits** across the BZea NILs.
This note records the assessment of ELAI (Guan 2014) as a reference-based ancestry caller,
the admixture-mapping methods from Yang et al. 2023 ("Two teosintes made modern maize",
Science 382, eadg8940; `data/yang2023.pdf`), the conceptual clarifications established in
discussion, a verification of feasibility against the data we actually have, and a concrete
**proposed methodology**. Companion to `16-parameter-search-and-rigidity-mle.md` and the
caller-comparison work. Memory: `bzea-skim-bulked-bc2s2`, `b73-control-impostors`,
`wideseq-jaccard-donor-id`, `nilhmm-counts-path`.

---

## 1. Objective and the key reframe

Goal: map teosinte-introgression effects on traits. The crucial reframe established in
discussion:

> **Admixture mapping uses LOCAL ANCESTRY as the predictor, not SNP genotypes.**

Two distinct paradigms:
- **reliable genotypes → ordinary GWAS** (SNP allele → trait);
- **local ancestry → admixture mapping** (donor-ancestry state → trait). ← *this is what
  we have and want.*

So the worry "I only have a binary ancestry call, not reliable genotypes" is misplaced: the
chromosome paintings are the **native input** to admixture mapping, not a degraded SNP
substitute. The only thing binary costs vs. a dosage is the HET/HOM additive resolution —
which on 0.4× skim is unidentifiable anyway (depth-1 degeneracy), and which the **binary
recurrent-vs-donor collapse deliberately discards**.

### 1.1 Conceptual clarifications (for the record)
- **"Locus" = a position on the ancestry track, not each SNP.** Local ancestry is
  piecewise-constant (changes only at breakpoints); all SNPs inside one block carry the
  identical dosage. Effective resolution and the multiple-testing burden are set by the
  number of distinguishable ancestry **segments**, not the SNP count.
- **The HMM smooths noisy genotypes into runs** — this is the whole value. The transition
  prior (ELAI's switch rates, governed by γ = generations) penalises frequent switching, so
  an isolated mis-genotyped SNP cannot flip the call; ancestry at any position borrows
  strength from the whole surrounding block. Ancestry inference therefore **supersedes noisy
  per-SNP genotypes** — which is why it works at 0.4× where per-SNP GWAS would not.
  *Caveats:* (i) the run-length prior is only right if **γ is set to our generation**, not a
  default; (ii) smoothing fixes noise, not **non-identifiability** — the depth-1 HET/HOM
  degeneracy is an emission-level limit no HMM can undo.
- **Bins are a legitimate test unit; you do not fabricate SNPs.** The ancestry dosage is a
  step function defined along the whole chromosome. Testing at a bin (or bin midpoint) is
  *sampling that function*; the value there is real because it is the bin's aggregated
  read signal, not a phantom genotype. Mapping resolution then = bin/segment size (~Mb),
  which is the intrinsic ceiling of admixture mapping anyway.
- **Prefer posterior dosage over a hard 0/1.** Every caller (RTIGER, nilHMM, Skim-BIN, ELAI)
  yields a posterior P(donor); use the expected dosage as the covariate to propagate calling
  uncertainty. Collapse to binary only as the robust fallback for skim.

---

## 2. ELAI (Guan 2014) — the reference-based caller

ELAI is a two-layer HMM (upper = source population/ancestry; lower = haplotype cluster).
Full model below (Genetics 196:625). It differs from our existing callers (RTIGER / Skim-BIN
/ nilHMM) by using **reference panels for both sources** and an explicit **haplotype-cluster**
layer and **generation (γ) prior** in the transition — genuinely more information than the
per-marker "B73 vs not-B73" count callers, hence plausibly more precise for the 2-way
maize/teosinte call.

### 2.1 The two-layer HMM (LaTeX)

```latex
% ELAI two-layer HMM (Guan 2014, Genetics 196:625) — amsmath, bm
\textbf{States.} For individual $i$ at marker $m$:
\[ X_m^{(i)}\in\{1,\dots,S\}\ \text{(upper: source population)},\quad
   Y_m^{(i)}\in\{1,\dots,K\}\ \text{(lower: haplotype cluster)},\quad Z_m=(X_m,Y_m). \]

\textbf{Emission.} With cluster allele freq $\theta_{mk}$ and error term $\mu$,
$t_{mk}=\theta_{mk}(1-\mu)+(1-\theta_{mk})\mu$. Diploid genotype $g\in\{0,1,2\}$ from
clusters $j,k$:
\[ p(g\mid j,k)=\begin{cases} t_{mj}t_{mk}, & g=2\\
   t_{mj}(1-t_{mk})+(1-t_{mj})t_{mk}, & g=1\\ (1-t_{mj})(1-t_{mk}), & g=0.\end{cases}\]

\textbf{Initial.} $p(X_1=s,Y_1=k)=\alpha_s^{(i)}\,\beta_{1sk}$, with admixture proportions
$\alpha_s^{(i)}$ and lower-given-upper probs $\beta_{msk}$.

\textbf{Transition} ($c_m$ = inter-marker genetic distance):
\[ j_m=\exp(-t_m^{(u)})\ \text{(upper no-switch)},\quad t_m^{(u)}=4N_ec_m\delta_u;\qquad
   r_m=1-\exp(-t_m^{(l)})\ \text{(lower switch)},\quad t_m^{(l)}=4N_ec_m\delta_l,\ \delta_l\approx1/K. \]
\[ p(s,k\mid s',k')=\underbrace{(1-j_m)\alpha_s^{(i)}\beta_{msk}}_{\text{upper switches}}
   +\underbrace{j_m r_m\beta_{msk}\mathbb{1}(s{=}s')}_{\text{lower switches}}
   +\underbrace{j_m(1-r_m)\mathbb{1}(s{=}s')\mathbb{1}(k{=}k')}_{\text{neither}}. \]
(``if the upper layer switches, the lower must switch'' is enforced by the first term.)

\textbf{Admixture timing.} $\sum_m t_m^{(u)}=\gamma\sum_m c_m$, so mean tract length (cM)
$\lambda=M/\sum_m t_m$ obeys $\gamma\lambda=100$ (more generations $\Rightarrow$ shorter tracts).

\textbf{Output.} Local-ancestry dosage $D_m^{(i)}(s)=\sum_{c\in\{1,2\}}\Pr(X_m^{(i,c)}=s\mid g^{(i)},\xi^\ast)\in[0,2]$
(Yang's mexicana dosage; ``ELAI score $>1.8$'' $\equiv$ dosage $>0.9$/chromosome).
```

> **Note on the convention:** the paper defines $j_m=\exp(-t^{(u)})$ (a *no-switch* prob)
> but $r_m=1-\exp(-t^{(l)})$ (a *switch* prob); the transition above places both consistently.

---

## 3. The Yang 2023 admixture-mapping methods

Three analyses, in increasing scope:

1. **Per-locus admixture mapping.** Regress phenotype on **ancestry dosage** at each marker,
   with kinship to control structure; filter loci below a minimum donor frequency (Yang: ≥5%
   in inbreds, ≥1% in landraces). This is the core — start here.
2. **Multivariate / G×E (JointGWAS).** A generalised-least-squares F-test per locus testing
   the ancestry effect jointly across response dimensions, with a pre-estimated covariance:
   - *Inbred panel:* joint across **33 traits** (covariance via **MegaLMM**: genetic +
     non-genetic; kinship from IBS SNPs). One test per locus for "ancestry affects any trait."
   - *Across environments:* joint across **13 trials** (covariance via MegaLMM, factors =
     #trials−1) → a **G×E admixture scan**.
   Collapses to (1) if you have a single trait in one environment.
3. **Variance partitioning.** Two relationship matrices — a background SNP GRM and an
   **ancestry GRM** built from the dosage matrix (OSCA `--efile … --make-bod`) — into REML
   (LDAK) to split heritability into *introgression component* vs *background* vs error.
   Validated by a phenotype-simulation grid (h² × introgression-fraction); code at
   `rossibarra/maize_origins/simulate_phenotypes.Rmd`.

---

## 4. Feasibility on the SNP50K data — verified

| Requirement | Status (verified 2026-06-26) |
|---|---|
| Teosinte reference panel | ✅ `data/donor_id/refpanel_gt.tsv` — **49,002 SNPs × 217 accessions**, dense full GT |
| NIL 50K data | ✅ `data/rtiger_50K/counts/<donor>/<name>.tsv` — 82 donors, ~52K sites, `chr pos ref refcount alt altcount` |
| **Shared SNP set** | ✅ both start at `chr1 29860` — same 50K sites; intersect to common ~49K |
| Maize/B73 source | ✅ already staged — see §5 |
| Genetic map (for position file) | ✅ consensus map from the cM work (`docs/01`–`03`) |

Why a reference-based call is plausibly **more precise** here, and robust:
- Only a **2-way maize/teosinte** call is wanted → **pool all 217 accessions into one
  teosinte source**. Species-level sparsity (which sank the per-accession Jaccard donor-ID,
  `wideseq-jaccard-donor-id`) is irrelevant; we get maximum reference depth.
- The **binary back-conversion sidesteps the depth-1 HET/HOM degeneracy**: collapsing
  "teosinte present vs B73" merges HET and teosinte-hom anyway, so the degeneracy that makes
  *dosage* unreliable on skim does **not** affect the *binary* target.
- ELAI is built for **sparse/missing** genotypes (imputes through the HMM via the reference
  haplotypes) — the 73%-missing skim regime is what it was designed for.

---

## 5. The maize source is B73 — not a panel

B73 is the reference genome and the 50K sites are ascertained as teosinte-vs-B73
polymorphisms, so **B73's alt-allele frequency ≈ 0 at every site by construction**. There is
nothing to estimate — exactly how nilHMM / the wideseq painting code already fix the REF
emission mean at 0. All discriminating signal lives in the teosinte source; the maize source
is a degenerate fixed all-reference haplotype (the ideal 2-way setup).

ELAI is supervised, so B73 must still be *provided* as source-2, but that is trivial and
already staged:
- **Idealised B73** — the `is_reference=TRUE` `B73` row (group SS) in the panel; matches how
  the other callers treat REF; zero contamination risk. **Default choice.**
- **Clean real B73 field rep** — to capture the platform noise floor; use a verified-clean
  rep only. 12 skim B73 field reps exist (incl. `PN10_SID893`), and `b73_checks_qc` found
  **no impostors among skim B73s** (worst 0.054; impostors were all BRB-seq,
  `b73-control-impostors`).

---

## 6. PROPOSED METHODOLOGY (given our data)

### Stage A — produce a binary/dosage ancestry matrix via ELAI
1. **Install ELAI** (`yongtaoguan/elai`, compiled binary).
2. **Source 1 (teosinte):** all 217 accessions from `refpanel_gt.tsv`, pooled.
   **Source 2 (maize):** the idealised `B73` reference row (default), optionally a clean field
   rep for noise modelling.
3. **NIL input:** intersect to the common ~49K sites; convert read counts → ELAI/BIMBAM
   genotypes per site: `alt≥1 & ref=0 → 2`, `ref≥1 & alt=0 → 0`, `both≥1 → 1`, `0 reads →
   missing`. (Or BIMBAM mean-genotype dosages; on skim this stays ~0/2, so it adds little.)
4. **Position file** per chromosome with genetic distances from the consensus map (feeds $c_m$).
5. **Run per chromosome:** `-C 2 -c ~5–10 -s ≥20`, and — critically — **`-mg` = the NIL
   backcross generation (small, ~2–4), NOT Yang's 6,000.** Wrong `mg` → wrong tract-length
   prior → over/under-segmentation. This is the single most important parameter for our
   long-tract NIL regime (skim = bulked **BC2S2**, `bzea-skim-bulked-bc2s2`).
6. **Output** `ps21.txt` = teosinte dosage 0–2 per NIL per SNP → **threshold to binary**
   (e.g. dosage >0.9 → teosinte; else B73) → **sample on the bin grid** (the painting test
   unit). This is the NIL × locus ancestry matrix.

### Stage B — admixture mapping
1. **Ancestry GRM** from the dosage matrix (OSCA) for structure control; optionally a
   background SNP GRM and/or donor-taxon covariate (donor identity correlates with which
   segments a NIL carries — avoid confounding "which donor" with "which locus").
2. **Per-locus mixed model:** `trait ~ ancestry_dosage(locus) + random(ancestry GRM)`
   (GEMMA/EMMAX/LDAK). **Frequency-filter** loci (per-locus donor frequency is often low;
   power concentrates where multiple NILs share donor at a locus).
3. **If multiple traits / environments:** the JointGWAS GLS F-test (MegaLMM covariance) for a
   multivariate or G×E scan; otherwise this reduces to step 2.
4. **Variance partition:** SNP GRM + ancestry GRM into REML (LDAK) → fraction of trait
   heritability attributable to teosinte introgression genome-wide; validate via the
   simulation grid.

### Stage C — robustness (leverage the caller-comparison harness)
Run the mapping on each caller's ancestry track (ELAI, RTIGER, Skim-BIN, nilHMM) and keep
peaks that **replicate across callers**. A peak in only one caller = likely a calling
artifact. This turns the existing multi-caller infrastructure into the mapping's validation
layer.

---

## 7. Caveats / open items
- **`mg` calibration:** set to the NIL generation; if uncertain, sweep a small range and
  check tract-length sanity against the BC2S2 fragment-size expectation (`docs/16`,
  `bc2s3-fragment-size-expectation`).
- **Mapping resolution** = bin/tract size (~Mb), intrinsic to admixture mapping, not a defect.
- **Depth-1 HET/HOM** stays unidentifiable on skim → the binary collapse is the honest call;
  dosage-resolved mapping needs higher depth (BRB-seq ~2.8×).
- **Whether ELAI actually beats the count callers** for the 2-way call is **empirical** —
  benchmark head-to-head (Stage C), do not assume.
- **Setup costs:** install ELAI; intersect/convert counts→BIMBAM; build the position/map file.
