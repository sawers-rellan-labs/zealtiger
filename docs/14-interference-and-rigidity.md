# m (interference), the genetic map, and the true rigidity r — analysis notes

**Date:** 2026-06-17
**Context:** Redo of `molbreeding_vs_skim`. Reasoning about where the Stahl
interference parameter `m` comes from, why it is *not* derivable from the map,
and how it relates to the truth-optimal RTIGER rigidity `r`. Companion to
`13-single-locus-validation.md`.

---

## 1. Where m=10 comes from — NOT a fit

m=10 is the **simcross (Broman) default**, adopted because positive interference
is biologically defensible for maize. It is **not estimated** from any data here.
Labelled as such throughout: `run_all.R:24`, `nil_introgression_1400.qmd:19`,
`R/02_simulate.R:140`, `nil_introgression_1400.qmd:388`.

The codebase does a **sensitivity analysis** instead of a fit: it runs both
m=10 (primary) and m=0 (no interference) — `run_all.R:168-169` — and compares
segment-size tails. Result (`docs/05`, `docs/04`): interference preserves the
*mean* CO count (m=10: 33.3 vs m=0: 33.4) and only re-spaces crossovers (m=0 has
a fatter short-segment tail). m=0 is the conservative ceiling.

---

## 2. m is orthogonal to the map, not an addition to it

Two different orders of structure:

- **The genetic map** = first-order / **marginal** crossover intensity. The raw
  observable is the **recombination fraction r̂** (recombinants / total — needs
  genotyping, not a priori). The **map distance d (cM)** is a model-based
  transform of r̂ via a mapping function (Haldane = no interference; Kosambi =
  some), and a Morgan is *defined* as E[# crossovers]. So the map is *observed*
  AND it estimates an *expectation* — estimator/estimand, not a contradiction.
- **m** = second-order **spacing/interference** — the *joint* law of multiple
  crossovers in one meiosis (coincidence). A list of cM positions carries no
  second-order signal: two populations with identical maps can have different m.

**Decisive code fact** (`R/02_simulate.R:151-152`): `sim_from_pedigree(ped, L=L,
m, p)` consumes **only the total cM length per chromosome** — not marker
positions, not internal cM↔bp spacing. Crossovers are a *homogeneous* Stahl
gamma renewal in genetic coordinates; the map's internal shape enters only later
(cM→bp, marker placement, `R/05_make_rtiger_input.R`). So m **cannot** "fit to
the map": the map's shape is not even an input to the crossover generator — only
its length is. Proof m is not additive: m=10 and m=0 give the same CO count
(map sets that); only spacing differs.

m IS estimable from **raw multipoint crossover spacing** (adjacent-interval
coincidence) — e.g. from the **GBTS NILs** (100×, ~50K sites, near-complete
calls). NOT from the published consensus cM map.

---

## 3. Broman / simcross assumptions (as used here)

`sim_from_pedigree(ped, L, m=10, p=0)`, `obligate_chiasma=FALSE` (default):

1. **Stahl / chi-square model** — crossovers are a stationary **gamma renewal
   process** in genetic distance, shape ν = m+1 = 11. p=0 ⇒ pure interference,
   no Poisson escape pathway.
2. **Homogeneous intensity in Morgans** — rate set so E[#CO] = L; uniform per
   unit genetic distance (physical non-uniformity lives in cM↔bp, applied after).
3. **Chromosomes independent** — no inter-chromosome interference.
4. **No sex-specific recombination** (`ignore_sex=TRUE`; one autosomal map).
5. **No obligate chiasma** (default) — CO count not conditioned on ≥1;
   immaterial for >100 cM maize chromosomes.
6. **Two fully-inbred homozygous founders** (recurrent=1, donor=2); no residual
   heterozygosity; genotyping error added only later.
7. **Fixed BC2S3 pedigree, no selection, Mendelian, no segregation distortion.**
8. **Map taken as truth** — fixed known parameter; its estimation uncertainty
   ignored.

---

## 4. m relates to the true rigidity r (via segment size, regime-dependent)

m fixes the true segment-size distribution → a minimum real segment length
`L_floor(m)` in Mb (biology, same on every platform). RTIGER's r sets the minimum
*detectable* segment. The truth-optimal r matches them:

    true_r(platform) ≈ L_floor(m) × informative_marker_density(platform)

Which floor binds depends on the platform — and the binding variable is
**informative density**, set by *depth-driven missingness*, NOT site count:

- **Skim (Source 2):** ~50K sites, 0.4×, ~70% missing → ~**6.6** informative
  markers/Mb (`docs/09`). Detection floor (~0.5–1.5 Mb) is coarse; r tracks
  density, m barely matters. Sweep optimum r=2–3.
- **GBTS / MolBreeding (Source 5):** *target capture*, **same ~50K grid**, 100×,
  ~0% missing → ~**23** informative markers/Mb (~3.5× the skim). NOT WGS, NOT a
  finer grid — just no missingness + confident het.

Because informative density differs ~3.5×, the two platforms should **NOT share
an r**:

| platform | inf. density | r=8 ⇒ physical floor | r for a 1.21 Mb floor |
|---|---|---|---|
| skim | ~6.6 /Mb | 8/6.6 = **1.21 Mb** | 8 |
| GBTS | ~23 /Mb | 8/23 = **0.34 Mb** | ~28 |

**`r=8 for both platforms` (current `molbreeding_vs_skim.qmd:46`) does not give a
matched physical floor.** At r=8 GBTS resolves to ~0.34 Mb — ~3.5× finer than the
skim. If `L_floor(m=10)` ≈ 1 Mb, **GBTS at r=8 resolves below the interference
floor → manufactures spurious short segments.** This is very likely part of the
"GBTS reads more het / apparent over-fragmentation" observed — not all real fine
structure.

**Principled fix:** set r per platform to match a common physical floor:
1. Measure `L_floor(m=10)` from the n=1500 sim (small-segment quantile of the
   true donor-segment-size distribution); m=0 gives the conservative bound.
2. `true_r(platform) = L_floor × informative_density(platform)`.
3. Optionally replace default m with m estimated from GBTS crossover spacing
   (§2), then recompute.

---

## 5. Cross-refs

- Simulation: `R/02_simulate.R`, `run_all.R` (n_nil at `:22`; m=10/m=0 at `:168-169`).
- CO counts / interference sensitivity: `docs/05-crossovers.md`, `docs/04-distribution-and-interference.md`.
- Rigidity sweep & autotune failure: `docs/09-rigidity-selection.md`,
  `docs/09-rigidity-selection.md` (autotune over-rigidifies).
- RTIGER frozen at 649cbf6 (ell_eff reverted): memory ``rtiger-fork``.
- Genome-fraction / ergodicity companion: `13-single-locus-validation.md`.