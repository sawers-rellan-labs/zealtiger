# BC2S3 NIL Introgression-Size Simulation — Specification

**Handover document for a Claude Code session.**
**Author:** Fausto Rodríguez-Zapata
**Confirmed environment:** R 4.5.2, macOS arm64, `simcross` 0.8 (CRAN binary),
`sim_from_pedigree` / `sim_from_pedigree_allchr` / `check_pedigree` /
`convert2geno` all present.

---

## 0. Objective

Forward-simulate **1400 maize BC2S3 near-isogenic lines (NILs)** and estimate
the **distribution of introgression (donor) segment sizes in megabases**, using
a maize genetic map compiled from MaizeGDB consensus-map files and anchored to
the **B73 NAM v5 (`Zm-B73-REFERENCE-NAM-5.0`)** physical assembly.

The job has four stages, run in order:

1. **Compile the map** — parse MaizeGDB consensus files, extract cM + v5 bp,
   clean to a monotonic Marey map, thin to ~2500 markers.
2. **Build the pedigree + simulate** — BC2S3 pedigree, 1400 replicate NILs.
3. **Convert segments cM → Mb** — per-chromosome monotone interpolation.
4. **Summarize + plot** the size distribution.

Stage 1 is the gate: if the cM↔Mb map is wrong, every Mb number downstream is
wrong. Spend QC effort there.

---

## 1. Inputs

- Directory `maizegdb_consensus_map/` containing `consensus_chr01.txt` …
  `consensus_chr10.txt`.
- File structure (verified from the sample):
  - **Line 1** is a free-text title (`MaizeGDB: Details of Map Genetic N`) and
    must be skipped.
  - **Line 2** is the tab-separated column header.
  - Columns are: `Locus`, `Coordinate` (cM), `Bin`, then **repeating 4-column
    blocks** per assembly: `<assembly>_gene_model`, `<assembly>_chr`,
    `<assembly>_start`, `<assembly>_end`, ending in a `Sequence` column.

### Critical parsing gotchas

- **Column sets differ between chromosome files.** chr01 carries B104, PH207,
  W22, Ky21; chr02 carries B97, CML333, Oh7B, Tx303, Tzi8, etc. The v5 block is
  **not at a fixed column index.** Locate it by *name match*, never by position.
- **Header whitespace is dirty** (e.g.
  `Zm-B73-REFERENCE-NAM-5.0_gene_model    ` with trailing spaces). `str_squish()`
  every column name before matching.
- **v5 `chr` values are lowercase** (`chr1`), unlike v3/v4 (`Chr1`). Don't rely
  on case.
- Many rows (telomeres, `K1S2`, genetic-only markers) have **empty physical
  coordinates**. Drop any row lacking a v5 `start`/`end`.
- Rows are **ragged** (trailing empty fields). Read as TSV by name; do not
  assume column count.

### Target v5 columns (match after squish)

```
Zm-B73-REFERENCE-NAM-5.0_gene_model
Zm-B73-REFERENCE-NAM-5.0_chr
Zm-B73-REFERENCE-NAM-5.0_start
Zm-B73-REFERENCE-NAM-5.0_end
```

---

## 2. Environment

Already installed. The session should only verify and add helpers:

```r
# Confirmed present: simcross 0.8
pkgs <- c(
  "tidyverse", "data.table", "fs",
  "furrr", "future",      # parallel replicate simulation
  "scales", "patchwork"   # plotting
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
```

The X11 warning seen at install is unrelated (missing XQuartz); all figures go
to PNG/PDF via ggplot2, so ignore it. Do **not** open base-graphics X11 devices.

---

## 3. R style (strict — follow throughout)

- 2-space indent, 80-char lines, snake_case, tidyverse idioms.
- Exported functions get roxygen2 + `stopifnot()` input validation at the top.
- **Consolidate diagnostics** into one `cat()` block at the end of each chunk;
  never interleave `print()`/`cat()` through computation.
- Set a global seed (`set.seed(20260609)`) and record it in outputs for
  reproducibility.

---

## 4. Stage 1 — Compile the map

**Purpose.** Produce a clean, monotonic per-chromosome Marey map and a thinned
~2500-marker grid.

**Approach.**

1. Read each file with `read_tsv(file, skip = 1)`; `str_squish()` names.
2. Select `Locus`, `Coordinate`, and the four v5 columns by name match.
3. Coerce `Coordinate`, `start`, `end` to numeric; drop rows with `NA` in any.
4. Physical anchor = midpoint: `bp = (start + end) / 2`.
5. Standardize chromosome to integer `1:10` from the v5 `chr` field.
6. **Enforce monotonicity** (the core QC step): within a chromosome, sort by
   `bp`, then keep only markers whose `Coordinate` is non-decreasing along the
   physical axis (running-max filter). This removes mismapped/paralogous
   markers that would invert the map. Report how many are dropped per chr.
7. **Thin to ~2500 markers total**, allocated per chromosome in proportion to
   its cM span, selecting markers as evenly spaced in cM as possible. If a
   chromosome has fewer clean markers than its quota, keep all and report the
   shortfall — 2500 is a *cap/target*, not a guarantee.

**Outputs.**
- `output/maize_map_v5_clean.rds` — all clean monotonic anchors (used for
  cM→Mb interpolation; keep the full set, not just the thinned grid).
- `output/maize_map_v5_grid2500.rds` — thinned marker grid (used as the
  simcross genetic map).
- `output/marey_maps.pdf` — per-chr cM-vs-Mb scatter with the retained
  monotone series overlaid. **Eyeball this before proceeding.**

```r
#' Parse one MaizeGDB consensus-map chromosome file to v5 anchors
#'
#' Reads a single `consensus_chrNN.txt`, locates the B73 NAM v5 columns by
#' name (robust to per-file column differences and dirty whitespace), and
#' returns tidy (locus, cM, chr, bp) anchors with physical coordinates.
#'
#' @param path Path to a consensus_chrNN.txt file.
#'
#' @return A tibble with columns locus, cm, chr (int), bp (numeric).
#' @export
parse_consensus_v5 <- function(path) {
  stopifnot(fs::file_exists(path))

  raw <- readr::read_tsv(path, skip = 1, show_col_types = FALSE)
  names(raw) <- stringr::str_squish(names(raw))

  v5 <- "Zm-B73-REFERENCE-NAM-5.0"
  col_chr   <- paste0(v5, "_chr")
  col_start <- paste0(v5, "_start")
  col_end   <- paste0(v5, "_end")

  stopifnot(all(c("Locus", "Coordinate", col_chr, col_start, col_end)
                %in% names(raw)))

  raw %>%
    dplyr::transmute(
      locus = .data$Locus,
      cm    = suppressWarnings(as.numeric(.data$Coordinate)),
      chr   = suppressWarnings(
        as.integer(stringr::str_remove(.data[[col_chr]], "(?i)^chr"))),
      start = suppressWarnings(as.numeric(.data[[col_start]])),
      end   = suppressWarnings(as.numeric(.data[[col_end]]))
    ) %>%
    dplyr::filter(
      !is.na(.data$cm), !is.na(.data$chr),
      !is.na(.data$start), !is.na(.data$end)
    ) %>%
    dplyr::mutate(bp = (.data$start + .data$end) / 2) %>%
    dplyr::select(locus, cm, chr, bp)
}

#' Enforce a monotone Marey map within each chromosome
#'
#' Sorts markers by physical position and retains only those whose genetic
#' position is non-decreasing (running-max filter), removing map inversions
#' from mismapped or paralogous loci.
#'
#' @param anchors Tibble from `parse_consensus_v5()` (possibly many chrs).
#'
#' @return Tibble of the same columns, monotone in (bp, cm) per chr, plus a
#'   logical `kept` is implied by filtering (dropped rows are removed).
#' @export
enforce_monotone <- function(anchors) {
  stopifnot(is.data.frame(anchors),
            all(c("chr", "cm", "bp") %in% names(anchors)))

  anchors %>%
    dplyr::arrange(.data$chr, .data$bp) %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::mutate(cm_runmax = cummax(.data$cm)) %>%
    dplyr::filter(.data$cm >= .data$cm_runmax) %>%
    dplyr::select(-cm_runmax) %>%
    dplyr::ungroup()
}
```

**Suggested (implement only if requested):** `report_map_qc()` returning
per-chr marker counts before/after monotonicity, cM span, Mb span, and
mean cM/Mb. Consider also an isotonic-regression variant (`stats::isoreg`) as a
sensitivity check against the running-max filter.

---

## 5. Stage 2 — Pedigree + simulation

**Purpose.** Simulate 1400 independent BC2S3 NILs from a 2-founder cross.

**Pedigree.** Founder 1 = recurrent parent (B73, allele 1); founder 2 = donor
(teosinte, allele 2). Recurrent founder is reused for both backcrosses; selfing
is encoded by `mom == dad`.

```
id mom dad sex gen   meaning
 1   0   0   0   0   recurrent parent (B73)
 2   0   0   1   0   donor parent (teosinte)
 3   1   2   0   1   F1
 4   1   3   1   2   BC1  (recurrent x F1)
 5   1   4   0   3   BC2  (recurrent x BC1)
 6   5   5   1   4   BC2S1 (self)
 7   6   6   0   5   BC2S2 (self)
 8   7   7   1   6   BC2S3 (self)  <- the NIL, id 8
```

Validate with `check_pedigree()` before simulating. Sex coding is required by
simcross (0 = female, 1 = male); for autosomal maize it's immaterial but must
be present and consistent.

**Simulation loop.** Replicate the pedigree 1400 times; each call yields one
NIL. The documented API is `sim_from_pedigree(pedigree, L = chr_cm_lengths,
m, p)`, where `L` is a **vector of per-chromosome lengths in cM** (the max cM
per chromosome from the clean map) — *not* a marker map. It returns the
**continuous mosaic**: a list per individual, each with `$mat`/`$pat`, each
haplotype carrying `$alleles` (founder label per interval) + `$locations`
(right-endpoint cM). Extract individual `8` and read those fields directly —
this gives exact cM breakpoints, which is what we want for true segment sizes.
The vignette in hand confirms this structure, so the `donor_segments_cm()`
extractor below is correct as written.

Markers do **not** enter `sim_from_pedigree`. The ~2500-marker grid is used only
(a) as cM↔Mb anchors (Stage 3) and (b) optionally via `convert2geno()` for a
marker-resolution sensitivity check. `sim_from_pedigree_allchr()` is a
genotype-returning wrapper; it is not needed for the continuous-mosaic primary
analysis.

Parallelize replicates with `furrr::future_map()` (`plan(multisession)`); 1400
runs is seconds-to-minutes.

```r
#' Build the BC2S3 NIL pedigree (single replicate)
#'
#' Two inbred founders (1 = recurrent, 2 = donor), F1, two backcrosses to the
#' recurrent founder, then three generations of selfing. Final NIL is id 8.
#'
#' @return A data frame with columns id, mom, dad, sex, gen.
#' @export
bc2s3_pedigree <- function() {
  tibble::tibble(
    id  = 1:8,
    mom = c(0, 0, 1, 1, 1, 5, 6, 7),
    dad = c(0, 0, 2, 3, 4, 5, 6, 7),
    sex = c(0, 1, 0, 1, 0, 1, 0, 1),
    gen = c(0, 0, 1, 2, 3, 4, 5, 6)
  )
}

#' Extract donor (founder-2) segments for one simulated NIL
#'
#' Takes the simcross result for the final individual and returns donor
#' intervals (union across both haplotypes) in cM, per chromosome.
#'
#' @param nil_chr One chromosome's result for the NIL: list with `$mat` and
#'   `$pat`, each having `$alleles` and `$locations`.
#' @param donor_allele Founder label for the donor (default 2).
#'
#' @return Tibble: hap, start_cm, end_cm for donor-carrying segments.
#' @export
donor_segments_cm <- function(nil_chr, donor_allele = 2L) {
  stopifnot(all(c("mat", "pat") %in% names(nil_chr)))

  one_hap <- function(hap, label) {
    locs <- hap$locations
    alle <- hap$alleles
    left <- c(0, head(locs, -1))
    tibble::tibble(
      hap     = label,
      start_cm = left,
      end_cm   = locs,
      allele   = alle
    ) %>%
      dplyr::filter(.data$allele == donor_allele) %>%
      dplyr::select(-allele)
  }

  dplyr::bind_rows(one_hap(nil_chr$mat, "mat"),
                   one_hap(nil_chr$pat, "pat"))
}
```

> **Note on the segment data model.** `donor_segments_cm()` returns segments
> *per haplotype*. An "introgression" as physically detected is the **union of
> donor intervals across both haplotypes** (a region is introgressed if either
> haplotype carries the donor allele). Collapse overlapping mat/pat intervals
> per chromosome (e.g. interval-union via `dplyr` + a sweep, or `IRanges` if you
> prefer) **before** measuring Mb. Decide and document whether you also want the
> zygosity-resolved version (het vs homozygous-donor span), since BC2S3 retains
> some heterozygous donor segments (~3% of loci).

---

## 6. Stage 3 — cM → Mb conversion

**Purpose.** Turn donor cM intervals into physical Mb.

**Approach.** Per chromosome, build a monotone interpolator from the **full
clean anchor set** (Stage 1), mapping cM → bp:

```r
make_cm_to_bp <- function(chr_anchors) {
  stopifnot(all(c("cm", "bp") %in% names(chr_anchors)))
  # collapse cM ties to mean bp so approxfun has a function, not a relation
  d <- chr_anchors %>%
    dplyr::group_by(cm) %>%
    dplyr::summarise(bp = mean(bp), .groups = "drop") %>%
    dplyr::arrange(cm)
  stats::approxfun(d$cm, d$bp, rule = 2)  # rule=2 clamps at ends
}
```

Segment Mb = `(bp_fun(end_cm) - bp_fun(start_cm)) / 1e6`. Linear interpolation
between anchors is the defensible baseline; a monotone spline
(`splinefun(method = "monoH.FC")`) is an optional smoother variant — report
both if results are sensitive near recombination-cold pericentromeres, where
small cM intervals map to large Mb spans (the dominant driver of large physical
introgressions).

---

## 7. Stage 4 — Distribution + outputs

Per-NIL summaries: total donor Mb, donor-segment count, per-segment Mb, and the
**target-locus segment** (Inv4m, chr4) if conditioning is enabled (see §8).

**Figures** (`output/figures/`, PNG + PDF):
- Histogram/density of **total donor Mb per line**.
- Histogram of **individual donor-segment Mb** (pooled).
- ECDF of segment Mb with median + 90th-percentile annotated.
- Optional: per-chromosome facets.

**Tables:** tidy `output/nil_segments.csv` (one row per donor segment, with
nil_id, chr, start_cm, end_cm, start_bp, end_bp, mb) and
`output/nil_summary.csv` (one row per NIL). Also `.rds` of both. Record seed,
marker count, `simcross` version, and interference settings in a
`output/run_metadata.json`.

---

## 8. Modeling decisions to surface (ask Fausto, don't silently pick)

1. **Selection — resolved: none.** Fausto confirmed no selection. Report the
   genome-wide donor-segment distribution from the neutral pedigree as-is. Do
   **not** condition on the target locus; some simulated NILs will carry little
   or zero donor genome, and that is the correct unconditioned result. This is
   exactly the regime where simcross is the right tool and SLiM would be
   overkill.
2. **Crossover interference — pick deliberately.** simcross uses the Stahl
   model: proportion `p` of crossovers from a no-interference process, the rest
   from a chi-square model with interference parameter `m`. **simcross's own
   default is `m = 10, p = 0`** (strong positive interference) — not zero. No
   interference = `m = 0, p = 0` (equivalently `p = 1`). Maize shows positive
   interference, so the package default is biologically defensible and will
   *lengthen* retained donor segments relative to a Poisson model. Recommend
   running the package default, recording `m`/`p` in metadata, and optionally
   comparing against `m = 0` to bound the effect on the tail.
3. **Map total length.** Confirm the compiled map's total cM is in the expected
   maize range (~1400–1600 cM). If the consensus map's max cM per chr looks off
   after cleaning, flag before simulating.
4. **2500-marker grid vs continuous mosaic.** Use the **continuous mosaic** from
   `sim_from_pedigree` for the primary distribution — it gives exact segment
   boundaries. The 2500-marker grid (via `convert2geno`) reflects genotyping
   *detectability* and is a sensitivity analysis only; it is not the simulation
   input.

---

## 9. Validation (must pass)

- **Genome-wide donor fraction** across the 1400 NILs should converge to the
  analytical BC2S3 expectation **0.125** (= f₂ + ½·f₁ = 0.109375 + ½·0.03125).
  This ties back to the matrix-derived frequencies; if the simulated mean donor
  fraction departs from ~0.125 by more than Monte Carlo error, the pedigree or
  allele labeling is wrong — stop and debug.
- Donor-segment **count per genome** should be modest (single digits typical for
  BC2S3) and roughly Poisson in shape.
- No segment should exceed its chromosome's physical length; assert this.

---

## 10. Repo layout + deliverables

```
nil_introgression/
├── R/
│   ├── 01_compile_map.R      # Stage 1 functions + script
│   ├── 02_simulate.R         # Stage 2
│   ├── 03_segments_to_mb.R   # Stage 3
│   └── 04_summarize.R        # Stage 4
├── output/
│   ├── maize_map_v5_clean.rds
│   ├── maize_map_v5_grid2500.rds
│   ├── marey_maps.pdf
│   ├── nil_segments.csv / .rds
│   ├── nil_summary.csv / .rds
│   ├── run_metadata.json
│   └── figures/
├── nil_introgression.qmd     # narrative notebook tying stages together
└── renv.lock                 # pin simcross 0.8 + deps
```

**Final deliverable to Fausto:** the rendered `nil_introgression.qmd` showing
the introgression-size distribution figures, the validation check against 0.125,
and the summary table — plus the saved CSV/RDS for downstream use.

---

## 11. Execution order for the Claude Code session

1. Confirm `?sim_from_pedigree` (0.8) returns the `$mat`/`$pat` +
   `$alleles`/`$locations` structure the extractor assumes; adapt if not.
2. Stage 1; **inspect `marey_maps.pdf`**; report per-chr marker counts and total
   cM. Pause for Fausto if anything looks off.
3. Stage 2 on a small batch (n = 20) first; verify donor fraction ≈ 0.125 in
   miniature before scaling to 1400.
4. Full run (1400, parallel), Stages 3–4.
5. Render the notebook; surface the §8 decisions and §9 validation results.
