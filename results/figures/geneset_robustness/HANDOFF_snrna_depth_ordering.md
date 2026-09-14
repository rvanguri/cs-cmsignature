# Handoff — snRNA cross-disease ordering and why depth is not the confound it is in Visium

**Question addressed.** (1) Does the cross-disease ordering of the CS composite
signature hold in the single-nucleus RNA (snRNA) cohort? (2) Why is sequencing
depth not the confound there that it is in the spatial (Visium) cohort?

**Bottom line.** In the snRNA cohort the composite orders the diseases cleanly
and **CS is the single most elevated cardiomyopathy** — the opposite of the
Visium, where the cohort was underpowered to order the diseases at all (see
`HANDOFF_three_depth_checks.md`, Check 2). Depth is not a meaningful confound in
snRNA for three structural reasons plus a direct covariate test.

---

## Data and method

- **Cohort.** Per-patient **CM pseudobulk** built by the existing pipeline
  (`scripts/pseudobulk_cm.py`): raw counts summed across all cardiomyocyte
  nuclei per sample. Tables used (already in repo):
  `results/cm_pseudobulk_counts.tsv` (18,789 genes × 142 samples),
  `results/cm_pseudobulk_meta.tsv` (disease, study, n_cells per sample).
- **Diseases / n.** After collapsing CS multi-region samples (LV/RV/SP) to true
  patients: **CS 21, DCM 52, NF (non-failing / normal) 24, HCM 16, ARVC 8.**
  Sources: CS = Neyazi; DCM/ARVC/NF and others = Reichart atlas + chin/larson.
  (Contrast: the Visium had 1–2 sections per disease — hence unorderable.)
- **Scoring.** CPM-log normalize each pseudobulk sample, then composite =
  mean across the four signature genes (*GJB7, TNNI3K, MLIP, PANK1*) of their
  z-scored (across samples) log-CPM. Per-patient score = mean of its samples.
- **Depth test.** OLS `composite ~ disease` vs `composite ~ disease + log10(library size)`,
  reference = NF (normal). Mirrors the Visium per-spot-UMI covariate check.

---

## Result 1 — the ordering holds, and CS is highest

| disease | n patients | composite (mean) | β vs normal | β depth-adjusted | p (adj) |
|---|---|---|---|---|---|
| **CS** | 21 | **+0.73** | **+1.20** | **+1.20** | 1.2e-19 |
| ARVC | 8 | +0.20 | +0.67 | +0.67 | 1.9e-05 |
| DCM | 52 | +0.01 | +0.48 | +0.52 | 1e-07 |
| NF (normal) | 24 | −0.47 | 0 (ref) | 0 (ref) | — |
| HCM | 16 | −1.31 | −0.84 | −0.69 | 2.4e-07 |

CS sits ~1.2 composite units above normal and ~9 SEM clear of the next disease.
This is **not** the Visium point-estimate ordering (which placed HCM "top" and
CS mid-pack) — confirming the Visium ordering was within-noise. With replication
(21 CS patients vs 1 CS section), the snRNA resolves an ordering the Visium
could not, and CS is unambiguously first. HCM sits below normal.

Collapsing CS to true patients (41 samples → 21) does not change the estimates
(CS β +1.20 → +1.20), so within-patient pseudoreplication is not driving it.

## Result 2 — why depth is not the confound

1. **Pseudobulk aggregation removes the dropout mechanism.** The Visium confound
   was a *detection-rate* artifact: shallow spots carry more zeros, so a
   cross-study depth gap deflated detection fractions and manufactured
   fold-changes. Pseudobulk sums thousands of nuclei per patient, so signature
   genes are essentially always detected at the sample level — no zero-inflation
   for depth to act on.
2. **Library-size (CPM) normalization divides depth out by construction.**
   Pseudobulk library size spans **1,275×** across samples (67K → 86M counts);
   CPM-log removes that scaling before scoring (limma-voom applies the same
   normalization in the DE pipeline).
3. **One uniformly-processed integrated atlas, not a cross-study raw comparison.**
   All diseases pass through the same DecontX + scANVI pipeline. There is no
   Foong-CS (8,786 UMI) vs Kuppe-normal (4,094 UMI, 2.15×) depth asymmetry like
   the one that drove the spatial confound.
4. **Empirical covariate test.** Adding per-sample log-depth barely moves the
   disease coefficients (CS +1.20 → +1.20; model R² 0.72 → 0.74); depth's own
   effect is small (β_depth = +0.18). Contrast: in Visium depth explained ~90%
   of the explainable per-spot detection variance.

---

## Caveat to state in the manuscript

CS in the snRNA cohort comes from a **single study (Neyazi)**, so disease and
study are partly confounded. That is a **batch** concern — the scANVI
integration is meant to absorb it, and the pipeline's `cm_disease_distance.py`
flags it with a study-colored panel — and is *separate* from depth. The depth
question specifically is answered cleanly: with depth in the model, CS is the
most elevated cardiomyopathy state and normal/HCM the lowest.

## Files (all in `results/figures/geneset_robustness/`)

- `snrna_disease_order_depth.png` — panel a: composite by disease (per-patient
  points, mean ± SEM, normal reference line); panel b: β vs normal, baseline vs
  depth-adjusted, showing depth barely moves the coefficients.
- `snrna_disease_composite.tsv` — per-disease n, composite mean/SEM/median,
  β vs normal, depth-adjusted β, p, and β_logdepth.
- Inputs: `results/cm_pseudobulk_counts.tsv`, `results/cm_pseudobulk_meta.tsv`.

---

# Robustness of the snRNA ordering — procurement and control-set specificity

The ordering above raises two confound questions. Both were tested; both firm up
the CS-specificity claim rather than undermine it.

## Structural caveat: disease is nearly collinear with study

In this pseudobulk cohort each disease comes from essentially one source:
**CS = Neyazi, DCM = Reichart, ARVC = Reichart, NCC = Reichart; only HCM
(chin2022 + larson) and NF (chin2021 + chin2022 + Reichart) span >1 study.**
So disease and study/procurement cannot be fully separated for CS/DCM/ARVC — any
absolute cross-disease ranking inherits a study/batch caveat (the pipeline's
`cm_disease_distance.py` flags this with a study-colored panel). The two checks
below are the identifiable pieces of that question.

## Check A — HCM procurement split (the two HCM sources do NOT diverge)

HCM has two sources: **larson (GSE174691, myectomy-type) and chin2022
(GSE181764, explant-type)**. If HCM's low position were a tissue-source artifact,
the myectomy and explant HCM should differ. They do not:

| HCM source | n | signature (original4) | cmstructural | log10 depth |
|---|---|---|---|---|
| chin2022 (explant) | 7 | −1.26 | +0.31 | 6.23 |
| larson (myectomy) | 9 | −1.34 | +0.30 | 6.68 |

Both HCM sources give an essentially identical, very low signature despite a
depth difference — HCM's low position is not procurement-driven. Fitting
`composite ~ disease + procurement` (procurement = donor/myectomy/explant)
attenuates the CS signature coefficient from **+1.19 → +0.86** (still p=2e-19):
procurement explains part of the cross-disease spread but nowhere near all of
the CS elevation. (Because explant ≈ CS/DCM/ARVC by construction, this
adjustment is only partial deconfounding — treat +0.86 as a conservative floor.)
File: `snrna_hcm_source_split.tsv`.

## Check B — control sets through the identical pseudobulk ordering (decisive)

Scoring `negctrl53`, `objective53`, and `cmstructural` exactly as the signature,
per patient, same cohort/normalization:

| gene set | CS rank | CS mean | ordering | CS β vs NF | β + procurement |
|---|---|---|---|---|---|
| **original4** (signature) | **1/5** | +0.73 | CS > ARVC > DCM > NF > HCM | +1.19 | +0.86 |
| **objective53** (broader signature) | **1/5** | +1.15 | CS > HCM > NF > ARVC > DCM | +1.64 | +1.19 |
| negctrl53 (negative control) | 2/5 | +0.11 | ARVC > CS > DCM > NF > HCM | +0.22 | +0.14 |
| **cmstructural** (CM-structural / severity) | **5/5** | −0.27 | HCM > DCM > ARVC > NF > CS | −0.22 | −0.30 |

**Interpretation.** CS tops the two real signatures but is **last** on the
generic CM-structural program and only mid-pack (rank 2, small β) on the random
negative control. A generic study/severity/procurement shift would lift
`cmstructural` in CS too — instead CS is the *lowest* disease there. So the CS
elevation is not a global expression shift; it is specific to the signature
genes. This is the single most reassuring result for the specificity claim, and
it holds after procurement adjustment (signature +0.86, cmstructural stays
negative). Figure: `snrna_specificity_procurement.png`;
table: `snrna_specificity_procurement.tsv`.

## Check C — compositional robustness: TMM/RLE + control-bin scoring + permutation null (decisive)

The Check-B scoring (mean of z-scored log-CPM) is the *less* protected method: plain
CPM is not composition-robust, and a z-mean has no abundance-matched subtraction.
The Visium arm used `score_genes` (control-bin subtraction); the pseudobulk arm
did not — so the CS-top arm was the vulnerable one. This check makes both arms
consistent and neutralizes any compositional lift three ways at once:

1. **Composition-robust normalization** — recomputed with edgeR **TMM** and
   **RLE (median-of-ratios)** size factors in place of plain CPM
   (`snrna_compositional_null.tsv`, factors in `normfactors.tsv`; range 0.58–1.38,
   so composition is only mildly skewed to begin with).
2. **`score_genes`-style scoring** — mean(set) minus mean of abundance-matched
   control genes (25 expression bins, 50 controls/bin), the exact routine the
   Visium arm used. This subtraction is what removes a global lift.
3. **Permutation null** — 2,000 random 4-gene and 1,000 random 53-gene
   CM-expressed sets scored identically; the observed sets are placed against
   that null (standardized CS-vs-other-diseases effect *d*).

| gene set | *d* (TMM) | CPM / RLE | CS rank | null percentile | p (right) |
|---|---|---|---|---|---|
| **original4** (signature) | **+2.30** | +2.29 / +2.29 | 1/5 | **100th** | **5e-4** |
| **objective53** | **+2.29** | +2.26 / +2.31 | 1/5 | **100th** | **<1e-4** |
| negctrl53 (neg ctrl) | +1.21 | +1.32 / +1.26 | 1/5 | **95th** | 0.05 |
| cmstructural (CM struct.) | −0.96 | −0.96 / −0.96 | 5/5 | 8.5th | 0.92 |

**Verdict — the compositional explanation is finished.**
- CS stays **rank 1 on the signature under TMM and RLE** and control-bin scoring;
  the effect barely moves across normalizations (+2.29 → +2.30 → +2.29), so
  library-size composition is not driving it.
- The signature sits in the **extreme right tail of the permutation null**
  (100th percentile, p≈5e-4) — far beyond what random CM gene sets produce.
- **negctrl53 lands exactly at the null's 95th percentile (p=0.05).** Its
  nominal "rank 1" is a boundary chance draw of one random 53-gene set, not a
  specific signal — its effect (+1.21) is about half the signature's and equals
  the null's own 95th percentile. Read the *magnitude against the null*, not the
  rank: the negative control behaves like a random set, exactly as it should.
- **cmstructural is in the left tail (8.5th percentile, CS-low).** A global
  compositional/severity lift would push CS *up* here; instead CS is the lowest
  disease on the sarcomere program under every method. The global-shift scenario
  is excluded, and CS-low sarcomere expression is itself interpretable
  (dedifferentiation-like).

Figure: `snrna_compositional_null.png` (a: permutation null with the four sets
marked; b: CS effect stable across CPM/TMM/RLE, null 5–95th band shaded).
Table: `snrna_compositional_null.tsv`.

**Specificity claim: established.** The CS elevation is signature-specific,
composition-robust, and beyond the permutation null; the only residual caveat is
the single-study origin of CS (a batch, not a compositional, concern — best
addressed by the within-study Liu CS-vs-ICM contrast).
