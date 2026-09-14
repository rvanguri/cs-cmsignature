# Handoff — three sequencing-depth checks on the CS spatial signature

Date: 2026-09-08. All analyses ran on bigpurple against **raw counts**;
cardiomyocyte-dominant spots only; scoring reproduces the pipeline exactly
(`normalize_total(1e4)`→`log1p`→`score_genes(use_raw=False, random_state=0)`;
single-gene sets use the log-normalized value directly). Common depth = binomial
thinning of raw counts to a shared T, retaining only spots with total ≥ T so
every retained spot sits at exactly T. Scripts, tables, and figures are in
`results/figures/geneset_robustness/`.

Composite signature = **GJB7, TNNI3K, MLIP, PANK1** (CS-cardiomyocyte set).
Cohorts: CS = Foong Visium (12 sections); normal = Kuppe control (4 sections).

---

## Check 1 — CS vs normal at common depth (the upstream check)

**Question.** The defining "elevated over normal myocardium" claim compares CS
Visium to a normal reference from a *different study*. GJB7 349/33,628 CS vs
0/19,775 normal is a pure detection-rate statistic across two studies, and the
*TNNI3K* ~3–4× and *PANK1* ~6–12× folds are cross-study too. If the normal
reference is shallower, the claim inherits the same depth confound as the
cross-disease work. Run the same common-depth downsampling.

**Answer — the confound is real, but the elevation survives it.**

- The confound exists: CS cardiomyocyte-dominant spots are **median 8,786 UMI
  vs 4,094 for Kuppe normal — a 2.15× depth asymmetry.** The raw cross-study
  folds were inflated.
- At matched depth (thinning both cohorts to shared T; retaining spots ≥ T
  discards the shallowest ~40% of normal spots, which biases *normal upward* and
  makes the surviving fold **conservative**), the elevation persists and is
  stable across T = 1,000–3,000:

  | gene | detection fold CS/normal, native | at common depth | verdict |
  |---|---|---|---|
  | *TNNI3K* | 4.4× | **~2.3×** | survives |
  | *MLIP* | 3.2× | **~2.8×** | survives |
  | *PANK1* | 12.4× | **~4.5×** | survives (attenuated) |
  | GJB7 | 204× | unstable (≈0 in both) | **artifact — down-weight** |

- The **composite score is higher in CS than normal at every matched depth**
  (e.g. T=2000: CS +0.062 vs normal −0.188; the CS>normal ordering holds at
  T=1k/1.5k/2k/3k and native).
- **GJB7 is the exception.** Its "349 vs 0" contrast is a depth + spot-count
  artifact; at equal depth it is essentially undetectable in *both* cohorts
  (0 control spots at some T). It is not a load-bearing "elevated over normal"
  gene and should be reported for completeness only.

**Bottom line for the manuscript.** Report depth-adjusted (common-depth) folds,
not the raw cross-study folds. State that CS Visium is ~2× deeper than the
normal reference and that after equalizing depth the composite and *TNNI3K /
MLIP / PANK1* remain elevated over normal myocardium — roughly half the raw
log-fold was depth, but a real ~2–4.5× elevation remains. This is consistent
with Check 2: the signature marks **stressed/cardiomyopathic myocardium vs
healthy**, not **CS vs other cardiomyopathy**.

Files: `cs_vs_normal_commondepth.png`, `cs_vs_normal_pergene.tsv`,
`cs_vs_normal_composite.tsv`, `cs_vs_normal_fold.tsv`; script
`cs_vs_normal_depth.py`.

---

## Check 2 — composite across diseases at common depth

**Question.** commondepth_grid.tsv already has the original four-gene medians at
T=3,000. If the composite is flat across diseases at equal depth, that is the
shared-program finding stated cleanly.

**Answer — flat, and the cohort is underpowered to order the diseases.** At the
shared downsampling depth the four-gene composite median spans a narrow,
near-zero range across all cardiomyopathies. Disease means at T=3000:

| disease | composite @ T=3000 | n sections |
|---|---|---|
| HCM | +0.053 | 1 |
| **CS** | **+0.045** | 1 |
| Chagas | +0.037 | 1 |
| ARVC | −0.006 | 2 (−0.029, +0.017) |
| LMNA | −0.033 | 2 (−0.080, +0.014) |

The point estimates place CS among HCM/Chagas rather than LMNA/ARVC — but **this
ordering is not statistically determinable and should not be interpreted.** The
design has only one section for CS, HCM, and Chagas (no replication → within-
disease variance unestimable) and two *discordant* sections for LMNA and ARVC.
Where within-disease spread can be measured it already ≈ the between-disease
spread: **LMNA's two sections span 0.094 and ARVC's span 0.045, vs a 0.086 range
across all five disease means.** So the between-disease ranking is within noise;
the Visium cohort is not powered to resolve any ordering of the mimics.

The correct statement is therefore the *absence* of separation, not a ranking:
combined with Check 1, the signature is elevated in diseased/stressed myocardium
relative to normal (Check 1) but shows **no depth-resolved separation among
cardiomyopathies** that this cohort can detect (Check 2). A shared cardiomyopathy
program; distinguishing CS from other cardiomyopathies would require a replicated
multi-section cohort per disease.

Files: `commondepth_grid.tsv` (prior session); `commondepth_disease_order.png`
(disease ordering with the within- vs between-disease spread annotated).

---

## Check 3 — per-spot UMI covariate in the distance-to-lesion model

**Question.** Fibrotic/granulomatous tissue yields less RNA, so a reviewer will
ask whether the distance-from-lesion gradient is a depth gradient. The opposing
directional controls mostly answer it; adding per-spot UMI as a covariate makes
it airtight.

**Answer — the gradient is not a depth gradient.** Depth does rise toward
preserved myocardium (corr(distance_z, UMI_z) = **+0.34**), so the objection is
legitimate — but after adding within-sample z(log10 UMI) to the mixed model the
anti-paracrine gradient holds (β per SD of distance; + = signature rises in
preserved myocardium away from lesions; patient random effect):

| gene set | β(dist) baseline | β(dist) **+UMI** | p(+UMI) | β(UMI) |
|---|---|---|---|---|
| composite (original4) | +0.044 | **+0.029** | 2.6e-106 | +0.044 |
| *TNNI3K* | +0.081 | +0.034 | 6.3e-28 | +0.139 |
| *MLIP* + *PANK1* | +0.050 | +0.034 | 2.2e-60 | +0.049 |
| CM-structural (pos. control) | +0.162 | +0.138 | ~0 | +0.070 |
| inflammatory (opposite control) | −0.0025 | −0.0023 | 1.4e-7 | −0.0004 |

The composite distance slope stays **positive and highly significant** after
adjustment (attenuated ~34%, not abolished), and the inflammatory control stays
**negative** in both models. The reviewer's depth-gradient objection is answered
directly, not only via the directional controls.

Files: `distance_umi_covariate.png`, `distance_umi_covariate.tsv`, per-spot
`distance_umi_perspot.tsv.gz`; script `distance_umi_covariate.py`.

---

## One-paragraph summary

CS Visium is ~2× deeper than the normal reference and deeper than most of the
cardiomyopathy mimics, so all raw cross-study folds were depth-confounded. After
equalizing depth: (1) the composite and *TNNI3K/MLIP/PANK1* remain elevated in
CS over **normal** myocardium (GJB7 does not — it is a detection artifact);
(2) the composite is **flat across cardiomyopathies**, so it is not CS-specific;
(3) the spatial anti-paracrine gradient survives per-spot-UMI adjustment and is
therefore not a depth gradient. Net: a robust cardiomyopathy/stress program with
a genuine spatial gradient, honestly bounded — not a CS-specific diagnostic at
the level of the raw fold-changes.
