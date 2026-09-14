# Handoff — CS spatial signature: robustness & the sequencing-depth caveat

**For:** Claude writing the cardiac-sarcoidosis (CS) manuscript.
**You have:** all figures (`.png`) and tables (`.tsv`) in
`results/figures/geneset_robustness/`. This note tells you what each file is,
what it shows, and — critically — **how the spatial results must be framed** so
the manuscript doesn't over-claim.

---

## 1. What the signature is

The CS-cardiomyocyte signature is the four-gene hand-curated set
**`original4` = *GJB7*, *TNNI3K*, *MLIP*, *PANK1*** (published `CM_intrinsic`).
It is scored through scanpy `score_genes(use_raw=False, random_state=0)` on
log-normalized (`normalize_total(1e4)` → `log1p`) expression. The single gene
*TNNI3K* is used both as the headline spatial marker and as a one-gene reduction
of the signature; for one gene the "score" is just its log-normalized value.

Gene-set families referenced in the files:

| set | genes | role |
|---|---|---|
| `original4` | *GJB7, TNNI3K, MLIP, PANK1* | the published 4-gene signature |
| `tnni3k` | *TNNI3K* | single-gene reduction / spatial marker |
| `mlip_pank1` | *MLIP, PANK1* | pairwise reduction |
| `objective53` | 53-gene broad CM program | objective/unbiased comparator |
| `cmstructural` | *TNNT2, MYH7, MYH6, TTN, ACTC1, MYL2, TNNI3, TPM1, ACTN2, MYBPC3* | high-abundance structural-CM positive control |
| `negctrl53` | 53 random genes | negative control (must fail) |
| `inflammation4`, `highconf111` | inflammation/fibrosis programs | directional controls (must behave oppositely) |

---

## 2. TL;DR — what holds, what needs a caveat

**Holds (safe to claim):**
- The CM signature is **depleted toward lesions and in granulomatous/fibrotic
  zones** — an anti-paracrine gradient. `original4` and `mlip_pank1` PASS the
  distance-to-lesion mixed model (β>0, 8/8 patients) and the histologic-zone
  gradient. Directional controls (inflammation/fibrosis programs) move
  **oppositely**, confirming the contrast is real biology, not a scoring
  artifact. This is the core spatial story and it is well-validated across 8–12
  CS sections. Full narrative already written in
  `geneset_robustness_result.md` (the A/B/C/D acceptance matrix).
- The **composite signature score is robust to sequencing depth** — see §3.

**Needs a caveat (do NOT over-claim):**
- The **cross-disease *TNNI3K* comparison** ("*TNNI3K* higher in CS/LMNA than in
  other cardiomyopathies") is **confounded by sequencing depth** and does not
  survive depth equalization. Treat it as depth, not disease specificity. This
  is the load-bearing new finding — details in §3.
- **Every disease has ~1 section** (CS 1, LMNA 2, ARVC 2, HCM 1, Chagas 1), so
  disease is aliased with section/patient. All cross-disease statements are
  **hypothesis-generating**. Keep this caveat on any cross-disease figure.

---

## 3. The sequencing-depth confound (this session's new result)

The cross-disease *TNNI3K* pattern seen at native depth is an artifact of
per-section sequencing depth, established three independent ways that all agree.

**3a. Diagnostic — detection tracks depth, not disease** (`crossdisease_depth_confound.png`,
`depth_diagnostic.tsv`). Across the 7 sections, the fraction of
cardiomyocyte-dominant spots with any *TNNI3K* count rises with median UMI/spot;
the two sections that read "elevated" (CS 10.9k UMI, LMNA-DC 16.2k UMI) are
simply the deep ones. **Where *TNNI3K* is detected, its level is uniformly high
in every disease** (detected-spot median 0.65–1.2 across all 7) — so the
apparent cross-disease difference is a **detection-rate** effect (zeros from
shallow sampling), not a difference in expression level.

**3b. Depth-covariate model** (`crossdisease_depth_model.png`,
`perspot_depth_tnni3k.tsv`, 25,122 CM-dom spots, 7 sections). Logistic model of
*TNNI3K* detection:
- Depth dominates: **OR ≈ 55 per 10× UMI**; depth alone gives pseudo-R² 0.20,
  adding disease raises it only to 0.22 (~90% of explainable variance is depth).
- **After adjusting for depth, CS is not the top disease.** Model-predicted
  detection at common (median) depth: **Chagas 0.72 > LMNA 0.62 ≈ CS 0.61 >
  ARVC 0.47 > HCM 0.41.** CS's own disease coefficient is *below* the Chagas
  reference (OR 0.62). The native-depth ranking that put CS/LMNA on top does not
  hold once depth is controlled.

**3c. Common-depth downsampling** (`crossdisease_commondepth.png`,
`commondepth_grid.tsv`). Every section's CM-dom spots downsampled to a shared
**3,000 UMI/spot** (binomial thinning; spots ≥3,000 retained). At equal depth
the **median *TNNI3K* is 0 in all 7 sections** — CS (0.60→0) and LMNA-DC
(0.85→0) collapse to match the flat ones — and the detection ranking scrambles
(Chagas ties CS; HCM/ARVC lowest). The elevated-vs-flat grid structure is a
depth ordering.

**3d. Aggregation buffers depth — the reassuring part** (`depth_robustness_sweep.png`,
`depth_sweep.tsv`). On a *fixed* cohort of deep spots per section (≥8,000 UMI)
thinned across a depth ladder (same spots throughout), score drift from the
8,000-UMI baseline down to 2,000:

| gene set | drift 8k→2k | behavior |
|---|---|---|
| ***TNNI3K*** alone | **−0.35** (median → 0) | collapses; detection 50%→17% |
| 4-gene CS signature | +0.07 | effectively flat |
| 53-gene objective set | −0.02 | effectively flat |
| 10 structural CM genes | +0.23 | flat/high |

**The depth vulnerability is specific to displaying one low-abundance transcript
as a raw per-section median.** Averaging over genes buffers the dropout, so the
composite classifier the manuscript proposes is far less depth-exposed than the
single-gene spatial panel implied.

---

## 4. File inventory

All paths relative to `results/figures/geneset_robustness/`.

### 4a. Gene-set robustness validation (prior work — the core spatial story)
- `geneset_robustness_result.md` — **authoritative narrative** for analyses A–D
  and the acceptance matrix. Read this first for the anti-paracrine/zone story.
- `summary.tsv` — acceptance matrix, one row per gene set (A/B/C/D verdicts + detail).
- `A_distance_to_lesion_by_geneset.png` / `A_distance_stats_by_geneset*.tsv` —
  distance-to-lesion mixed model (β per SD, 95% CI, patient slope counts).
  `*_ALL7`/`*_ALL8` = variants including additional sets.
- `B_zone_gradient_by_geneset.png` / `B_zone_stats_by_geneset*.tsv` —
  histologic-zone gradient (preserved/granulomatous/fibrotic; spot- and
  patient-level p-values).
- `C_snrna_zone_by_geneset*.tsv` — snRNA pseudobulk zone-flatness (stage
  independence) across 22 CS samples.
- `D_crossdisease_by_geneset.png` / `D_crossdisease_stats_by_geneset*.tsv` —
  cross-disease median ordering per gene set (CS vs controls/ARVC/LMNA/HCM/Chagas).
  `D_crossdisease_perspot.tsv.gz` — per-spot scores backing panel D.
- `acceptance_matrix.png` — the summary.tsv rendered as a figure.
- `panelB_crossdisease.png`, `crossdisease_tnni3k_maps.png`,
  `crossdisease_map_data.tsv` (per-spot: disease, barcode, cm_dom,
  tnni3k_lognorm, px, py) — cross-disease spatial maps.

> **Note for panel D / cross-disease figures:** these are computed at **native
> depth**. Frame per §3 — the CS/LMNA "elevation" for `tnni3k` is depth-driven.
> `original4`/`objective53` medians in panel D are near-zero and do not separate
> CS from mimics within CM-dom spots either (see D_crossdisease_stats).

### 4b. Sequencing-depth analysis (this session)
- `crossdisease_depth_confound.png` + `depth_diagnostic.tsv` — §3a. Columns:
  disease, n_cm, cm_pct, med_UMI, det_pct (detection %), med_all (all-spot
  median), med_det (detected-spot median), mean_det.
- `crossdisease_depth_model.png` + `perspot_depth_tnni3k.tsv` — §3b. Columns:
  disease, section, umi, tnni3k (log-norm value); 25,122 CM-dom spots.
- `crossdisease_commondepth.png` + `commondepth_grid.tsv` — §3c. Columns:
  section, disease, T (`full` or `3000`), n, med_umi, tnni3k_detect,
  med_tnni3k, med_tnni3k_detected, and per-set medians (original4, cmstructural,
  objective53).
- `depth_robustness_sweep.png` + `depth_sweep.tsv` — §3d. Columns: section,
  disease, T (depth 2000–8000), n_fixed, tnni3k_detect, and per-set medians.

### 4c. H&E deep-dives (histology-grounded, single-section illustrations)
- `cs_tnni3k_he_deepdive.png` + `cs_lineage.tsv` — CS_103: *TNNI3K* on H&E,
  preserved-myocardium vs granuloma zoom windows (myocardium 95% CM-dom, mean
  *TNNI3K* 0.83; granuloma myeloid-dominant, no *TNNI3K*). Honest illustration of
  the anti-paracrine contrast **within one section** (not depth-confounded — it's
  a within-section spatial contrast).
- `lmna_tnni3k_he_twosection.png` + `lmna_dc_lineage.tsv` + `lmna_1_lineage.tsv`
  — the two LMNA sections; shows detected-spot medians match (0.92 vs 0.88) while
  detection rate differs (87% vs 49%) — the same depth/detection point at
  single-section resolution. Lineage tsv columns: barcode, dom_lineage, cm_dom,
  tnni3k_lognorm, sc_* lineage scores, px, py.

---

## 5. Manuscript-framing recommendations

**Lead with the within-CS spatial gradient (analyses A/B + H&E deep-dive).** This
is the validated, depth-robust core: the CM signature is depleted toward lesions
and in granulomatous/fibrotic zones, directional controls confirm the contrast.
Use `cs_tnni3k_he_deepdive.png` as the histology-anchored figure.

**For the cross-disease panel, do one of:**
- Reframe as **detection maps with an explicit depth caveat**, and add the
  common-depth panel (`crossdisease_commondepth.png`) or the depth-covariate
  model (`crossdisease_depth_model.png`) as a supplementary control showing the
  native-depth ranking does not survive depth equalization; **or**
- Retire the per-section median cross-disease comparison for *TNNI3K* and make
  the cross-disease point at the **composite-signature** level, where the score
  is depth-robust (§3d).

**Do NOT claim** *TNNI3K* is specifically elevated in CS (or CS+LMNA) versus
other cardiomyopathies. At equal depth that pattern disappears and Chagas
matches CS.

**Do claim** the composite CS signature score is robust to sequencing-depth
variation (`depth_robustness_sweep.png`) — this strengthens the classifier's
credibility and pre-empts a depth-artifact reviewer critique.

**Always caption cross-disease figures** with "single section(s) per disease —
hypothesis-generating; disease aliased with section/patient."

**Style:** gene names italicized (*TNNI3K*, *GJB7*, *MLIP*, *PANK1*).

---

## 6. Methods notes (for reproducibility text)
- Interpretation restricted to **cardiomyocyte-dominant spots** (argmax over
  per-spot lineage scores: Cardiomyocyte / Fibroblast / Endothelial / Myeloid /
  Lymphoid), to keep the spatial analysis faithful to the pooled analysis.
- Scoring: `normalize_total(target_sum=1e4)` → `log1p` → scanpy
  `score_genes(use_raw=False, random_state=0)`; single-gene sets use the
  log-normalized value directly. `score_genes` bins all genes by mean expression
  and samples control genes, so it was always run on the **full** gene matrix.
- Downsampling: `sc.pp.downsample_counts(counts_per_cell=T, random_state=0)` on
  raw counts, restricted to spots with total ≥ T so all retained spots hit
  exactly T (true common depth); binomial thinning.
- Depth model: logistic regression of *TNNI3K* detection on log10(UMI) +
  categorical disease, 25,122 CM-dom spots pooled across 7 sections.
- Sections: CS_103; LMNA-DC & LMNA-2; ARVC_1 & ARVC_2; HCM; Chagas.

---

## ADDENDUM (2026-09-08) — the CS-vs-normal comparison and two more depth checks

Three checks requested after the cross-disease depth work. All three ran on
bigpurple against raw counts; scripts + tables + figures are in this folder.

### 1. CS-vs-normal at common depth — THE UPSTREAM CHECK (passes, with a caveat)

The defining "elevated over normal myocardium" claim compares CS Visium (Foong,
12 sections) to a normal reference (Kuppe, 4 sections) — two different studies.
The raw statistics (GJB7 349/33,628 CS vs ~0/19,775 normal; *TNNI3K* ~3–4×;
*PANK1* ~6–12×) were **cross-study and depth-confounded**, exactly like the
cross-disease numbers.

- **The confound is real:** CS CM-dominant spots are **median 8,786 UMI vs 4,094
  for Kuppe normal — a 2.15× depth asymmetry.**
- **The elevation SURVIVES depth equalization.** Binomial thinning both cohorts
  to a shared depth T (retaining only spots ≥ T, which discards the shallowest
  ~40% of normal spots → conservatively biases normal *upward*):
  - *TNNI3K* detection fold **4.4× → ~2.3×** (stable across T=1k–3k)
  - *MLIP* **3.2× → ~2.8×**
  - *PANK1* **12.4× → ~4.5×**
  - composite score (GJB7/TNNI3K/MLIP/PANK1) is **higher in CS than normal at
    every matched depth** (e.g. T=2000: +0.062 vs −0.188).
- **GJB7 is the exception — down-weight it.** The "349 vs 0" contrast is a
  depth + spot-count artifact; at equal depth GJB7 is essentially undetectable
  in *both* cohorts (0 control spots at some T; fold unstable). Do not lead with
  GJB7 as an "elevated over normal" gene.

**Recommended manuscript language.** Report depth-adjusted (common-depth) folds,
not the raw cross-study folds. State that CS Visium is ~2× deeper than the
normal reference and that after equalizing depth the composite and *TNNI3K /
MLIP / PANK1* remain significantly elevated over normal myocardium (roughly half
the raw log-fold was depth; a real ~2–4.5× elevation remains). Frame GJB7 as a
low-abundance transcript reported for completeness, not a load-bearing marker.

This reconciles cleanly with the cross-disease result: the signature marks
**stressed/cardiomyopathic myocardium vs healthy** (CS > normal at equal depth)
but does **not** distinguish CS from other cardiomyopathies at equal depth
(flat across mimics) — a shared cardiomyopathy program, honestly stated.

Files: `cs_vs_normal_commondepth.png` (fig), `cs_vs_normal_pergene.tsv`
(detection/mean per gene × cohort × depth), `cs_vs_normal_composite.tsv`
(composite score + median UMI), `cs_vs_normal_fold.tsv` (CS/normal folds).
Script: `cs_vs_normal_depth.py`.

### 2. Composite at common depth across diseases — flat; cohort underpowered to order

From `commondepth_grid.tsv`: at the shared downsampling depth (T=3000) the
original four-gene composite median spans a narrow near-zero range across all
diseases. Disease means: HCM +0.053, CS +0.045, Chagas +0.037, ARVC −0.006,
LMNA −0.033. The point estimates place CS near HCM/Chagas rather than LMNA/ARVC,
**but that ordering is not statistically determinable** — only one section for
CS/HCM/Chagas (no replication) and two *discordant* sections for LMNA and ARVC,
whose within-disease spread (LMNA 0.094, ARVC 0.045) ≈ the 0.086 range across
all five disease means. So do **not** state a ranking; state the absence of
depth-resolved separation among cardiomyopathies that this cohort can detect —
the shared-program finding. Distinguishing CS from other cardiomyopathies would
require a replicated multi-section cohort per disease.
Fig: `commondepth_disease_order.png`.

### 3. Per-spot UMI covariate in the distance-to-lesion model (airtight)

Distance-from-lesion mixed model (preserved CM spots, patient random effect),
baseline vs adding within-sample z(log10 UMI) as a covariate. β is per SD of
distance; positive = signature rises in preserved myocardium away from lesions
(anti-paracrine / cell-autonomous).

| gene set | β(dist) baseline | β(dist) **+UMI** | p(+UMI) | β(UMI) |
|---|---|---|---|---|
| composite (original4) | +0.044 | **+0.029** | 2.6e-106 | +0.044 |
| *TNNI3K* | +0.081 | +0.034 | 6.3e-28 | +0.139 |
| *MLIP*+*PANK1* | +0.050 | +0.034 | 2.2e-60 | +0.049 |
| CM-structural (pos ctrl) | +0.162 | +0.138 | ~0 | +0.070 |
| inflammatory (opposite ctrl) | −0.0025 | −0.0023 | 1.4e-7 | −0.0004 |

corr(dist_z, umi_z) = **+0.34** — depth does rise toward preserved myocardium,
so the raw slope is partly depth. But after adjusting for per-spot UMI the
composite's distance slope stays **positive and highly significant**
(attenuated ~34%, not abolished), and the inflammatory control stays negative in
both models. **The anti-paracrine gradient is not a depth gradient.** Reviewer's
depth-gradient objection is answered directly, not only by the directional
controls.

Files: `distance_umi_covariate.png` (fig), `distance_umi_covariate.tsv`,
per-spot `distance_umi_perspot.tsv.gz`. Script: `distance_umi_covariate.py`.

### Method notes (all three)
- CS = Foong `load_geo_visium`; normal = Kuppe h5ads (raw counts in `.raw`,
  Ensembl→symbol via `ensure_symbols`). CM-dominant selection = lineage
  `score_genes` argmax, same as pipeline.
- Common depth = binomial thinning of raw counts to shared T (spots ≥ T only),
  then `normalize_total(1e4)`→`log1p`→`score_genes(use_raw=False,
  random_state=0)`. Detection = raw count > 0.
- Distance model reproduces `run_within_cs`; UMI captured from raw X before
  normalization; `umi_z` = within-sample z of log10(UMI+1).
