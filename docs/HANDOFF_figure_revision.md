# Handoff — Figure 1 revision for the JACC: Advances Brief Report

**Prepared for:** Claude Science (figure regeneration on the cluster, where the DE tables and Visium
objects live).
**Goal:** revise the single Brief-Report figure to (1) address Reviewer #2's figure comments and
(2) fold in the new cross-disease / per-gene-split story, within the hard limit of **one figure,
maximum two panels**.

Decision already made by the corresponding author: **Choice 1B** — Panel A becomes a grouped
cross-disease dot plot that shows the CS-leaning vs shared split; Panel B stays the spatial
localization panel.

---

## 1. Reviewer #2 figure comments to satisfy (checklist)

- [ ] Panel A: show **statistical significance**, not just log fold-change; color+number both showing
  logFC is redundant → use a **dot plot** (color = logFC, dot size = significance).
- [ ] Center the disease (x-axis) labels.
- [ ] Panel B: add **dashed outlines/shading** to highlight the key regions (granuloma vs preserved).
- [ ] Put a **title directly above each panel**.

Keep it to **two panels** (Brief Report rule). No supplement is permitted.

---

## 2. Panel A — grouped cross-disease dot plot (the redesign)

**Layout.** Genes on the y-axis, **grouped into two labeled tiers**, with a visual separator:
- **CS-leaning:** TNNI3K, GJB7
- **Shared arrhythmogenic-remodeling:** MLIP, PANK1

Diseases on the x-axis, **grouped by modality** with a header band:
- **snRNA (discovery):** DCM, ARVC, HCM, non-failing (CS is the reference).
- **Spatial (Visium):** ARVC, dilated-LMNA, pre-dilated-LMNA, HCM, Chagas (CS is the reference).

**Encodings.**
- **snRNA columns:** dot **color = CS-vs-mimic log2FC** (diverging, colorblind-safe, e.g. blue–white–red,
  symmetric about 0), dot **size = −log10(adjusted P)**. This is the reviewer's requested dot plot.
- **Spatial columns:** n is 1–2 per disease with no valid p-value, so **do not** encode significance as
  size. Represent spatial as **direction glyphs / tiles** (e.g. filled triangle up = CS-higher,
  down = mimic-higher; neutral = comparable), color-matched to the logFC scale but clearly marked as
  directional. A footnote/legend must state spatial is directional (n=1–2), not a significance test.
- Keep the two modalities visually distinct (separate blocks, a divider, and a sub-legend) so no one
  reads spatial glyphs as p-values.

**Annotations (the story):**
- The grouping itself should make the split obvious: CS-leaning tier is CS-elevated across the board;
  shared tier is CS-elevated only over HCM/non-failing (snRNA) and comparable-or-lower in the
  arrhythmogenic mimics.
- Mark the two exceptions explicitly: **TNNI3K** is CS-higher everywhere **except dilated-LMNA**;
  **GJB7** is CS-higher except the **second ARVC section (ARVC_2)** (both near the detection floor).
- Optional small caption cue that HCM appears at two stages (myectomy = discovery/Foong;
  explant = added mimic) and that MLIP/PANK1 rise with that stage.

**Panel A title (above panel):** e.g. *"Cardiomyocyte signature across cardiac sarcoidosis and
cardiomyopathy comparators (snRNA-seq and spatial)."*

---

## 3. Panel B — spatial localization (retain, polish)

Keep the representative CS Visium section: **H&E | immune-cell density (inferred granuloma) |
cardiomyocyte-signature score**.

- [ ] Add **dashed outlines or shading** delineating granuloma vs preserved regions (Reviewer #2).
- [ ] **Panel title above** the panel, e.g. *"Signature localizes to preserved, cardiomyocyte-rich
  myocardium, not granulomas."*
- Keep the **already-fixed CM-intrinsic colormap** (the green/blue/green version was flagged
  hard-to-read in prior review and replaced — do not revert).
- Frame the takeaway as **cardiomyocyte-localized / not paracrine**, consistent with the manuscript's
  softened wording (the localization is largely cardiomyocyte-content; do not label it
  "cell-autonomous" or CS-specific).

---

## 4. Data sources (in-repo) + what must be computed

**Already available:**
- snRNA logFC + adjusted P per mimic: `results/de/de_cm/DE_CS_vs_DCM.tsv`, `DE_CS_vs_ARVC.tsv`,
  `DE_CS_vs_HCM.tsv`, `DE_CS_vs_NF.tsv` (columns `logFC`, `adj.P.Val`). Use the **raw pairwise** values
  for the dot plot, not the confounded model logFC (state this in the legend/source data).
- Spatial per-gene means (original 4 mimics + CS + control):
  `results/figures/cm_spatial_crossdisease/cm_dominant_pergene.tsv`.

**Must be emitted (Claude Science):**
- Per-gene (GJB7, TNNI3K, MLIP, PANK1) **mean expression in CM-dominant spots for all six spatial
  comparators** — the four original mimics plus **HCM (explant)** and **Chagas** — so the spatial
  columns are complete. The geneset-robustness run scored composites and single-gene TNNI3K, but the
  per-gene HCM/Chagas values for GJB7/MLIP/PANK1 still need to be written out (same
  `cm_spatial_crossdisease` scoring, just per gene, for the two added sections).
- A tidy **source-data table** underlying the figure (gene × disease × modality × value × stat).

**Scripts to adapt:** `scripts/brief_report_figure.py` and/or `scripts/composite_figure.py` (the
existing Figure 1 builders); `scripts/cm_spatial_crossdisease.py` for the per-gene spatial emission.

---

## 5. Style / formatting

- Colorblind-safe diverging palette for logFC; symmetric limits; label the scale "CS vs mimic log2FC."
- Center all x-axis labels; horizontal, legible.
- Panel titles above each panel; consistent font with the manuscript.
- Export **both PDF (vector) and PNG (≥300 dpi)** to `figures/`, replacing/augmenting
  `Figure1_JACC_briefreport.{pdf,png}`.
- Two panels only; no third panel, no supplement.

---

## 6. Deliverables

1. Revised `figures/Figure1_JACC_briefreport.pdf` and `.png` (2 panels: A dot plot, B localization).
2. Source-data table (`figures/Figure1_source_data.tsv`).
3. A revised **figure legend** matching the new Panel A (grouped dot plot, snRNA dots vs spatial
   direction glyphs, the two annotated exceptions) and Panel B (shaded regions) — return as text so it
   can be pasted into the manuscript.
4. One-paragraph note confirming which reviewer comments each change addresses.

---

## 7. Caveats to encode in the figure/legend (integrity)

- snRNA dot sizes derive from a **cohort-confounded** discovery contrast (single CS cohort); the dot
  plot uses raw pairwise logFC and the split is the interpretable feature, not any single fold-change.
- Spatial glyphs are **directional at n=1–2 per disease** — explicitly not a powered specificity test.
- The split (CS-leaning TNNI3K/GJB7 vs shared MLIP/PANK1) is the figure's message; it is
  hypothesis-generating, consistent with the manuscript's framing.
