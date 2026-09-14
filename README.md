# Metabolic gene expression differs in non-lesional myocardium between cardiac sarcoidosis and other cardiomyopathies

Analysis code and derived tables for the brief report of the same title.

Spatial transcriptomics (Visium) of cardiac sarcoidosis (CS) and comparator cardiomyopathies
is restricted to **cardiomyocyte-enriched spots at least 0.5 mm from any immune-enriched
spot** — non-lesional myocardium, with granulomas and diffuse immune infiltration excluded.
Differential expression between CS and the comparators in that compartment is then compared
with published single-nucleus data and with published cardiomyocyte comparisons against
non-failing donor myocardium.

## Key result

In cardiomyocyte-enriched non-lesional myocardium, a substrate-handling and mitochondrial
metabolism program — **SLC2A4, ACADVL, IDH2, ACADM, IDH3B, SDHB, HSD17B4**, with additional
ion handling (SLC4A3) and adhesion (ITGA7) genes — is **lower in CS than in comparator
cardiomyopathies**. Of 985 genes tested, 192 differ at FDR 5% (139 lower, 53 higher in CS);
direction is preserved in every leave-one-section-out and leave-one-disease-out refit, and
141 hits remain significant in every leave-one-section-out refit. Direction replicates in
independent single-nucleus data while effect magnitude does not (Spearman rho = 0.03 across
134 genes testable in both), so the claim is about direction, not effect size.

## Reproducing the manuscript numbers

```bash
# R stages need limma + edgeR; Python stages need pandas, scipy, matplotlib, adjustText, openpyxl
RSCRIPT=/path/to/r-env/bin/Rscript PYTHON=/path/to/py-env/bin/python \
  bash scripts/spatial_first/run_all.sh
```

This runs the six analysis stages, rebuilds the figure, and finishes by verifying every
quantitative claim in the manuscript against the regenerated tables. It starts from
**committed inputs** (`data/spatial_first/pseudobulk_T2000.tsv.gz` and friends), because
stage `00_export_from_raw.py` needs the raw spatial objects, which are not redistributed
here — see `docs/data_availability.md`.

`python scripts/manuscript_results.py` can also be run on its own. It prints each claim as
*manuscript value | recomputed value | status*, writes `results/manuscript_results.tsv`,
and **exits non-zero if any reproducible claim fails**. As of this commit: **30 claims
reproduce, 8 mismatch, 2 are not addressable from committed data.** All 8 mismatches are in
the published-non-failing-comparison paragraph and follow from one unresolved question;
`docs/DISCREPANCIES.md` records each one, what the pipeline computes instead, and what was
ruled out. The direction of that finding is unchanged and the recomputed evidence is
stronger than the text claims.

## Pipeline

| Stage | What it does | Main outputs |
|---|---|---|
| `scripts/spatial_first/00_export_from_raw.py` | spot lineage scoring, lesion-distance exclusion, depth downsampling, per-section pseudobulk | `data/spatial_first/*` (committed) |
| `01_limma_q1q2.R` | per-section pseudobulk limma-voom; CS vs comparators (Q1) and diseased vs normal (Q2) | `q1_raw_T2000_d25.tsv`, `q2_raw_T2000_d25.tsv` |
| `02_leave_one_out.R` | refits dropping one comparator section, then one comparator disease | `loo/q1_raw_drop*.tsv` |
| `03_locus_qc.py` | flags segmental duplications, length-model deviation, missing annotation | `q2_locus_qc.tsv` |
| `04_snrna_pool.py` | precision-weighted DCM+ARVC pooling, sign concordance, permutation null | `q1_core_snrna_pergene.tsv`, `q1_core_snrna_concordance.tsv` |
| `05_controls.py` | probe-panel version, panel-matched comparators, cardiomyocyte-content adjustment, mitochondrial-gene removal | `q1_control_comparison.tsv`, `q1_noMT_*.tsv` |
| `06_published_reference.py` | CS-lower genes against published DCM/HCM vs non-failing cardiomyocyte tables | `q1_published_nf_*.tsv` |
| `scripts/figure1_build.py` | the two-panel manuscript figure from committed panel-value tables | `figures/Figure1.{png,pdf}` |
| `scripts/manuscript_results.py` | verifies every manuscript claim; gates the commit | `results/manuscript_results.tsv` |

## Figure

`figures/Figure1.png` / `.pdf`. Panel **a**: spatial vs pooled single-nucleus log2 fold
change for the 134 core genes testable in both platforms, open markers failing locus QC,
shaded quadrants concordant — no line is fitted because the two axes are not calibrated to
each other. Panel **b**: the 18 genes concordant and significant in both arms, with the
published DCM-vs-non-failing and HCM-vs-non-failing effect sizes alongside; open squares are
genes the published table does not report. Both panels use a symmetric-log colour and
magnitude scale so the small single-nucleus effects and the large spatial effects are
readable on one key. Panel values are committed as `results/tables/figure1_panel{A,B}_values.tsv`
and the legend text as `results/tables/figure1_legend.md`.

## Repository layout

- `data/spatial_first/` — committed pseudobulk, detection fractions, section metadata, region
  selection, gene spans, segmental-duplication track, Ensembl annotation
- `data/snrna/` — single-nucleus comparator contrasts and the patient cohort table used by the
  cross-platform arm
- `data/published/` — cardiomyocyte rows distilled from the published supplementary tables
- `results/spatial_first/` — every table the manuscript's numbers come from
- `results/tables/` — figure panel values and legend
- `docs/DISCREPANCIES.md` — where text and code disagree, and why
- `docs/data_availability.md` — accessions and what is not redistributed here

## Earlier arm

An earlier single-nucleus arm of this project nominated a four-gene cardiomyocyte panel
(GJB7, TNNI3K, MLIP, PANK1). That framing is **superseded** by the spatial-first analysis
above and is not what the current manuscript reports; `docs/STATUS.md` records its state and
its limitations, and `results/DE_CS_vs_*.tsv` predate the requantification used here.
