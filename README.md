# Metabolic gene expression differs in non-lesional myocardium between cardiac sarcoidosis and other cardiomyopathies

Spatial transcriptomics (Visium) of cardiac sarcoidosis (CS) and comparator cardiomyopathies
is restricted to **cardiomyocyte-enriched spots at least 0.5 mm from any immune-enriched
spot** - non-lesional myocardium, with granulomas and diffuse immune infiltration excluded.
Differential expression between CS and the comparators in that compartment is then compared
with published single-nucleus data and with published cardiomyocyte comparisons against
non-failing donor myocardium.

## Key result

In cardiomyocyte-enriched non-lesional myocardium, a substrate-handling and mitochondrial
metabolism program — **SLC2A4, ACADVL, IDH2, ACADM, IDH3B, SDHB, HSD17B4**, with additional
ion handling (SLC4A3) and adhesion (ITGA7) genes — is **lower in CS than in comparator
cardiomyopathies**. Of 985 genes tested, 192 differ at FDR 5% (139 lower, 53 higher in CS);
direction is preserved in every leave-one-section-out and leave-one-disease-out refit, and
141 of the hits remain significant across all section refits (117 across disease refits). Direction replicates in
independent single-nucleus data.

## Reproducing the manuscript numbers

```bash
mamba env create -f environment.yml      # Python stages, figure, verification
mamba env create -f environment-R.yml    # limma-voom stages
RSCRIPT=$(conda run -n cs-cmsig-r which Rscript) \
PYTHON=$(conda run -n cs-cmsig-py which python) \
  bash scripts/spatial_first/run_all.sh
```

This runs the six analysis stages and rebuilds the figure. It starts from
**committed inputs** (`data/spatial_first/pseudobulk_T2000.tsv.gz` and friends).

## Pipeline

| Stage | What it does | Main outputs |
|---|---|---|
| `scripts/spatial_first/00_export_from_raw.py` | spot lineage scoring, lesion-distance exclusion, depth downsampling, per-section pseudobulk | `data/spatial_first/*` (committed) |
| `01_limma_q1q2.R` | per-section pseudobulk limma-voom; CS vs comparators (Q1) and diseased vs normal (Q2) | `q1_raw_T2000_d25.tsv`, `q2_raw_T2000_d25.tsv` |
| `02_leave_one_out.R` | refits dropping one comparator section, then one comparator disease | `loo/q1_raw_drop*.tsv` |
| `03_locus_qc.py` | flags segmental duplications, length-model deviation, missing annotation | `q2_locus_qc.tsv` |
| `04_snrna_pool.py` | precision-weighted DCM+ARVC pooling, sign concordance, permutation null, gene-length scan | `q1_core_snrna_pergene.tsv`, `q1_core_snrna_concordance.tsv`, `snrna_gene_length_regression.tsv` |
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
- `docs/data_availability.md` — accessions, provenance of the committed inputs, constraints on use

