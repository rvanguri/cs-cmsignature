# Cardiomyocyte transcription in cardiac sarcoidosis

Analysis code, derived inputs and result tables for the ESC Heart Failure short report of the same
title. Every number in the manuscript text and figure legend is regenerated from the committed inputs and
checked by `scripts/manuscript_results.py`.

## Data and design

- **Spatial transcriptomics (Visium, probe-based).** 12 cardiac sarcoidosis (CS) samples from 9 patients
  (GSE314910) and 6 comparator samples from 6 patients: 2 arrhythmogenic (ARVC), 2 lamin A/C-associated dilated
  (LMNA DCM), 1 hypertrophic (HCM) and 1 Chagas cardiomyopathy. Each spot is thinned to 2,000 UMIs and assigned to
  its dominant lineage (argmax over five lineage gene-set scores). Samples with at least 50
  cardiomyocyte-dominant spots are analysed (11 CS samples, 8 patients).
- **Differential expression.** limma-voom on per-sample pseudobulk of cardiomyocyte-dominant spots, CS vs
  comparators, 985 genes passing the 25% detection rule. Refits: patient pseudobulk, patient blocking,
  exclusion of spots within 0.5 mm of an immune-dominant spot, and leave-one-sample-out and leave-one-disease-out.
  Locus-level quantification QC on every tested gene.
- **Distance to inflammation.** Per-patient slopes of the score of the genes lower in CS, and of single genes,
  against distance to the nearest immune-dominant spot.
- **Single-nucleus RNA-seq.** CS (Neyazi) against pooled DCM and ARVC (Reichart) and against non-failing
  myocardium, with a pre-stated ±3 log2 screen on the pooled snRNA effect, random-gene-set expectations, and an
  explant-only subset.
- **Published reference.** Chaffin et al. 2022 cardiomyocyte DCM vs non-failing and HCM vs non-failing tables.

Accessions, the provenance of each committed input, and constraints on use are in `docs/data_availability.md`.
Two sets of publisher supplements (Neyazi Tables S1 and S2; Hamed supplementary files) are not redistributed. Download them into
`data/external/` as listed there before running the pipeline.

## Reproducing the manuscript numbers

```bash
mamba env create -f environment.yml      # Python stages, figure, verification
mamba env create -f environment-R.yml    # limma-voom and mixed-model stages
PYTHON=$(conda run -n cs-cmsig-py which python) \
RSCRIPT=$(conda run -n cs-cmsig-r which Rscript) \
  bash scripts/eschf/run_pipeline.sh all
```

This runs every analysis stage from the committed exports, rebuilds Figure 1 and the graphical abstract, writes
`results/eschf/SUMMARY.md`, and finishes by checking every manuscript number against the regenerated tables. It takes
about two minutes on a laptop and, once the supplements are downloaded, needs no network access. Logs go to `$LOGDIR` (default `logs/`).

The three export stages that read the raw Visium objects (`00_export_from_raw.py`, `a0_primary_export.py`,
`a2a6_spot_export.py`) ran on the cluster next to the raw data, and their outputs are committed under `data/`. The
raw objects are not redistributed. `run_pipeline.sh cluster1|cluster2|cluster3` documents those runs.

## Verification

`python scripts/manuscript_results.py` prints each claim as *manuscript value | recomputed value | status*, writes
`results/manuscript_results.tsv`, and **exits non-zero on any mismatch**. As of this commit **865 claims reproduce,
none mismatch, and 2 are NOT_RUN**: the patient-blocked refits of the ≥0.5 mm exclusion fit, where the
within-patient correlation is not estimable. The `ShortReport` block (62 claims) encodes each number as printed in
the manuscript and compares it at printed precision. The `ESCHF[...]` blocks (805 claims, including the 2 NOT_RUN) recompute the underlying
tables independently. `docs/DISCREPANCIES.md` describes both blocks, the two NOT_RUN claims and the revision history.

## Layout

| Path | Contents |
|---|---|
| `scripts/eschf/` | short-report stages, the driver `run_pipeline.sh`, `derived_inputs.py` |
| `scripts/spatial_first/` | the limma-fit stages both fits share (limma fits, leave-one-out, locus QC, snRNA pooling, controls, published reference) |
| `scripts/manuscript_results.py` | the verification harness |
| `data/eschf/` | primary-fit exports (`primary_v2/`), per-spot scores, H&E JPEGs for Figure 1A (`img/`) |
| `data/spatial_first/` | ≥0.5 mm-fit pseudobulk, detection, section metadata, gene spans, segmental duplications, Ensembl annotation |
| `data/snrna/`, `data/published/`, `data/external/` | snRNA contrasts and patient cohort; Chaffin tables; GEO sample listing (publisher supplements downloaded here) |
| `results/eschf/` | every short-report table, `SUMMARY.md`, `a7_panel_values.tsv` (every plotted value) |
| `results/spatial_first/` | ≥0.5 mm-fit tables, locus QC and its annotation caches |
| `figures/` | `eschf_figure1.{png,pdf}`, `eschf_figure1_legend.md`, `eschf_graphical_abstract.png` |
| `docs/` | `DISCREPANCIES.md`, `data_availability.md` |

Gene-label placement in Figure 1C is not deterministic, so a rerun can move a few labels. Plotted values do not
change.

## Earlier versions

- Tag `brief-report-spatial-first`: the spatial-first brief report (*Metabolic gene expression differs in
  non-lesional myocardium between cardiac sarcoidosis and other cardiomyopathies*), with its driver, figure and
  claim block.
- Tag `pre-slim-2026-09-14`: the superseded single-nucleus arm (four-gene panel GJB7, TNNI3K, MLIP, PANK1).

```bash
git show brief-report-spatial-first --stat
git checkout brief-report-spatial-first -- <path>
```
