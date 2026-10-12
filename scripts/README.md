# Analysis scripts

Every script here is part of the short-report pipeline. Run them through `bash scripts/eschf/run_pipeline.sh all`
from the repo root (see `README.md` for the interpreter variables). The driver gives the exact arguments and order.
Rows follow its three rounds; within the last round, `a9_screen` and `a10_levels` run after the per-spot stages.

| Stage (driver label) | Script | Produces |
|---|---|---|
| cluster | `spatial_first/00_export_from_raw.py` | ≥0.5 mm export: per-sample cardiomyocyte pseudobulk beyond 0.5 mm of immune-dominant spots (`data/spatial_first/`) |
| derived_meta | `eschf/derived_inputs.py meta` | section metadata the fits read |
| 01_sens, 02_sens_* | `spatial_first/01_limma_q1q2.R`, `02_leave_one_out.R` | ≥0.5 mm exclusion fit and its leave-one-sample/disease-out refits |
| 01_prim, 02_prim_* | the same scripts on `data/eschf/primary_v2/` | primary fit (all cardiomyocyte-dominant spots) and refits |
| 03_qc | `spatial_first/03_locus_qc.py` | locus QC over the union of both fits' tested genes (`q2_locus_qc.tsv`) |
| 04_*, 05, 06 | `spatial_first/04_snrna_pool.py`, `05_controls.py` (calls `01b_limma_twogroup.R`), `06_published_reference.py` | snRNA pooling and concordance, technical controls, published non-failing comparison |
| 00_inputs, derived_lists | `eschf/00_inputs.py`, `eschf/derived_inputs.py lists` | annotated hit table, gene lists |
| cluster | `eschf/a0_primary_export.py` | primary-fit export (`data/eschf/primary_v2/`) |
| a0_refit, a0_compare | `eschf/a0_primary_refit.R`, `a0_compare.py` | primary refit, universe, comparison with the ≥0.5 mm fit |
| a1 | `eschf/a1_patient_refit.R` | patient pseudobulk and patient-blocked refits |
| a3a4_*, a5_hamed_*, a5b_* | `eschf/a3a4_qc_nonfailing.py`, `a5_hamed.py`, `a5b_procurement.R`, `a5b_pool.py` | QC-consistent statistics, non-failing comparisons, Hamed proteomics, procurement and explant-only snRNA |
| a9_means, a9_screen | `eschf/a9_study_means.R`, `a9_snrna_screen.py` | per-study snRNA means, ±3 screen and source-absence check |
| a11_* | `eschf/a11_raw_donor_check.py`, `a11_snrna_cohort.py` | snRNA CS cohort accounting |
| a10_* | `eschf/a10_patient_log2cpm.R`, `a10_pathway_levels.py` | a priori pathway levels |
| cluster | `eschf/a2a6_spot_export.py` | per-spot program scores and distances (`data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz`) |
| a0_perspot, a2a6_*, a2_model_*, a2_perpatient | `eschf/a0_perspot_scores.py`, `a2a6_analysis.py`, `a2_model.R`, `a2_perpatient.py` | distance bins, mixed model, per-patient slopes, tissue fibrosis |
| a14 | `eschf/a14_pathway_hits_removed.py` | pathway scores with fit hits removed |
| a7 | `eschf/a7_figure.py` | Figure 1, legend, graphical abstract, `a7_panel_values.tsv` |
| a8 | `eschf/a8_summary.py` | `results/eschf/SUMMARY.md` |
| harness | `manuscript_results.py` | `results/manuscript_results.tsv`; non-zero exit on any mismatch |

The `_old` and `_v2` suffixes on driver labels mark the same script run with `ESCHF_PRIMARY=old` (the ≥0.5 mm
exclusion fit, outputs in `results/eschf/`) and `ESCHF_PRIMARY=v2` (the primary fit, outputs in
`results/eschf/primary_v2/`). `eschf/eschf_paths.py` holds the locations for both. `eschf/v3_old_vs_new.py` writes the
`results/eschf/primary_v3/` comparison of the rerun without normal donors against the earlier tree
(`run_pipeline.sh compare`).

## Cluster stages

The three export scripts read the raw Visium objects, which are not redistributed (see `docs/data_availability.md`).
They ran on the cluster, and their outputs are committed. `run_pipeline.sh all` starts from those outputs, and
`run_pipeline.sh cluster1|cluster2|cluster3` prints each job script.

## Options

`03_locus_qc.py --refresh` re-fetches the Ensembl and UCSC annotation instead of using the committed caches.
`05_controls.py --skip-fits` rebuilds its comparison tables without refitting.
