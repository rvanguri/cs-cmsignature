# Analysis scripts

## Manuscript pipeline (spatial-first arm)

These produce every number in the current brief report. Run them with
`bash scripts/spatial_first/run_all.sh` from the repo root; see `README.md` for the
interpreter overrides the two-environment split needs.

| Script | Inputs | Output |
|---|---|---|
| `spatial_first/00_export_from_raw.py` | raw Visium objects (not redistributed) | `data/spatial_first/*` |
| `spatial_first/01_limma_q1q2.R` | pseudobulk, detection fractions, section metadata | `q1_raw_T2000_d25.tsv`, `q2_raw_T2000_d25.tsv` |
| `spatial_first/01b_limma_twogroup.R` | same, plus a two-group metadata file | any single two-group contrast (used by 05) |
| `spatial_first/02_leave_one_out.R` | same, `--mode section\|disease` | `loo/q1_raw_drop*.tsv` |
| `spatial_first/03_locus_qc.py` | Q2 results, gene spans, segdup track, annotation | `q2_locus_qc.tsv` |
| `spatial_first/04_snrna_pool.py` | Q1 results, refits, `data/snrna/snrna_CS_vs_*.tsv` | `q1_stable_down_115.tsv`, `q1_core_snrna_pergene.tsv`, `q1_core_snrna_concordance.tsv` |
| `spatial_first/05_controls.py` | pseudobulk + the control metadata files | `q1_control_comparison.tsv`, `q1_noMT_*.tsv`, `ctrl_panelB_vs_panelA.tsv` |
| `spatial_first/06_published_reference.py` | `data/published/*`, Q1 results | `q1_published_nf_pergene.tsv`, `q1_published_nf_summary.tsv` |
| `figure1_build.py` | `results/tables/figure1_panel{A,B}_values.tsv` | `figures/Figure1.{png,pdf}` |
| `manuscript_results.py` | everything above | `results/manuscript_results.tsv`, non-zero exit on any mismatch |

Conventions for this arm: R stages take `--pb/--det/--meta/--thresh/--out` and write
tab-separated tables with a `gene` column; Python stages default their paths to the
committed locations and need no arguments. `05_controls.py --skip-fits` rebuilds the
comparison tables without refitting, and `03_locus_qc.py --refresh` re-fetches the Ensembl
annotation instead of using the committed copy.

## Earlier single-nucleus arm

The scripts below implement the v2.1 scientific logic of the earlier arm (see
`docs/STATUS.md`); the `slurm/*.sbatch` wrappers call them with cluster
paths. Real library calls for the tractable parts (knee/manifest, DecontX, scanpy QC, scANVI,
pseudobulk + limma-voom/pydeseq2, gseapy); a handful of `TODO(you)` markers flag the spots that need a
dataset id, credential, or sample→disease map only you have. All Python files pass `py_compile`; R
files are syntax-balanced (run on the cluster where celda/limma are installed).

| Script | Called by | Key inputs | Output | TODO before real run |
|---|---|---|---|---|
| `_common.py` | all py |: | logging, gene panels |: |
| `make_synthetic_data.py` | `smoke_test.sbatch` | `--out` dir | synthetic raw mtx + decontx h5ads + manifest meta |: (dry-run only) |
| `download_geo.py` | 01 | `raw/accessions.tsv` (dataset⇥gse⇥supp_url?) | matrices in `raw/<dataset>/` | confirm each GSE's supp filenames |
| `download_reichart_cxg.py` | 01 | `CXG_DATASET_ID` or `--h5ad_url` | `reichart_cxg.h5ad` + provenance.md | set CELLxGENE dataset id |
| `build_manifest.py` | 01 | `raw/` + optional `raw/sample_meta.tsv` | `manifest.tsv` (+ knee expected_cells) | fill `sample_meta.tsv` (disease/procurement/region) |
| `run_decontx.R` | 03 | manifest + CellBender/CELLxGENE filtered matrices | `<ds>_decontx_counts.h5ad` (carries `sample_id` in obs) | manifest-driven sample selection (done); Reichart reads `raw.X` (done) |
| `run_qc.py` | 04 | `<ds>_decontx_counts.h5ad` + `sample_meta.tsv` | `<ds>_qc.h5ad` (obs has disease/anatomy/sex/individual) | merges sample_meta onto obs by `sample_id`, non-destructive (done) |
| `run_scanvi.py` | 05 | five `*_qc.h5ad` | `integrated_scanvi.h5ad` + model | confirm Reichart label column name |
| `run_gates.py` | 06 | integrated h5ad | gate TSVs + UMAPs | wire raw CellRanger counts for dual-count concordance |
| `run_pseudobulk_de.R` | 07 | integrated h5ad | `DE_<contrast>.tsv` | confirm obs cols: individual/disease/anatomy/sex |
| `run_liu_standalone.py` | 08 | `liu_qc.h5ad` | `DE_Liu_CS_vs_ICM.tsv` | run Liu through 03+04 first (dataset='liu') |
| `run_gsea_figures.py` | 09 | DE + Liu TSVs | volcanoes, GSEA, high-confidence set |: |

## Conventions
- Every script uses `argparse`/`optparse`; the wrappers already pass the right flags.
- `obs` column names assumed: `study`, `disease`, `individual`, `anatomy`, `sex`,
  `predicted_cell_type`. If your h5ads differ, set the `--*_col` flags (DE/gates expose them)
  or rename in `run_scanvi.py` once.
- Non-fatal degradation: missing scib / scrublet warn and continue so the chain still
  produces the open-arm results.
- `manifest.tsv` column order is contractual with `02_cellbender.sbatch` (raw_path = col 7,
  expected_cells = col 8). Don't reorder without updating that script.

## What is genuinely stubbed (needs your input, not just compute)
1. Dataset identifiers / URLs that aren't public-by-convention (Reichart CELLxGENE id).
2. Per-sample disease/procurement/region (`sample_meta.tsv`): required by the DE model.
3. Dual-count concordance raw-counts wiring in `run_gates.py` (Plan Step 6b).
Everything else runs as-is given the conda env and staged data.
