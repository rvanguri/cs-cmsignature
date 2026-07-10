# Analysis scripts

These implement the v2.1 scientific logic; the `slurm/*.sbatch` wrappers call them with BigPurple
paths. Real library calls for the tractable parts (knee/manifest, DecontX, scanpy QC, scANVI,
pseudobulk + limma-voom/pydeseq2, gseapy); a handful of `TODO(you)` markers flag the spots that need a
dataset id, credential, or sample→disease map only you have. All Python files pass `py_compile`; R
files are syntax-balanced (run on the cluster where celda/limma are installed).

| Script | Called by | Key inputs | Output | TODO before real run |
|---|---|---|---|---|
| `_common.py` | all py | — | logging, gene panels | — |
| `make_synthetic_data.py` | `smoke_test.sbatch` | `--out` dir | synthetic raw mtx + decontx h5ads + manifest meta | — (dry-run only) |
| `download_geo.py` | 01 | `raw/accessions.tsv` (dataset⇥gse⇥supp_url?) | matrices in `raw/<dataset>/` | confirm each GSE's supp filenames |
| `download_reichart_cxg.py` | 01 | `CXG_DATASET_ID` or `--h5ad_url` | `reichart_cxg.h5ad` + provenance.md | set CELLxGENE dataset id |
| `build_manifest.py` | 01 | `raw/` + optional `raw/sample_meta.tsv` | `manifest.tsv` (+ knee expected_cells) | fill `sample_meta.tsv` (disease/procurement/region) |
| `run_decontx.R` | 03 | manifest + CellBender/CELLxGENE filtered matrices | `<ds>_decontx_counts.h5ad` (carries `sample_id` in obs) | manifest-driven sample selection (done); Reichart reads `raw.X` (done) |
| `run_qc.py` | 04 | `<ds>_decontx_counts.h5ad` + `sample_meta.tsv` | `<ds>_qc.h5ad` (obs has disease/anatomy/sex/individual) | merges sample_meta onto obs by `sample_id`, non-destructive (done) |
| `run_scanvi.py` | 05 | five `*_qc.h5ad` | `integrated_scanvi.h5ad` + model | confirm Reichart label column name |
| `run_gates.py` | 06 | integrated h5ad | gate TSVs + UMAPs | wire raw CellRanger counts for dual-count concordance |
| `run_pseudobulk_de.R` | 07 | integrated h5ad | `DE_<contrast>.tsv` | confirm obs cols: individual/disease/anatomy/sex |
| `run_liu_standalone.py` | 08 | `liu_qc.h5ad` | `DE_Liu_CS_vs_ICM.tsv` | run Liu through 03+04 first (dataset='liu') |
| `run_gsea_figures.py` | 10 | DE + Liu TSVs | volcanoes, GSEA, high-confidence set | — |

## Conventions
- Every script uses `argparse`/`optparse`; the wrappers already pass the right flags.
- `obs` column names assumed: `study`, `disease`, `individual`, `anatomy`, `sex`,
  `predicted_cell_type`. If your h5ads differ, set the `--*_col` flags (DE/gates expose them)
  or rename in `run_scanvi.py` once.
- Non-fatal degradation: missing HeartMap atlas / scib / scrublet warn and continue so the chain still
  produces the open-arm results.
- `manifest.tsv` column order is contractual with `02_cellbender.sbatch` (raw_path = col 7,
  expected_cells = col 8). Don't reorder without updating that script.

## What is genuinely stubbed (needs your input, not just compute)
1. Dataset identifiers / URLs that aren't public-by-convention (Reichart CELLxGENE id; SCP3689 auth).
2. Per-sample disease/procurement/region (`sample_meta.tsv`) — required by the DE model.
3. Dual-count concordance raw-counts wiring in `run_gates.py` (Plan Step 6b).
Everything else runs as-is given the conda env and staged data.
