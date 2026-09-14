#!/usr/bin/env bash
# Regenerate every manuscript number and the figure from committed inputs, then verify.
#
# Run from the repo root. Requires the two environments in environment-R.yml (limma, edgeR)
# and environment.yml (pandas, scipy, statsmodels, matplotlib, adjustText).
#
# Stage 00 is NOT run here: it exports the per-section pseudobulk matrix from the raw
# spatial objects, which are not in this repository. Its output is committed as
# data/spatial_first/pseudobulk_T2000.tsv.gz and the pipeline starts from there.
#
# The R and Python stages usually live in different conda environments. Point the
# interpreters at them explicitly if `Rscript` and `python` on PATH are not both right:
#   RSCRIPT=$(conda run -n cs-cmsig-r which Rscript) \
#   PYTHON=$(conda run -n cs-cmsig-py which python) bash scripts/spatial_first/run_all.sh
set -euo pipefail

RSCRIPT=${RSCRIPT:-Rscript}
PYTHON=${PYTHON:-python}

PB=data/spatial_first/pseudobulk_T2000.tsv.gz
DET=data/spatial_first/detection_fraction_T2000.tsv.gz
DETSEC=data/spatial_first/detection_by_section_T2000.tsv.gz
META=data/spatial_first/meta_T2000.tsv
OUT=results/spatial_first

echo "== 01  primary contrasts (CS vs comparators, CS vs normal)"
$RSCRIPT scripts/spatial_first/01_limma_q1q2.R \
  --pb "$PB" --det "$DET" --meta "$META" --thresh 0.25 --tag _T2000_d25 --out "$OUT"

echo "== 02  leave-one-section-out and leave-one-disease-out refits"
for mode in section disease; do
  $RSCRIPT scripts/spatial_first/02_leave_one_out.R \
    --pb "$PB" --det "$DET" --meta "$META" --thresh 0.25 --mode "$mode" --out "$OUT/loo"
done

echo "== 03  locus-level quality control"
$PYTHON scripts/spatial_first/03_locus_qc.py

echo "== 04  cross-platform pooling and sign concordance"
$PYTHON scripts/spatial_first/04_snrna_pool.py

echo "== 05  technical controls (panel, cardiomyocyte content, mitochondrial genes)"
$PYTHON scripts/spatial_first/05_controls.py

echo "== 06  published comparison against non-failing myocardium"
$PYTHON scripts/spatial_first/06_published_reference.py

echo "== figure"
$PYTHON scripts/figure1_build.py

echo "== verify against the manuscript"
$PYTHON scripts/manuscript_results.py
