#!/usr/bin/env bash
# Reproduce the ESC HF short-report analysis (Foong CS + six comparator Visium sections; no normal donors)
# from the raw matrices to results/manuscript_results.tsv.
#
# Three exports need the raw Visium matrices and run on the cluster (BigPurple, SLURM; header below).
# Each round needs outputs from the round before, so the order is:
#   cluster1 -> local1 -> cluster2 -> local2 -> cluster3 -> local3
#
#   bash scripts/eschf/run_pipeline.sh cluster1            # print the round-1 job script (submit it on the cluster)
#   bash scripts/eschf/run_pipeline.sh local1  <round1_out> # install the discovery export, run the sensitivity arm
#   bash scripts/eschf/run_pipeline.sh cluster2            # print the round-2 job script + list of files to stage
#   bash scripts/eschf/run_pipeline.sh local2  <round2_out> # install the primary all-CM export, run the primary arm and A-scripts
#   bash scripts/eschf/run_pipeline.sh cluster3            # print the round-3 job script + list of files to stage
#   bash scripts/eschf/run_pipeline.sh local3  <round3_out> # install the per-spot export, run the per-spot A-scripts, figure, harness
#   bash scripts/eschf/run_pipeline.sh all                 # rerun every local stage from the committed exports
#                                                          # (localN without <roundN_out> also skips the install)
#   bash scripts/eschf/run_pipeline.sh compare <old_tree> <outdir>   # optional: claim-by-claim diff against another tree
#
# Interpreters: PYTHON (default: python; environment.yml) and RSCRIPT (default: Rscript; environment-R.yml).
# 05_controls.py launches R itself and reads RSCRIPT from the environment. Logs go to LOGDIR (default: logs/).
set -euo pipefail
PY=${PYTHON:-python}; RS=${RSCRIPT:-Rscript}; export RSCRIPT=$RS
SHIM=$(mktemp -d); ln -sf "$(command -v "$RS")" "$SHIM/Rscript"; export PATH="$SHIM:$PATH"
LOG=${LOGDIR:-logs}; mkdir -p "$LOG"
run(){ n=$1; shift; echo "== $n"; "$@" > "$LOG/$n.log" 2>&1 || { echo "$n FAILED (see $LOG/$n.log)"; exit 1; }; }

SLURM='#SBATCH --time=03:00:00
#SBATCH --partition=cpu_short
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
PY=/gpfs/scratch/vangur01/heartmap-cs/conda/heartmap/bin/python
COMP=/path/to/comparator_visium   # folder holding the six comparator Visium sections'

D=data/spatial_first; S=results/spatial_first; DP=data/eschf/primary_v2; V=results/eschf/primary_v2

# Publisher supplements are not redistributed; download them into data/external/ first (docs/data_availability.md).
need_external(){ miss=""
  for f in neyazi_table_s1.xlsx neyazi_table_s2.xlsx 44161_2026_880_MOESM4_ESM.xlsx 44161_2026_880_MOESM5_ESM.xlsx \
           44161_2026_880_MOESM7_ESM.xlsx; do [ -f "data/external/$f" ] || miss="$miss $f"; done
  [ -z "$miss" ] && return 0
  echo "missing in data/external/:$miss"
  echo "  neyazi_table_s1/s2.xlsx: Supplemental Tables S1 and S2, saved under these names, of Neyazi et al., Circulation 2026, doi:10.1161/CIRCULATIONAHA.126.079304"
  echo "  44161_2026_880_MOESM{4,5,7}_ESM.xlsx: supplementary files (publisher download names) of Hamed et al., Nat Cardiovasc Res 2026, doi:10.1038/s44161-026-00880-w"
  exit 2; }

case "${1:-}" in
cluster1)   # stage: scripts/spatial_first/00_export_from_raw.py, data/spatial_first/spatial_first_config.json
  printf '%s\nln -sfn $COMP comparators\n$PY 00_export_from_raw.py 2>&1 | tee run_log.txt\n# outputs: out/*\n' "$SLURM" ;;

local1)     # discovery export -> sensitivity (>=0.5 mm) arm
  if [ -n "${2:-}" ]; then O=$2
  for f in pseudobulk_T2000.tsv.gz detection_by_section_T2000.tsv.gz detection_fraction_T2000.tsv.gz region_selection.tsv; do cp "$O/$f" $D/$f; done; fi
  run derived_meta $PY scripts/eschf/derived_inputs.py meta
  run 01_sens $RS scripts/spatial_first/01_limma_q1q2.R --pb $D/pseudobulk_T2000.tsv.gz --det $D/detection_fraction_T2000.tsv.gz --meta $D/meta_T2000.tsv --thresh 0.25 --tag _T2000_d25 --out $S
  # metadata/tested_genes_985.txt keeps its historical file name; its content is the sensitivity arm's tested set
  $PY -c "import pandas as pd; q=pd.read_csv('$S/q1_raw_T2000_d25.tsv',sep='\t'); open('metadata/tested_genes_985.txt','w').write('\n'.join(q.gene)+'\n')"
  for mode in section disease; do
    run 02_sens_$mode $RS scripts/spatial_first/02_leave_one_out.R --pb $D/pseudobulk_T2000.tsv.gz --det $D/detection_fraction_T2000.tsv.gz --meta $D/meta_T2000.tsv --thresh 0.25 --mode $mode --out $S/loo
  done ;;

cluster2)   # primary all-cardiomyocyte arm export
  echo "# stage: scripts/eschf/a0_primary_export.py data/spatial_first/{spatial_first_config.json,pseudobulk_T2000.tsv.gz,region_selection.tsv} metadata/{tested_genes_985.txt,geneset_cmstructural.txt,geneset_ecm7.txt,geneset_kegg_oxphos_hsa00190.txt}"
  printf '%s\n$PY a0_primary_export.py spatial_first_config.json $COMP . pseudobulk_T2000.tsv.gz region_selection.tsv 2>&1 | tee run_log.txt\n# outputs: out/*\n' "$SLURM" ;;

local2)     # primary-arm export -> primary arm, locus QC, pooling, controls, A-scripts not needing per-spot data
  need_external
  if [ -n "${2:-}" ]; then O=$2; mkdir -p $DP; cp "$O"/* $DP/; fi
  run 01_prim $RS scripts/spatial_first/01_limma_q1q2.R --pb $DP/pseudobulk_allcm_T2000.tsv.gz --det $DP/detection_fraction_allcm_T2000.tsv.gz --meta $DP/meta_allcm_T2000.tsv --thresh 0.25 --tag _allcm --out $V
  for mode in section disease; do
    run 02_prim_$mode $RS scripts/spatial_first/02_leave_one_out.R --pb $DP/pseudobulk_allcm_T2000.tsv.gz --det $DP/detection_fraction_allcm_T2000.tsv.gz --meta $DP/meta_allcm_T2000.tsv --thresh 0.25 --mode $mode --out $V/loo
  done
  run 03_qc $PY scripts/spatial_first/03_locus_qc.py --universe $S/q1_raw_T2000_d25.tsv $V/q1_raw_allcm.tsv   # committed annotation caches; Ensembl REST only for genes in none of them
  run 04_sens $PY scripts/spatial_first/04_snrna_pool.py
  run 04_prim $PY scripts/spatial_first/04_snrna_pool.py --q1 $V/q1_raw_allcm.tsv --loo $V/loo --res $V
  run 05 $PY scripts/spatial_first/05_controls.py
  run 06 $PY scripts/spatial_first/06_published_reference.py
  run 00_inputs $PY scripts/eschf/00_inputs.py
  run derived_lists $PY scripts/eschf/derived_inputs.py lists
  run a0_refit $RS scripts/eschf/a0_primary_refit.R
  run a0_compare $PY scripts/eschf/a0_compare.py
  run a1 $RS scripts/eschf/a1_patient_refit.R
  for mode in old v2; do export ESCHF_PRIMARY=$mode
    run a3a4_$mode $PY scripts/eschf/a3a4_qc_nonfailing.py
    run a5_hamed_$mode $PY scripts/eschf/a5_hamed.py
    run a5b_proc_$mode $RS scripts/eschf/a5b_procurement.R
    run a5b_pool_$mode $PY scripts/eschf/a5b_pool.py
  done; unset ESCHF_PRIMARY
  run a9_means $RS scripts/eschf/a9_study_means.R
  run a11_raw $PY scripts/eschf/a11_raw_donor_check.py
  run a11_coh $PY scripts/eschf/a11_snrna_cohort.py
  run a10_cpm $RS scripts/eschf/a10_patient_log2cpm.R ;;

cluster3)   # per-spot program-score export (A2/A6)
  echo "# stage: scripts/eschf/a2a6_spot_export.py data/spatial_first/{spatial_first_config.json,pseudobulk_T2000.tsv.gz,region_selection.tsv} metadata/{tested_genes_985.txt,geneset_primary_program_qcpass_down.txt,geneset_cmstructural.txt,geneset_ecm7.txt,geneset_kegg_oxphos_hsa00190.txt}"
  printf '%s\n$PY a2a6_spot_export.py spatial_first_config.json $COMP . pseudobulk_T2000.tsv.gz region_selection.tsv 2>&1 | tee run_log.txt\n# outputs: out/* (images/ not needed locally)\n' "$SLURM" ;;

local3)     # per-spot export -> per-spot A-scripts, figure, summary, harness
  need_external
  if [ -n "${2:-}" ]; then O=$2; for f in a2_missing_genes.tsv eschf_export_meta.json eschf_perspot_T2000.tsv.gz; do cp "$O/$f" data/eschf/$f; done; fi
  run a0_perspot $PY scripts/eschf/a0_perspot_scores.py
  for mode in old v2; do export ESCHF_PRIMARY=$mode
    run a2a6_$mode $PY scripts/eschf/a2a6_analysis.py
    run a2_model_$mode $RS scripts/eschf/a2_model.R
  done
  export ESCHF_PRIMARY=v2; run a2_perpatient $PY scripts/eschf/a2_perpatient.py; unset ESCHF_PRIMARY
  run a9_screen $PY scripts/eschf/a9_snrna_screen.py
  run a10_levels $PY scripts/eschf/a10_pathway_levels.py
  run a14 $PY scripts/eschf/a14_pathway_hits_removed.py
  ESCHF_PRIMARY=v2 run a7 $PY scripts/eschf/a7_figure.py
  run a8 $PY scripts/eschf/a8_summary.py
  run harness $PY scripts/manuscript_results.py; tail -1 "$LOG/harness.log" ;;   # exits 1 on any MISMATCH

all)        # local1 -> local2 -> local3 from the committed exports (no cluster access needed)
  need_external
  for st in local1 local2 local3; do bash "$0" $st; done ;;
compare)
  $PY scripts/eschf/v3_old_vs_new.py "$2" . "$3" ;;

*) sed -n 2,22p "$0"; exit 2 ;;
esac
