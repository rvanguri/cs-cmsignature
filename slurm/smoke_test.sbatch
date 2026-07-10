#!/usr/bin/env bash
#SBATCH --job-name=hm_smoke
#SBATCH --partition=YOUR_PARTITION            # VERIFY. No GPU needed for the default smoke chain.
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Dry-run the WHOLE pipeline on ~450 synthetic nuclei before real data lands. Validates the scripts and
# their file handoffs end-to-end. Writes everything under $SCRATCH/smoke (never touches real $PROJ).
#
# By default it exercises: synth-data -> build_manifest -> QC -> scANVI(2 epochs) -> gates -> DE ->
# Liu -> projection(skips, no atlas) -> GSEA  — all REAL script calls.
# The GPU CellBender leg and the R DecontX leg are validated separately (opt-in):
#   SMOKE_DECONTX=1   sbatch slurm/smoke_test.sbatch     # also runs run_decontx.R on synthetic raw
#   SMOKE_CELLBENDER=1 sbatch slurm/smoke_test.sbatch    # also runs the CellBender .sif on 1 sample (needs GPU)
#
set -uo pipefail   # NOTE: not -e; we want to record per-step PASS/FAIL and keep going
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
SC="${HERE}/../scripts"

SMOKE="${SCRATCH_ROOT}/smoke"
rm -rf "${SMOKE}"; mkdir -p "${SMOKE}"/{qc,integrated,de/tables,de/gates,figures,heartmap/projection,cellbender,scanvi_work}

PASS=0; FAIL=0
step () {  # step "name" command...
  local name="$1"; shift
  echo "----- [SMOKE] ${name} -----"
  if "$@"; then echo "PASS: ${name}"; PASS=$((PASS+1));
  else echo "FAIL: ${name}"; FAIL=$((FAIL+1)); fi
}

# 0) generate synthetic data
step "make_synthetic_data" python "${SC}/make_synthetic_data.py" --out "${SMOKE}"

# 1) build_manifest (real, over synthetic raw 10x mtx — exercises knee estimate)
step "build_manifest" python "${SC}/build_manifest.py" --raw "${SMOKE}/raw" --out "${SMOKE}/raw/manifest.tsv"

# 1b) OPTIONAL: R DecontX leg on synthetic raw (validates celda wiring)
if [[ "${SMOKE_DECONTX:-0}" == "1" ]]; then
  for ds in neyazi larson chin2022 chin2021 reichart; do
    step "decontx_${ds}" Rscript "${SC}/run_decontx.R" --dataset "${ds}" \
      --cellbender_dir "${SMOKE}/cellbender" --reichart_dir "${SMOKE}/raw/reichart" \
      --out "${SMOKE}/decontx/${ds}_decontx_counts.h5ad"
  done
fi
# (default: generator already wrote ${SMOKE}/decontx/*_decontx_counts.h5ad)

# 1c) OPTIONAL: CellBender GPU leg on ONE synthetic sample (needs --partition gpu* + --gres; run separately)
if [[ "${SMOKE_CELLBENDER:-0}" == "1" ]]; then
  SAMP=$(sed -n '2p' "${SMOKE}/raw/manifest.tsv"); RP=$(echo "$SAMP" | cut -f7)
  step "cellbender_one" singularity exec --nv -B "${SCRATCH_ROOT}" "${CELLBENDER_SIF}" \
    cellbender remove-background --cuda --input "${RP}" \
    --output "${SMOKE}/cellbender/smoke_cellbender.h5" --expected-cells 30 \
    --total-droplets-included 1500 --fpr 0.01 --epochs 20
fi

# 2) QC on each dataset (thresholds relaxed for tiny synthetic cells)
for ds in neyazi larson chin2022 chin2021 reichart; do
  step "qc_${ds}" python "${SC}/run_qc.py" --dataset "${ds}" \
    --counts "${SMOKE}/decontx/${ds}_decontx_counts.h5ad" \
    --out "${SMOKE}/qc/${ds}_qc.h5ad" --min_genes 5 --max_genes 100000 --max_mt 100
done

# 3) scANVI integration (2 epochs)
step "scanvi" python "${SC}/run_scanvi.py" --qc_dir "${SMOKE}/qc" --reichart_labels reichart \
  --workdir "${SMOKE}/scanvi_work" --out_model "${SMOKE}/integrated/scanvi_model" \
  --out_h5ad "${SMOKE}/integrated/integrated_scanvi.h5ad" --max_epochs 2

# NOTE: a 2-epoch scANVI on ~230 toy cells can't learn real predictions, so 'predicted_cell_type'
# has no cardiomyocytes. For the DRY RUN we point gates+DE at the ground-truth 'cell_type' column the
# synthetic data carries. The real pipeline (06/07 sbatch) keeps the default 'predicted_cell_type'.

# 4) gates
step "gates" python "${SC}/run_gates.py" --integrated "${SMOKE}/integrated/integrated_scanvi.h5ad" \
  --cellbender_dir "${SMOKE}/cellbender" --out "${SMOKE}/de/gates" --figdir "${SMOKE}/figures" \
  --celltype_col cell_type

# 5) pseudobulk DE
step "pseudobulk_de" Rscript "${SC}/run_pseudobulk_de.R" \
  --integrated "${SMOKE}/integrated/integrated_scanvi.h5ad" --gates "${SMOKE}/de/gates" \
  --out "${SMOKE}/de" --celltype_col cell_type

# 6) Liu standalone co-primary
step "liu_standalone" python "${SC}/run_liu_standalone.py" --raw "${SMOKE}/raw/liu" \
  --out "${SMOKE}/de/liu_cs_vs_icm"

# 7) HeartMap projection (non-fatal: no atlas in smoke -> exits 0)
step "heartmap_projection" python "${SC}/run_heartmap_projection.py" \
  --integrated "${SMOKE}/integrated/integrated_scanvi.h5ad" \
  --heartmap_atlas "${SMOKE}/heartmap/atlas" --out "${SMOKE}/heartmap/projection"

# 8) GSEA + figures + high-confidence set
step "gsea_figures" python "${SC}/run_gsea_figures.py" --de "${SMOKE}/de" \
  --liu "${SMOKE}/de/liu_cs_vs_icm" --projection "${SMOKE}/heartmap/projection" \
  --figdir "${SMOKE}/figures" --tabledir "${SMOKE}/de/tables" --fdr 0.5 --lfc 0.1

echo
echo "================ SMOKE SUMMARY ================"
echo "PASS=${PASS}  FAIL=${FAIL}   outputs under ${SMOKE}"
echo "Key artifacts to eyeball:"
echo "  ${SMOKE}/raw/manifest.tsv"
echo "  ${SMOKE}/integrated/integrated_scanvi.h5ad"
echo "  ${SMOKE}/de/DE_CS_vs_*.tsv  ${SMOKE}/de/liu_cs_vs_icm/DE_Liu_CS_vs_ICM.tsv"
echo "  ${SMOKE}/de/tables/high_confidence_cs_set.tsv  ${SMOKE}/de/gates/decontam_qc_gate.tsv"
echo "==============================================="
[[ "${FAIL}" -eq 0 ]] && exit 0 || exit 1
