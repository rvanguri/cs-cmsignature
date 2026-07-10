#!/usr/bin/env bash
#SBATCH --job-name=hm_gates
#SBATCH --partition=YOUR_PARTITION          # VERIFY
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                       # ~1M-cell integrated object; 32G risks OOM on load + 6a
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
mkdir -p "${DE_OUT}" "${FIGURES}"

# 6a integration validation: UMAP by disease/study/celltype, batch mixing, kBET/LISI (scib).
# 6b DECONTAMINATION QC GATE incl. dual-count concordance (CellBender/DecontX vs CellRanger raw)
#    + background heuristic <0.5. Corrected expectations: TNNI3K/GJB7 persist;
#    SPECC1L/CALM2/TJP2 NOT expected; CM-dediff panel = technical control only.
# 6c PROCUREMENT-ROBUSTNESS GATE across explant-matched mimics (DCM, ICM-Simonson, ARVC, Chaffin-LV).
python "${HERE}/../scripts/run_gates.py" \
    --integrated "${INTEGRATED}/integrated_scanvi.h5ad" \
    --cellbender_dir "${CB_OUT}" \
    --out "${DE_OUT}/gates" \
    --figdir "${FIGURES}"
echo "Gates complete. REVIEW gate outputs before trusting any DEG."
