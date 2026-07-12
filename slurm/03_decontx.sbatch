#!/usr/bin/env bash
#SBATCH --job-name=hm_decontx
#SBATCH --partition=YOUR_PARTITION          # VERIFY. Reichart (880k nuclei) + Neyazi are memory-heavy.
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G                       # Reichart (880k nuclei) peaks here: DecontX itself is cheap, but
                                         # the final t()+reticulate write_h5ad holds transient R+Python
                                         # copies of the ~880k-cell sparse matrix. 128G OOM'd on the write;
                                         # 256G clears it (fits on a standard cpu_short node, 320G/node).
#SBATCH --time=08:00:00
#SBATCH --array=1-6                      # 6 datasets: Neyazi, Larson, Chin2022, Chin2021, Reichart, Liu(5')
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"

mkdir -p "${DECONTX_OUT}" "${DECONTX_WORK}"

# UNIFORM DecontX pass — identical params on the FILTERED matrix of every dataset:
#   - 4 GEO sets: DecontX on the CellBender-filtered matrix ($CB_OUT)
#   - Reichart:   DecontX on the CELLxGENE filtered matrix
# Output rounded to integers for pseudobulk. This makes the LAST decontamination op common to all.
# Liu (5') gets the SAME decontamination as the 3' sets (it's CellBendered via the manifest array);
# it stays out of scANVI co-embedding (chemistry), but its standalone DE needs identical decontam.
DATASETS=(neyazi larson chin2022 chin2021 reichart liu)
DS=${DATASETS[$((SLURM_ARRAY_TASK_ID - 1))]}

# idempotent: skip a dataset that already produced its DecontX output (lets a re-array redo only failures)
if [[ -f "${DECONTX_OUT}/${DS}_decontx_counts.h5ad" ]]; then
  echo "[$(date)] ${DS}: decontx output already exists -> skipping"
  exit 0
fi

Rscript "${HERE}/../scripts/run_decontx.R" \
    --dataset "${DS}" \
    --manifest "${MANIFEST}" \
    --cellbender_dir "${CB_OUT}" \
    --reichart_dir "${RAW}/reichart" \
    --out "${DECONTX_OUT}/${DS}_decontx_counts.h5ad"
echo "DecontX done: ${DS}"
