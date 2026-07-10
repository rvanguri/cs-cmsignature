#!/usr/bin/env bash
#SBATCH --job-name=hm_cmsig
#SBATCH --partition=YOUR_PARTITION           # VERIFY
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"

# CM-core specificity analyses.
# (1) Per-patient CM pseudobulk disease-distance: does CS cluster apart from other CMs?
python "${HERE}/../scripts/cm_disease_distance.py" \
    --pb_counts "${DE_OUT}/cm_pseudobulk_counts.tsv" \
    --pb_meta   "${DE_OUT}/cm_pseudobulk_meta.tsv" \
    --out "${FIGURES}/cm_distance"

# (2) CS-CM signature in CM-dominant Visium spots across diseases.
# Manifest (condition<TAB>path). CS/Foong is ready now; add control (Kuppe) / HCM rows as you stage
# them. Default manifest lists just the Foong CS data.
CM_MANIFEST="${CM_MANIFEST:-${PROJ}/external/visium_manifest.tsv}"
if [[ ! -f "${CM_MANIFEST}" ]]; then
    mkdir -p "$(dirname "${CM_MANIFEST}")"
    printf 'CS\t%s\n' "${RAW}/foong_spatial" > "${CM_MANIFEST}"
    echo "[cmsig] wrote default manifest (CS only) -> ${CM_MANIFEST}; add control/HCM rows to compare"
fi
python "${HERE}/../scripts/cm_spatial_crossdisease.py" \
    --manifest "${CM_MANIFEST}" \
    --out "${FIGURES}/cm_spatial_crossdisease"

# (3) CM-content-normalized per-gene comparison (removes the CM density/purity confound).
python "${HERE}/../scripts/cm_content_normalized.py" \
    --manifest "${CM_MANIFEST}" \
    --out "${FIGURES}/cm_content_normalized"

# (4) Within-Foong regional analysis (no cross-dataset batch): signatures by CS tissue zone +
#     CM-intrinsic distance-to-lesion curve (cell-autonomous vs paracrine).
python "${HERE}/../scripts/foong_regional.py" \
    --data "${RAW}/foong_spatial" \
    --out "${FIGURES}/foong_regional"

echo "CM-signature analyses complete -> ${FIGURES}/{cm_distance,cm_spatial_crossdisease,cm_content_normalized,foong_regional}"
