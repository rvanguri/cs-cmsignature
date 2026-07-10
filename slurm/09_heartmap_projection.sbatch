#!/usr/bin/env bash
#SBATCH --job-name=hm_project
#SBATCH --partition=YOUR_PARTITION          # VERIFY. CPU fallback: cpu_medium, drop --gres
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
mkdir -p "${HEARTMAP}/projection"

# Project CS CM/macrophage/fibroblast onto HeartMap SCP3689 (scArches/scANVI transfer).
# NO raw co-embed (avoids decontamination mismatch) — transfer only.
# Report CS position vs DCM/HCM/ICM/ARVC/control; COL22A1+ (DCM) vs TNC+ (ICM) fibroblast niches.
python "${HERE}/../scripts/run_heartmap_projection.py" \
    --integrated "${INTEGRATED}/integrated_scanvi.h5ad" \
    --heartmap_atlas "${HEARTMAP}/atlas" \
    --out "${HEARTMAP}/projection"
echo "HeartMap projection complete -> ${HEARTMAP}/projection"
