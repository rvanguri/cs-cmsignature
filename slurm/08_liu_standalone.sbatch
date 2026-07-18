#!/usr/bin/env bash
#SBATCH --job-name=cs_liu
#SBATCH --partition=YOUR_PARTITION           # VERIFY
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G                        # Liu ~57k cells; 8G OOM'd during marker typing + DESeq2
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"

# Liu CS-vs-ICM (5' v1.1) recomputed FRESH through QC+DecontX+pseudobulk+apeglm.
# Standalone: NOT co-embedded with the 3' arm. Only within-study, within-chemistry,
# explant-matched CS contrast -> co-primary. Report overlap (hypergeometric) + concordance heatmap.
# Reads the DECONTAMINATED + QC'd Liu produced by steps 03/04 (qc/liu_qc.h5ad), so Liu gets the SAME
# CellBender->DecontX->QC as the 3' sets: just not scANVI co-embedding.
python "${HERE}/../scripts/run_liu_standalone.py" \
    --raw "${QC_OUT}/liu_qc.h5ad" \
    --out "${DE_OUT}/liu_cs_vs_icm"
echo "Liu co-primary complete -> ${DE_OUT}/liu_cs_vs_icm"
