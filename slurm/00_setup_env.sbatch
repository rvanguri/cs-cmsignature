#!/usr/bin/env bash
#SBATCH --job-name=hm_setup
#SBATCH --partition=YOUR_PARTITION          # VERIFY: sinfo -s
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# One-time environment provisioning. Run on a node with internet egress
# (or run the conda/singularity-pull lines manually on a login/data-mover node).
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
module purge
module load singularity/3.9.8          # VERIFY
module load miniconda3                 # VERIFY

# 1) conda env on persistent storage
source "$(conda info --base)/etc/profile.d/conda.sh"
if [[ ! -d "${CONDA_ENV}" ]]; then
    mamba env create -f "${HERE}/../env/environment.yml" -p "${CONDA_ENV}"
fi

# 2) CellBender Singularity image
if [[ ! -f "${CELLBENDER_SIF}" ]]; then
    singularity pull "${CELLBENDER_SIF}" \
        docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0
fi

echo "Env ready: ${CONDA_ENV} + ${CELLBENDER_SIF}"
