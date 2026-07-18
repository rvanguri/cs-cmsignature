#!/usr/bin/env bash
# Lmod modules + conda activation, sourced by every sbatch script so the environment is identical
# everywhere. env/paths.sh must be sourced BEFORE this file (it exports PROJ, CONDA_ENV).
#
# The system conda on the cluster may be ancient (e.g. 4.13, no libmamba) and can't solve this env, so we install
# a project-local Miniforge (see PLAN / setup steps) and source THAT here. Falls back to a system
# conda module if Miniforge isn't present.

module purge 2>/dev/null || true

# --- Singularity (for CellBender). VERIFY the version with: module avail singularity ---
module load singularity/3.9.8 2>/dev/null || module load singularity 2>/dev/null || true

# --- conda: prefer project Miniforge, else a system module ---
if [[ -f "${PROJ}/env/miniforge/etc/profile.d/conda.sh" ]]; then
    source "${PROJ}/env/miniforge/etc/profile.d/conda.sh"
else
    module load miniconda3 2>/dev/null || module load anaconda3 2>/dev/null || true
    if command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
    fi
fi

# --- activate the project env if it exists ---
if [[ -n "${CONDA_ENV:-}" && -d "${CONDA_ENV}" ]]; then
    conda activate "${CONDA_ENV}"
fi

# R 'anndata' pkg reads/writes .h5ad through reticulate -> use the env's python (not a stray one)
[[ -n "${CONDA_PREFIX:-}" ]] && export RETICULATE_PYTHON="${CONDA_PREFIX}/bin/python"

# flush python stdout to the Slurm log in real time (otherwise progress is invisible until the job ends)
export PYTHONUNBUFFERED=1
