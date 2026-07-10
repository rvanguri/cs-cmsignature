#!/usr/bin/env bash
# Central path + identity config. Source this at the top of every sbatch script:
#   source "$(dirname "$0")/../env/paths.sh"
#
# EDIT THESE TWO LINES for your account:
LAB="CHANGEME_lab"          # your /gpfs/data/<lab> group
KID="${USER}"               # KerberosID (defaults to $USER)

# ---- persistent (deliverables; survives scratch purge) ----
export PROJ="/gpfs/data/${LAB}/heartmap-cs"
# ---- fast scratch (intermediates; PURGED periodically) ----
export SCRATCH_ROOT="/gpfs/scratch/${KID}/heartmap-cs"
# ---- conda env + singularity images ----
export CONDA_ENV="${PROJ}/env/conda/heartmap"
export CELLBENDER_SIF="${PROJ}/env/cellbender.sif"

# ---- subdirs (created by scaffold.sh) ----
export RAW="${PROJ}/raw"
export CONTROLLED="${PROJ}/controlled"
export CB_OUT="${PROJ}/cellbender"
export DECONTX_OUT="${PROJ}/decontx"
export QC_OUT="${PROJ}/qc"
export INTEGRATED="${PROJ}/integrated"
export DE_OUT="${PROJ}/de"
export HEARTMAP="${PROJ}/heartmap"
export FIGURES="${PROJ}/figures"
export LOGS="${PROJ}/logs"
export MANIFEST="${RAW}/manifest.tsv"

# ---- scratch work dirs ----
export CB_WORK="${SCRATCH_ROOT}/cellbender_work"
export DECONTX_WORK="${SCRATCH_ROOT}/decontx_work"
export SCANVI_WORK="${SCRATCH_ROOT}/scanvi_work"

# Singularity cache on home (avoid scratch purge of layers)
export SINGULARITY_CACHEDIR="${HOME}/.singularity"
