#!/usr/bin/env bash
#SBATCH --job-name=cs_cellbender
#SBATCH --partition=YOUR_PARTITION         # VERIFY name + MaxTime (sinfo -s). gpu8_* if you have it.
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-70%4                   # set N = total manifest rows (filtered samples self-skip); %4 = 4 concurrent GPU tasks
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"      # loads singularity

# --- pick this task's sample from the manifest (header on line 1) ---
ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
LINE=$(sed -n "${ROW}p" "${MANIFEST}")
# manifest cols: 1=sample_id 2=dataset 3=accession 4=chemistry 5=procurement 6=region 7=raw_path
#                8=expected_cells_seed 9=n_barcodes 10=cellbender 11=total_droplets
SAMPLE_ID=$(echo "${LINE}" | cut -f1)
RAW_PATH=$(echo  "${LINE}" | cut -f7)
EXP_CELLS=$(echo "${LINE}" | cut -f8)
CELLBENDER=$(echo "${LINE}" | cut -f10)   # 1 = raw (run CellBender) ; 0 = filtered (skip -> DecontX-only)
TOT_DROP=$(echo  "${LINE}" | cut -f11)    # adaptive droplet cap (<= n_barcodes)

# --- filtered samples (too few barcodes for ambient modeling) skip CellBender; DecontX handles them ---
if [[ "${CELLBENDER}" != "1" ]]; then
  echo "[$(date)] ${SAMPLE_ID}: filtered matrix (cellbender=0) -> skipping CellBender, DecontX-only"
  exit 0
fi

# --- idempotent: if this sample already produced output, skip (lets a re-array redo only failures) ---
if [[ -f "${CB_OUT}/${SAMPLE_ID}_cellbender_filtered.h5" ]]; then
  echo "[$(date)] ${SAMPLE_ID}: output already exists in ${CB_OUT} -> skipping (already done)"
  exit 0
fi

# Fresh per-sample working dir. CellBender writes its checkpoint (ckpt.tar.gz) to the WORKING
# DIRECTORY, so concurrent array tasks MUST NOT share a cwd or they race/overwrite each other's
# checkpoint ("Workflow hash does not match" / "ckpt.tar.gz.tmp -> ckpt.tar.gz" FileNotFound).
WORK="${CB_WORK}/${SAMPLE_ID}"
rm -rf "${WORK}"; mkdir -p "${WORK}" "${CB_OUT}"     # rm clears any stale checkpoint from a prior attempt
cd "${WORK}"
echo "[$(date)] CellBender ${SAMPLE_ID}  expected_cells=${EXP_CELLS}  total_droplets=${TOT_DROP}  cwd=${WORK}  src=${RAW_PATH}"

# --- resolve symlinks: $RAW (or a per-dataset subdir) may be symlinked to shared storage that is
#     OUTSIDE the container's default binds. Singularity only sees paths you -B mount, and it does
#     not follow a host symlink that escapes a mounted tree. Resolve RAW_PATH to its real location
#     and bind that real directory, so the matrix is visible inside the container regardless of where
#     the raw data physically lives. ---
REAL_RAW_PATH=$(readlink -f "${RAW_PATH}")
REAL_RAW_DIR=$(dirname "${REAL_RAW_PATH}")

# --- identical params across all RAW samples; --pwd pins the container cwd to this sample's dir ---
singularity exec --nv --pwd "${WORK}" -B "${SCRATCH_ROOT}","${PROJ}","${REAL_RAW_DIR}" "${CELLBENDER_SIF}" \
    cellbender remove-background --cuda \
      --input  "${REAL_RAW_PATH}" \
      --output "${WORK}/${SAMPLE_ID}_cellbender.h5" \
      --expected-cells "${EXP_CELLS}" \
      --total-droplets-included "${TOT_DROP}" \
      --fpr 0.01 --epochs 150

# --- persist durable output off scratch immediately ---
rsync -a "${WORK}/${SAMPLE_ID}_cellbender"* "${CB_OUT}/"
echo "[$(date)] persisted ${SAMPLE_ID} -> ${CB_OUT}"
