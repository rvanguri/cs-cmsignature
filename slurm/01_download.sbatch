#!/usr/bin/env bash
#SBATCH --job-name=hm_download
#SBATCH --partition=YOUR_PARTITION           # data_mover if available for large pulls; VERIFY
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"

mkdir -p "${RAW}"

# --- OPEN GEO datasets (no access gate) ---
#   Neyazi CS GSE319770 (41 raw_feature_bc_matrix H5)
#   Larson HCM GSE174691 ; Chin2022 GSE181764 ; Chin2021 GSE161921
# Use your preferred GEO fetch (wget on supplementary tarballs, or geofetch). Example:
python "${HERE}/../scripts/download_geo.py" --manifest "${RAW}/accessions.tsv" --out "${RAW}"

# --- Reichart CELLxGENE matrix (open) + provenance check ---
python "${HERE}/../scripts/download_reichart_cxg.py" --out "${RAW}/reichart"
#   -> writes ${RAW}/reichart_provenance.md : is the matrix CellBender output or other? DO NOT ASSUME.

# --- CONTROLLED arms (only if env/access_status.md says raw-reprocess) ---
#   dbGaP phs001539 via SRA Toolkit with .ngc key, into $CONTROLLED (restricted ACLs).
#   Run on data-mover nodes; coordinate with your HPC support team first. NOT executed by default.

# --- Build the sample manifest that every array job indexes into ---
python "${HERE}/../scripts/build_manifest.py" --raw "${RAW}" --out "${MANIFEST}"
echo "Manifest rows (samples w/ raw droplets): $(($(wc -l < "${MANIFEST}") - 1))"
