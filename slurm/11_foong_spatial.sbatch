#!/usr/bin/env bash
#SBATCH --job-name=hm_foong
#SBATCH --partition=YOUR_PARTITION           # VERIFY (download needs internet egress; data_mover if needed)
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"

# Foong cardiac-sarcoidosis Visium spatial data (GEO GSE314910: 12 CS myocardial Visium samples,
# CS-only). Download the RAW.tar of Space Ranger outputs, then render spatial maps of key CS marker
# programs + a TLS signature score per sample. (No CS-vs-HCM here: this series is CS-only; the HCM
# comparator lived in a separate Visium dataset used for the DE.)
export FOONG_GSE="${FOONG_GSE:-GSE314910}"
export FOONG_DIR="${FOONG_DIR:-${RAW}/foong_spatial}"

bash "${HERE}/../scripts/download_foong_spatial.sh"

python "${HERE}/../scripts/foong_spatial_figures.py" \
    --data "${FOONG_DIR}" \
    --out "${FIGURES}/foong_spatial" \
    ${FOONG_SAMPLE_MAP:+--sample_map "${FOONG_SAMPLE_MAP}"}

echo "Foong spatial figures complete -> ${FIGURES}/foong_spatial"
