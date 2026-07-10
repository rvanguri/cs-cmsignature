#!/usr/bin/env bash
#SBATCH --job-name=hm_qc
#SBATCH --partition=YOUR_PARTITION          # VERIFY. neyazi/reichart are big -> need headroom.
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                       # was 32G -> OOM on neyazi/reichart during QC
#SBATCH --time=04:00:00
#SBATCH --array=1-6                      # 6 datasets incl. Liu (5'); produces qc/liu_qc.h5ad for step 08
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
mkdir -p "${QC_OUT}"

DATASETS=(neyazi larson chin2022 chin2021 reichart liu)   # Liu QC'd here -> qc/liu_qc.h5ad (step 08 input)
DS=${DATASETS[$((SLURM_ARRAY_TASK_ID - 1))]}

# idempotent: skip a dataset already QC'd (lets a re-array redo only the OOM/failed ones)
if [[ -f "${QC_OUT}/${DS}_qc.h5ad" ]]; then
  echo "[$(date)] ${DS}: QC output already exists -> skipping"; exit 0
fi

# v2.1 QC verbatim: min_genes>=200, max_genes<=8000, MT%<5, CellBender cell-prob>0.5;
# scDblFinder/scrublet on decontaminated counts; genes in >=3 nuclei/sample; protein_coding+lncRNA.
python "${HERE}/../scripts/run_qc.py" \
    --dataset "${DS}" \
    --counts "${DECONTX_OUT}/${DS}_decontx_counts.h5ad" \
    --sample_meta "${RAW}/sample_meta.tsv" \
    --out "${QC_OUT}/${DS}_qc.h5ad"
echo "QC done: ${DS}"
