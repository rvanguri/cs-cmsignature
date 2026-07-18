#!/usr/bin/env bash
#SBATCH --job-name=cs_scanvi
#SBATCH --partition=YOUR_PARTITION           # VERIFY name/limits (sinfo -s). Needs a long GPU queue + big RAM.
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G                       # >1M nuclei (Reichart 880k dominates); 64G OOM'd on concat.
                                         # If the GPU partition can't give 192G, drop --gres and use a
                                         # high-mem CPU node (fn_*, --mem=256G): slower but works.
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
mkdir -p "${INTEGRATED}" "${SCANVI_WORK}"

# Steps 3+4+5: gene-ID harmonization, feature intersection (5 x 3' datasets), then scANVI.
#   .X = post-DecontX counts ; batch_key=study ; labels_key=cell_type_harmonized
#   unlabeled_category="Unknown" ; B-cell split from T/NK ; Liu (5') NOT co-embedded.
#   model.train(max_epochs=400, early_stopping=True). scvi-tools uses GPU if visible.
python "${HERE}/../scripts/run_scanvi.py" \
    --qc_dir "${QC_OUT}" \
    --reichart_labels reichart \
    --workdir "${SCANVI_WORK}" \
    --out_model "${INTEGRATED}/scanvi_model" \
    --out_h5ad "${INTEGRATED}/integrated_scanvi.h5ad"
echo "scANVI integration complete -> ${INTEGRATED}"
