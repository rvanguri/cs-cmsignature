#!/usr/bin/env bash
#SBATCH --job-name=hm_relabel
#SBATCH --partition=YOUR_PARTITION          # VERIFY
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"

# Relabel Ensembl var_names -> HGNC symbols on the objects consumed by symbol-space
# downstream code (gates/DE/gsea/foong/cm read integrated; liu_standalone reads liu_qc).
# The decontx/QC steps carry Ensembl IDs; the original run applied this map before gates.
# relabel_genes.py stashes the old Ensembl id in var['ensembl'] and disambiguates the
# handful of duplicate symbols, so the step is idempotent and safe to re-run.
FEAT="${RAW}/liu/GSM6222406_TWCM-11-290_features.tsv.gz"
echo "Using features map: ${FEAT}"

python "${HERE}/../scripts/relabel_genes.py" \
    --h5ad "${INTEGRATED}/integrated_scanvi.h5ad" \
    --features "${FEAT}"

python "${HERE}/../scripts/relabel_genes.py" \
    --h5ad "${QC_OUT}/liu_qc.h5ad" \
    --features "${FEAT}"

echo "Relabel complete -> integrated_scanvi.h5ad + liu_qc.h5ad now in symbol space"
