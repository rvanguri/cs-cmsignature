#!/usr/bin/env bash
#SBATCH --job-name=hm_de
#SBATCH --partition=YOUR_PARTITION          # VERIFY
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                       # loads the ~1.1M-cell integrated h5ad in R; 16G OOM'd
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
mkdir -p "${DE_OUT}"

# Extract CMs; pseudobulk per patient (sum -> integer).
# HeartMap DE model: Expr ~ 0 + disease + study + anatomy + sex  with limma-voom +
#   duplicateCorrelation(individual). apeglm shrinkage + dupCor + >=75% concordance + procurement gate.
# Primary: CS-vs-DCM, CS-vs-ICM(Simonson 3'), CS-vs-ARVC. Supporting: HCM(Chaffin LV / myectomy), NF, pooled.
# Filter panels: IEG/HSP + HeartMap cross-study DCM generic-HF signature.
# scANVI lumps all lymphocytes into one 'TNK' label (no separate B compartment). Split it into
# clean T/NK vs B by marker score FIRST, so the T-cell DE isn't contaminated by B cells (and vice
# versa). Writes an obs_name->{Tcell,Bcell} map that the pseudobulk step applies via --subtype_tsv.
SUBTSV="${DE_OUT}/tnk_subtype.tsv"
python "${HERE}/../scripts/subtype_tnk.py" \
    --integrated "${INTEGRATED}/integrated_scanvi.h5ad" --out "${SUBTSV}" \
  || echo "subtype_tnk failed — T/B lineages will find no cells and skip"

# Run the SAME pseudobulk -> limma-voom DE for each lineage. Format: "label:regex".
#   cm    = cardiomyocytes (headline)
#   tcell = T/NK cells (refined label 'Tcell' from subtype_tnk)
#   bcell = B/plasma cells (refined label 'Bcell' from subtype_tnk)
# The Python backed-read only materializes the matching lineage, so T/B passes are cheap.
LINEAGES=( "cm:ardiomyocyte" "tcell:Tcell" "bcell:Bcell" )
for spec in "${LINEAGES[@]}"; do
    label="${spec%%:*}"; pat="${spec#*:}"
    echo "===== lineage ${label}  (/${pat}/) ====="
    # Step 1 (Python, backed read): aggregate the lineage -> small genes x patients table.
    python "${HERE}/../scripts/pseudobulk_cm.py" \
        --integrated "${INTEGRATED}/integrated_scanvi.h5ad" \
        --out "${DE_OUT}" --label "${label}" --celltype_pattern "${pat}" --min_cells 10 \
        --subtype_tsv "${SUBTSV}" \
      || { echo "[${label}] pseudobulk failed (too few cells?) — skipping"; continue; }
    # Step 2 (R): limma-voom DE on the small pseudobulk table, into a per-lineage subdir.
    mkdir -p "${DE_OUT}/de_${label}"
    Rscript "${HERE}/../scripts/run_pseudobulk_de.R" \
        --pb_counts "${DE_OUT}/${label}_pseudobulk_counts.tsv" \
        --pb_meta   "${DE_OUT}/${label}_pseudobulk_meta.tsv" \
        --out "${DE_OUT}/de_${label}" \
      || echo "[${label}] DE failed — continuing"
done
# Keep the CM results at the top level too (back-compat: steps 08-10 read ${DE_OUT}/DE_CS_vs_*.tsv).
cp -f "${DE_OUT}/de_cm/"DE_*.tsv "${DE_OUT}/" 2>/dev/null || true
echo "Pseudobulk DE complete (cm/tcell/bcell) -> ${DE_OUT}"
