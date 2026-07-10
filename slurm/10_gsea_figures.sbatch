#!/usr/bin/env bash
#SBATCH --job-name=hm_gsea_fig
#SBATCH --partition=YOUR_PARTITION           # VERIFY
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                       # 8G OOM'd rendering per-term GSEA plots; no_plot=True now
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${HERE}/../env/paths.sh"
source "${HERE}/../env/modules.sh"
mkdir -p "${FIGURES}"

# Pre-ranked gseapy (Hallmark, GO BP, KEGG, Reactome) on procurement-matched primaries +
# high-confidence set (NOT the pooled contrast as headline). Then all v2.1 figures + tables:
# integration UMAPs, CM UMAP, volcanoes, DEG heatmap, GSEA dot plots, decontamination QC-gate fig,
# procurement-robustness fig (TNNI3K/GJB7 highlighted; CM-dediff panel labeled non-specific control),
# HeartMap positioning fig, Liu concordance heatmap, patient/region summary.
python "${HERE}/../scripts/run_gsea_figures.py" \
    --de "${DE_OUT}" \
    --liu "${DE_OUT}/liu_cs_vs_icm" \
    --figdir "${FIGURES}" \
    --tabledir "${DE_OUT}/tables"
echo "GSEA + figures complete -> ${FIGURES}"

# External CM validation: combine ALL CS-vs-mimic CM contrasts with independent Foong (Visium)
# replication -> de_cm/CM_foong_validation.tsv. Requires the Foong significant-DEG csv staged at $FOONG
# (this is the CS-vs-HCM spatial DE, i.e. the uploaded SIGDEG table), and the per-lineage de_cm/ dir.
FOONG="${FOONG:-${PROJ}/external/foong_cs_vs_hcm.csv}"
if [[ -f "${FOONG}" && -d "${DE_OUT}/de_cm" ]]; then
    python "${HERE}/../scripts/annotate_foong_validation.py" \
        --de_dir "${DE_OUT}/de_cm" --foong "${FOONG}" \
        --out "${DE_OUT}/de_cm/CM_foong_validation.tsv"
    # panel robustness sweep (Concern 4.3) reads the same de_cm tables + Foong csv
    python "${HERE}/../scripts/panel_robustness.py" \
        --de_dir "${DE_OUT}/de_cm" --foong "${FOONG}" \
        --out "${DE_OUT}/de_cm/panel_robustness.tsv" || true
else
    echo "Foong csv not found at ${FOONG} (or de_cm missing) — skipping external CM validation."
    echo "Stage it:  mkdir -p ${PROJ}/external && cp <SIGDEG.csv> ${FOONG}"
fi
