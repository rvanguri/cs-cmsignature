#!/usr/bin/env bash
# Launch the full pipeline as a Slurm dependency chain. Open-data arms run end-to-end;
# controlled arms fold in later (see env/access_status.md).
#
#   cd $PROJ/slurm && bash submit_all.sh
#
# Each stage waits for clean (afterok) completion of its prerequisite. CellBender (02) is a
# job array; downstream waits for the WHOLE array via afterok on the array job id.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${HERE}"

# Send logs to $PROJ/logs
source "${HERE}/../env/paths.sh"
mkdir -p "${LOGS}"
LOGDIR_FLAG="--chdir=${LOGS}"   # %x_%j.out files land in $LOGS; remove if you prefer cwd

jid1=$(sbatch --parsable ${LOGDIR_FLAG} 01_download.sbatch)
echo "01 download        : $jid1"
jid2=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid1 02_cellbender.sbatch)
echo "02 cellbender(array): $jid2"
jid3=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid2 03_decontx.sbatch)
echo "03 decontx         : $jid3"
jid4=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid3 04_qc.sbatch)
echo "04 qc              : $jid4"
jid5=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid4 05_scanvi.sbatch)
echo "05 scanvi          : $jid5"
jid6=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid5 06_validate_gates.sbatch)
echo "06 gates           : $jid6"
# 07, 08 both depend on the gates (06) and run in parallel:
jid7=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid6 07_pseudobulk_de.sbatch)
echo "07 pseudobulk DE   : $jid7"
jid8=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid6 08_liu_standalone.sbatch)
echo "08 liu standalone  : $jid8"
# 10 waits for both:
jid10=$(sbatch --parsable ${LOGDIR_FLAG} --dependency=afterok:$jid7:$jid8 10_gsea_figures.sbatch)
echo "10 gsea + figures  : $jid10"

echo
echo "Submitted. Monitor: squeue -u \$USER   |   sacct -j <jid> --format=JobID,State,Elapsed,MaxRSS"
echo "NOTE: edit 02_cellbender.sbatch --array=1-N%4 so N == manifest sample count before launch."
