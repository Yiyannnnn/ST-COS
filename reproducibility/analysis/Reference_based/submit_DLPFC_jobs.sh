#!/usr/bin/env bash
## ---------------------------------------------------------------
## Submit one SLURM job per DLPFC slide to run ST-COS, SRTsim,
## and scDesign3 for all expressed non-MT genes.
##
## Usage:
##   bash submit_DLPFC_jobs.sh              # submit all 12 slides
##   bash submit_DLPFC_jobs.sh 151673       # submit a single slide
## ---------------------------------------------------------------

set -euo pipefail

RSCRIPT="./Reference_based/DLPFC_run_methods_by_slide.R"
OUT_DIR="./Reference_based/DLPFC"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${LOG_DIR}"

ALL_SLIDES=(151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676)

# Determine which slides to run
if [[ $# -ge 1 ]]; then
  SLIDES=("$@")
else
  SLIDES=("${ALL_SLIDES[@]}")
fi

echo "[submit] Submitting ${#SLIDES[@]} job(s)..."

for SLIDE in "${SLIDES[@]}"; do
  JOB_NAME="DLPFC_${SLIDE}"
  LOG_OUT="${LOG_DIR}/${SLIDE}_%j.out"
  LOG_ERR="${LOG_DIR}/${SLIDE}_%j.err"

  JOB_ID=$(sbatch \
    --job-name="${JOB_NAME}" \
    --partition=general \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=128G \
    --time=60:00:00 \
    --output="${LOG_OUT}" \
    --error="${LOG_ERR}" \
    --wrap="module load r/4.4.0 && \
            export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}" && \
            echo \"[$(date)] Starting slide ${SLIDE}\" && \
            Rscript ${RSCRIPT} ${SLIDE} && \
            echo \"[$(date)] Finished slide ${SLIDE}\"" \
    | awk '{print $NF}')

  echo "[submit] Slide ${SLIDE} -> job ${JOB_ID}  (log: ${LOG_DIR}/${SLIDE}_${JOB_ID}.out)"
done

echo ""
echo "[submit] Monitor with:  squeue -u \$USER"
echo "[submit] Output files:  ${OUT_DIR}/<slide_id>_STCOS_simulated_counts.csv"
echo "                        ${OUT_DIR}/<slide_id>_SRTsim_counts.csv"
echo "                        ${OUT_DIR}/<slide_id>_scDesign3_counts.csv"
echo "                        ${OUT_DIR}/<slide_id>_simulated_coordinates.csv"
echo "                        ${OUT_DIR}/<slide_id>_timing.csv"
