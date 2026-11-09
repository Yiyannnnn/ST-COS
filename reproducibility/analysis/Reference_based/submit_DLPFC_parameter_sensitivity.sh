#!/usr/bin/env bash
# Submit after the prepare branch has produced DLPFC_151673_sensitivity_input.rds.

set -euo pipefail

ROOT=.
SCRIPT=${ROOT}/Reference_based/DLPFC_parameter_sensitivity.R
OUT=${ROOT}/Reference_based/DLPFC_parameter_sensitivity
mkdir -p "${OUT}/logs"

submit_branch() {
  local branch=$1
  local mem=$2
  local time=$3
  sbatch --parsable -p general -N 1 -n 1 --mem="${mem}" -t "${time}" \
    -J "sens_${branch}" \
    -o "${OUT}/logs/${branch}_%j.out" \
    -e "${OUT}/logs/${branch}_%j.err" \
    --wrap="module load r/4.4.0; export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"; Rscript ${SCRIPT} ${branch}"
}

submit_branch matrix 24g 08:00:00
submit_branch dispersion 24g 06:00:00
submit_branch spatial_4 32g 24:00:00
submit_branch spatial_6 32g 24:00:00
submit_branch spatial_8 32g 24:00:00
submit_branch spatial_12 32g 24:00:00
