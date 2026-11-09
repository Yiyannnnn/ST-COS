#!/usr/bin/env bash
set -euo pipefail

repro_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_dir="${repro_dir}/analysis"
cd "${analysis_dir}"

scripts=(
  "Data/S1_gene/simulated/S1_gene_simulate.R"
  "Data/S2_gene/simulated/S2_ctcomp_simulate.R"
  "Data/S3_gene/simulated/S3_gene_ctcomp_simulate.R"
  "Data/S4_gene/simulated/S4_gene_coexp_simulate.R"
  "Data/S5_gene/simulated/S5_ctcomp_coexp_simulate.R"
  "Data/S6_gene/simulated/S6_gene_ctcomp_coexp_simulate.R"
  "Data/S7_gene/simulated/S7_gene_STCOS_simulation.R"
)

for script in "${scripts[@]}"; do
  for replicate in 1 2 3 4 5; do
    printf 'Running %s, replicate %s\n' "${script}" "${replicate}"
    Rscript "${script}" "${replicate}"
  done
done

printf 'Completed S1-S7 with five replicates per scenario.\n'
