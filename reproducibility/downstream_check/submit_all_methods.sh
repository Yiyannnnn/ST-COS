#!/usr/bin/env bash
set -euo pipefail

repository=$(cd "$(dirname "$0")/../.." && pwd)
regenerated_root=${1:-"$repository/reproducibility/results/package_regeneration"}
check_root=${2:-$regenerated_root/downstream_checks}
runner="$repository/reproducibility/downstream_check/run_method_task.sh"
log_dir="$check_root/logs"

bash "$repository/reproducibility/downstream_check/prepare_check_tree.sh" \
  "$repository" "$regenerated_root" "$check_root"

methods=(
  BASS BAYESSPACE MAPLE STAGATE GRAPHST SPAGCN
  CARD RCTD STDECONVOLVE TANGRAM POLARIS STDGCN
  HEARTSVG SPARK SPATIALDE_MSC SPATIALDE2
  CELINA CSIDE CTSV STANCE SPVC
  SPMOCA CSSN STIHC BINSPECT HDWGCNA SMOOTHIE SPACEX
)

job_ids=()
job_table="$check_root/metadata/submitted_jobs.tsv"
printf 'method\tjob_id\ttime\tmemory\tconcurrency\n' > "$job_table"

for method in "${methods[@]}"; do
  memory=32G
  walltime=2-00:00:00
  concurrency=10
  cpus=1
  case "$method" in
    POLARIS|STDGCN|SPMOCA|HDWGCNA) memory=64G; walltime=4-00:00:00; concurrency=5 ;;
    BASS|BAYESSPACE|MAPLE|STAGATE|GRAPHST|SPAGCN) memory=48G; walltime=3-00:00:00; concurrency=8 ;;
    CELINA|CSIDE|CTSV|STANCE|SPVC|CSSN|STIHC) memory=48G; walltime=3-00:00:00; concurrency=6 ;;
  esac
  job=$(sbatch --parsable -p general -N 1 -n 1 --cpus-per-task="$cpus" \
    --mem="$memory" -t "$walltime" --array="0-34%${concurrency}" \
    -J "PKG_${method}" -o "$log_dir/${method}_%A_%a.out" -e "$log_dir/${method}_%A_%a.err" \
    "$runner" "$method" "$repository" "$check_root")
  job_ids+=("$job")
  printf '%s\t%s\t%s\t%s\t%s\n' "$method" "$job" "$walltime" "$memory" "$concurrency" >> "$job_table"
done

printf 'check_root=%s\n' "$check_root" > "$check_root/metadata/latest_submission.txt"
cat "$job_table"
echo "Check root: $check_root"
