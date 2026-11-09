#!/usr/bin/env bash
set -euo pipefail

ROOT="."
WORKER="${ROOT}/Independent_simulator_anchor/run_S7SRT_heavy_method.sh"

PARTITION="${PARTITION:-general}"
TIME_LIMIT="${TIME_LIMIT:-7-}"
MEM_LIMIT="${MEM_LIMIT:-50g}"
ARRAY_SPEC="${ARRAY_SPEC:-1-5%5}"

timestamp="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${ROOT}/logs_S7SRT_heavy_downstream_${timestamp}"
mkdir -p "$LOG_DIR"

submitted="${LOG_DIR}/submitted_jobs.tsv"
printf "method\tjob_id\tarray_spec\tdependency\tlog_prefix\n" > "$submitted"

submit_array() {
  local method="$1"
  local name="$2"
  local mem="${3:-$MEM_LIMIT}"
  local time="${4:-$TIME_LIMIT}"
  local cpus="${5:-1}"
  local dependency="${6:-}"
  local dep_arg=()

  if [[ -n "$dependency" ]]; then
    dep_arg=(--dependency="afterok:${dependency}")
  fi

  local jid
  jid="$(sbatch --parsable \
    -p "$PARTITION" \
    -N 1 \
    --ntasks=1 \
    --cpus-per-task="$cpus" \
    --mem="$mem" \
    -t "$time" \
    --array="$ARRAY_SPEC" \
    --job-name="$name" \
    --output="${LOG_DIR}/${name}_%A_%a.out" \
    --error="${LOG_DIR}/${name}_%A_%a.err" \
    "${dep_arg[@]}" \
    --wrap="bash -lc 'bash ${WORKER} ${method} \${SLURM_ARRAY_TASK_ID}'")"

  printf "%s\t%s\t%s\t%s\t%s\n" "$method" "$jid" "$ARRAY_SPEC" "${dependency:-NA}" "${LOG_DIR}/${name}" >> "$submitted"
  echo "$jid"
}

submit_array SEG_BASS SEG_BASS_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array SEG_BayesSpace SEG_BayesSpace_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array SEG_MAPLE SEG_MAPLE_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array SEG_STAGATE SEG_STAGATE_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
submit_array SEG_GraphST SEG_GraphST_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
submit_array SEG_SpaGCN SEG_SpaGCN_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null

submit_array DEC_CARD DEC_CARD_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array DEC_RCTD DEC_RCTD_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array DEC_STdeconvolve DEC_STdeconvolve_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array DEC_Tangram DEC_Tangram_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
submit_array DEC_POLARIS DEC_POLARIS_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
stdgcn_prep_id="$(submit_array DEC_STdGCN_PREP DEC_STdGCN_PREP_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2)"
stdgcn_prep_id="${stdgcn_prep_id%%;*}"
submit_array DEC_STdGCN_RUN DEC_STdGCN_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 "$stdgcn_prep_id" >/dev/null

submit_array SVG_HEARTSVG SVG_HEARTSVG_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array SVG_SPARK SVG_SPARK_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array SVG_SpatialDE_MSC SVG_SpatialDE_MSC_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
submit_array SVG_SpatialDE2 SVG_SpatialDE2_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null

submit_array CTSVG_CELINA CTSVG_CELINA_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array CTSVG_CELINA_TRUEPROP CTSVG_CELINA_TRUEPROP_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array CTSVG_CSIDE CTSVG_CSIDE_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array CTSVG_CTSV CTSVG_CTSV_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array CTSVG_STANCE CTSVG_STANCE_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array CTSVG_SPVC CTSVG_SPVC_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null

submit_array COEXP_SPMOCA COEXP_SPMOCA_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array COEXP_SPACEX COEXP_SPACEX_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array COEXP_CSSN COEXP_CSSN_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array COEXP_HDWGCNA COEXP_HDWGCNA_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
submit_array COEXP_SMOOTHIE COEXP_SMOOTHIE_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 4 >/dev/null
submit_array COEXP_BINSPECT COEXP_BINSPECT_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null
submit_array COEXP_STIHC COEXP_STIHC_S7SRT "$MEM_LIMIT" "$TIME_LIMIT" 2 >/dev/null

echo "Submitted S7SRT heavy downstream arrays."
echo "Log directory: ${LOG_DIR}"
echo "Job table: ${submitted}"
