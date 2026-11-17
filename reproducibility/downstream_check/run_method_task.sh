#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 METHOD REPOSITORY CHECK_ROOT" >&2
  exit 2
fi

method=$1
repository=$(cd "$2" && pwd)
check_root=$(cd "$3" && pwd)
resource_root=${STCOS_RESOURCE_ROOT:-"$(cd "$repository/.." && pwd)"}
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

if (( task_id < 0 || task_id > 34 )); then
  echo "SLURM_ARRAY_TASK_ID must be in 0..34" >&2
  exit 2
fi

scenario="S$((task_id / 5 + 1))"
replicate=$((task_id % 5 + 1))
suffix="_${replicate}"
analysis="$repository/reproducibility/analysis"
seed_table="$repository/reproducibility/METHOD_SEEDS.tsv"
scenario_number=${scenario#S}
base_seed=$(awk -F '\t' -v method="$method" \
  'NR > 1 && $1 == method { print $2; exit }' "$seed_table")
if [[ -z "$base_seed" ]]; then
  echo "No base seed recorded for method: $method" >&2
  exit 2
fi
method_seed=$((base_seed + 100 * scenario_number + replicate))
export STCOS_METHOD_SEED=$method_seed
export PYTHONHASHSEED=$method_seed

module load r/4.4.0
user_r_library=${STCOS_USER_R_LIBRARY:-${R_LIBS_USER:-}}
export R_LIBS_USER="$check_root/package_library${user_r_library:+:$user_r_library}"
export MPLCONFIGDIR="/tmp/stcos_check_mpl_${SLURM_JOB_ID:-$$}_${task_id}"
export XDG_CACHE_HOME="/tmp/stcos_check_cache_${SLURM_JOB_ID:-$$}_${task_id}"
export POT_BACKEND_DISABLE_TENSORFLOW=1
export R_PROFILE_USER="$repository/reproducibility/downstream_check/method_seed.Rprofile"
export PYTHONPATH="$repository/reproducibility/downstream_check${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

r_run() {
  Rscript "$@"
}

py=${STCOS_PYTHON_EVAL:-"$resource_root/env/STCOS_eval/bin/python"}
if [[ ! -x "$py" ]]; then
  echo "Missing Python environment: $py" >&2
  exit 1
fi

cd "$check_root"
seed_record="$check_root/metadata/runtime_seeds/${method}_${scenario}_${replicate}.tsv"
printf 'method\tscenario\treplicate\tbase_seed\teffective_seed\n' > "$seed_record"
printf '%s\t%s\t%s\t%s\t%s\n' \
  "$method" "$scenario" "$replicate" "$base_seed" "$method_seed" >> "$seed_record"
echo "method=$method scenario=$scenario replicate=$replicate seed=$method_seed host=$(hostname) started=$(date -Is)"

case "$method" in
  BASS) r_run "$analysis/Segmentation/BASS.R" "$scenario" "$suffix" ;;
  BAYESSPACE) r_run "$analysis/Segmentation/BayesSpace.R" "$scenario" "$suffix" ;;
  MAPLE) r_run "$analysis/Segmentation/MAPLE.R" "$scenario" "$suffix" ;;
  STAGATE) "$py" "$analysis/Segmentation/STAGATE.py" "$scenario" "$suffix" --graph_model KNN --k_cutoff 8 ;;
  GRAPHST) "$py" "$analysis/Segmentation/run_GraphST.py" "$scenario" "$suffix" ;;
  SPAGCN) "$py" "$analysis/Segmentation/run_SpaGCN.py" "$scenario" "$suffix" ;;

  CARD) r_run "$analysis/Deconvolution/CARD.R" "$scenario" "$suffix" ;;
  RCTD) r_run "$analysis/Deconvolution/RCTD.R" "$scenario" "$suffix" ;;
  STDECONVOLVE) r_run "$analysis/Deconvolution/STdeconvolve.R" "$scenario" "$suffix" ;;
  TANGRAM) "$py" "$analysis/Deconvolution/Tangram.py" "$scenario" "$suffix" ;;
  POLARIS)
    export POLARIS_PYTHON=${STCOS_PYTHON_POLARIS:-"$resource_root/env/AD_drug/bin/python"}
    "$py" "$analysis/Deconvolution/POLARIS.py" "$scenario" "$suffix" --epochs 200 --batch_size 256
    ;;
  STDGCN)
    "$py" "$repository/reproducibility/downstream_check/prepare_stdgcn_inputs.py" \
      "$check_root" "$scenario" "$suffix"
    export PYTHONPATH="${STCOS_STDGCN_ROOT:-$resource_root/Deconvolution/stdgcn}:$repository/reproducibility/downstream_check"
    "$py" "$analysis/Deconvolution/run_STdGCN.py" "$scenario" "$suffix"
    ;;

  HEARTSVG) r_run "$analysis/SVG_detection/HEARTSVG.R" "$scenario" "$suffix" ;;
  SPARK) r_run "$analysis/SVG_detection/SPARK.R" "$scenario" "$suffix" ;;
  SPATIALDE_MSC) "$py" "$analysis/SVG_detection/SpatialDE_MSC.py" "$scenario" "$suffix" ;;
  SPATIALDE2) "$py" "$analysis/SVG_detection/run_SpatialDE.py" "$scenario" "$suffix" ;;

  CELINA) r_run "$analysis/CT_SVG/CELINA.R" "$scenario" "$suffix" both ;;
  CSIDE) r_run "$analysis/CT_SVG/CSIDE.R" "$scenario" "$suffix" ;;
  CTSV) r_run "$analysis/CT_SVG/CTSV.R" "$scenario" "$suffix" ;;
  STANCE) r_run "$analysis/CT_SVG/STANCE.R" "$scenario" "$suffix" ;;
  SPVC) r_run "$analysis/CT_SVG/spVC.R" "$scenario" "$suffix" ;;

  SPMOCA) r_run "$analysis/coexpression/spMOCA.R" "$scenario" "$suffix" ;;
  CSSN) r_run "$analysis/coexpression/CSSN.R" "$scenario" "$suffix" ;;
  STIHC) r_run "$analysis/coexpression/stIHC.R" "$scenario" "$suffix" ;;
  BINSPECT) r_run "$analysis/coexpression/binSpect.R" "$scenario" "$suffix" ;;
  HDWGCNA)
    hd_work="$check_root/metadata/hdwgcna_work/${scenario}_${replicate}"
    mkdir -p "$hd_work"
    ln -sfn "$check_root/Data" "$hd_work/Data"
    ln -sfn "$check_root/coexpression" "$hd_work/coexpression"
    (cd "$hd_work" && r_run "$analysis/coexpression/hdWGCNA.R" "$scenario" "$suffix")
    ;;
  SMOOTHIE) "$py" "$analysis/coexpression/smoothie.py" "$scenario" "$suffix" ;;
  SPACEX)
    spacex_r=${STCOS_R_SPACEX:-"$resource_root/env/r410/bin/Rscript"}
    env -u R_LIBS_USER R_PROFILE_USER="$R_PROFILE_USER" \
      STCOS_METHOD_SEED="$STCOS_METHOD_SEED" \
      "$spacex_r" --no-site-file --no-environ --no-restore --no-save \
      "$analysis/coexpression/SpaceX.R" "$scenario" "$suffix"
    ;;
  *) echo "Unknown method: $method" >&2; exit 2 ;;
esac

echo "method=$method scenario=$scenario replicate=$replicate seed=$method_seed finished=$(date -Is)"
