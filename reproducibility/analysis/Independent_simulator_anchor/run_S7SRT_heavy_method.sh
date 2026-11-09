#!/usr/bin/env bash
set -euo pipefail

method="${1:?method name is required}"
rep="${2:?replicate index is required}"

ROOT="."
S="S7SRT"
suffix="_${rep}"

cd "$ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MPLCONFIGDIR="/tmp/mpl_${USER:-yiyan}_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-${rep}}"
mkdir -p "$MPLCONFIGDIR"

run_r44() {
  module load r/4.4.0
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
  Rscript "$@"
}

run_py_stcos() {
  module load anaconda
  module load r/4.4.0
  conda activate STCOS_eval
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
  export POT_BACKEND_DISABLE_TENSORFLOW=1
  export PYTHONNOUSERSITE=1
  export PYTHONPATH=""
  python "$@"
}

run_stdgcn_py() {
  module load anaconda
  module load r/4.4.0
  conda activate STCOS_eval
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
  export POT_BACKEND_DISABLE_TENSORFLOW=1
  export PYTHONNOUSERSITE=1
  export PYTHONPATH="./Deconvolution/stdgcn"
  python "$@"
}

case "$method" in
  SEG_BASS)
    run_r44 Segmentation/BASS.R "$S" "$suffix"
    ;;
  SEG_BayesSpace)
    run_r44 Segmentation/BayesSpace.R "$S" "$suffix"
    ;;
  SEG_MAPLE)
    run_r44 Segmentation/MAPLE.R "$S" "$suffix"
    ;;
  SEG_STAGATE)
    run_py_stcos Segmentation/STAGATE.py "$S" "$suffix"
    ;;
  SEG_GraphST)
    run_py_stcos Segmentation/run_GraphST.py "$S" "$suffix"
    ;;
  SEG_SpaGCN)
    run_py_stcos Segmentation/run_SpaGCN.py "$S" "$suffix"
    ;;
  DEC_CARD)
    run_r44 Deconvolution/CARD.R "$S" "$suffix"
    ;;
  DEC_RCTD)
    run_r44 Deconvolution/RCTD.R "$S" "$suffix"
    ;;
  DEC_STdeconvolve)
    run_r44 Deconvolution/STdeconvolve.R "$S" "$suffix"
    ;;
  DEC_Tangram)
    run_py_stcos Deconvolution/Tangram.py "$S" "$suffix"
    ;;
  DEC_POLARIS)
    ./env/AD_drug/bin/python Deconvolution/POLARIS.py "$S" "$suffix" --epochs 200 --batch_size 512
    ;;
  DEC_STdGCN_PREP)
    module load r/4.4.0
    export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
    bash Data/STdGCN_pre_data.sh "$S" "$suffix"
    ;;
  DEC_STdGCN_RUN)
    run_stdgcn_py Deconvolution/run_STdGCN.py "$S" "$suffix"
    ;;
  SVG_HEARTSVG)
    run_r44 SVG_detection/HEARTSVG.R "$S" "$suffix"
    ;;
  SVG_SPARK)
    run_r44 SVG_detection/SPARK.R "$S" "$suffix"
    ;;
  SVG_SpatialDE_MSC)
    run_py_stcos SVG_detection/SpatialDE_MSC.py "$S" "$suffix"
    ;;
  SVG_SpatialDE2)
    run_py_stcos SVG_detection/run_SpatialDE.py "$S" "$suffix"
    ;;
  CTSVG_CELINA)
    run_r44 CT_SVG/CELINA.R "$S" "$suffix"
    ;;
  CTSVG_CELINA_TRUEPROP)
    run_r44 CT_SVG/CELINA_trueprop.R "$S" "$suffix"
    ;;
  CTSVG_CSIDE)
    run_r44 CT_SVG/CSIDE.R "$S" "$suffix"
    ;;
  CTSVG_CTSV)
    run_r44 CT_SVG/CTSV.R "$S" "$suffix"
    ;;
  CTSVG_STANCE)
    run_r44 CT_SVG/STANCE.R "$S" "$suffix"
    ;;
  CTSVG_SPVC)
    run_r44 CT_SVG/spVC.R "$S" "$suffix"
    ;;
  COEXP_SPMOCA)
    run_r44 coexpression/spMOCA.R "$S" "$suffix"
    ;;
  COEXP_SPACEX)
    ./env/r410/bin/Rscript coexpression/SpaceX.R "$S" "$suffix"
    ;;
  COEXP_CSSN)
    run_r44 coexpression/CSSN.R "$S" "$suffix"
    ;;
  COEXP_HDWGCNA)
    run_r44 coexpression/hdWGCNA.R "$S" "$suffix"
    ;;
  COEXP_SMOOTHIE)
    run_py_stcos coexpression/smoothie.py "$S" "$suffix"
    ;;
  COEXP_BINSPECT)
    run_r44 coexpression/binSpect.R "$S" "$suffix"
    ;;
  COEXP_STIHC)
    run_r44 coexpression/stIHC.R "$S" "$suffix"
    ;;
  *)
    echo "Unknown S7SRT heavy method: $method" >&2
    exit 2
    ;;
esac
