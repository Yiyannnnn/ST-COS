#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 REPOSITORY STCOS_OUTPUT_ROOT CHECK_ROOT" >&2
  exit 2
fi

# The second argument is the STCOS_OUTPUT_ROOT used by the S1-S7 generators,
# which write <STCOS_OUTPUT_ROOT>/<scenario>_gene/simulated.

repository=$(cd "$1" && pwd)
regenerated_root=$(cd "$2" && pwd)
check_root=$3
resource_root=${STCOS_RESOURCE_ROOT:-"$(cd "$repository/.." && pwd)"}

required_resources=(
  Deconvolution/POLARIS_repo_20260223
  Deconvolution/stdgcn
)
for resource in "${required_resources[@]}"; do
  if [[ ! -d "$resource_root/$resource" ]]; then
    echo "Missing external method resource: $resource_root/$resource" >&2
    echo "Set STCOS_RESOURCE_ROOT to the directory containing the method resources." >&2
    exit 1
  fi
done

for scenario in S1 S2 S3 S4 S5 S6 S7; do
  source_dir="$regenerated_root/${scenario}_gene/simulated"
  if [[ ! -d "$source_dir" ]]; then
    echo "Missing regenerated inputs: $source_dir" >&2
    echo "Pass the STCOS_OUTPUT_ROOT used when running the S1-S7 generators." >&2
    exit 1
  fi
  # Each single-cell reference must have cells in rows: Cell_type followed by
  # exactly the gene columns of the matching spatial expression table.
  shopt -s nullglob
  expression_files=("$source_dir/${scenario}_gene_expression_"*.csv)
  shopt -u nullglob
  if (( ${#expression_files[@]} == 0 )); then
    echo "No expression tables in: $source_dir" >&2
    exit 1
  fi
  for expression in "${expression_files[@]}"; do
    replicate_csv=${expression##*_}
    reference="$source_dir/${scenario}_scref_${replicate_csv}"
    if [[ ! -f "$reference" ]]; then
      echo "Missing single-cell reference: $reference" >&2
      exit 1
    fi
    if [[ "Cell_type,$(head -n 1 "$expression" | tr -d '\r"')" != \
          "$(head -n 1 "$reference" | tr -d '\r"')" ]]; then
      echo "Reference columns do not match the expression genes: $reference" >&2
      echo "Expected Cell_type followed by the columns of: $expression" >&2
      exit 1
    fi
  done
done

mkdir -p "$check_root"/{Data,Deconvolution,Segmentation,SVG_detection,CT_SVG,coexpression,logs,metadata}
mkdir -p "$check_root/metadata/runtime_seeds"
mkdir -p "$check_root/coexpression"/{CSSN,hdWGCNA_Agent,smoothie_Agent,spMOCA_Agent,stIHC_Agent,binSpect_Agent}

for scenario in S1 S2 S3 S4 S5 S6 S7; do
  scenario_dir="$check_root/Data/${scenario}_gene"
  mkdir -p "$scenario_dir"
  ln -sfn "$regenerated_root/${scenario}_gene/simulated" "$scenario_dir/simulated"
done

# These links expose vendored method code without mixing check outputs with the
# manuscript-era result tree.
ln -sfn "$resource_root/Deconvolution/POLARIS_repo_20260223" "$check_root/Deconvolution/POLARIS_repo_20260223"
ln -sfn "$resource_root/Deconvolution/stdgcn" "$check_root/Deconvolution/stdgcn"
cp "$repository/reproducibility/METHOD_SEEDS.tsv" \
  "$check_root/metadata/method_seeds.tsv"
cp "$repository/reproducibility/RUNNER_DEPENDENCIES.tsv" \
  "$check_root/metadata/runner_dependencies.tsv"
cp "$repository/reproducibility/SOFTWARE_VERSIONS.tsv" \
  "$check_root/metadata/software_versions.tsv"

if git -C "$repository" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git -C "$repository" rev-parse HEAD > "$check_root/metadata/repository_head.txt"
  git -C "$repository" status --short > "$check_root/metadata/repository_status.txt"
  git -C "$repository" diff --binary > "$check_root/metadata/uncommitted.patch"
else
  cp "$repository/PACKAGE_GIT_COMMIT.txt" \
    "$check_root/metadata/repository_head.txt"
  cp "$repository/PACKAGE_GIT_STATUS.txt" \
    "$check_root/metadata/repository_status.txt"
  : > "$check_root/metadata/uncommitted.patch"
fi
(
  cd "$regenerated_root"
  find -H S[1-7]_gene/simulated -maxdepth 1 -type f -name '*.csv' -print0 |
    sort -z | xargs -0 sha256sum
) > "$check_root/metadata/regenerated_input_sha256.txt"
(
  cd "$repository"
  find reproducibility/analysis -type f \( -name '*.R' -o -name '*.py' \) \
    -print0 | sort -z | xargs -0 sha256sum
) > "$check_root/metadata/analysis_code_sha256.txt"

cat > "$check_root/metadata/run_scope.txt" <<EOF
Purpose: separate downstream validation on package-regenerated S1-S7 data
Repository: $repository
External method resource root: $resource_root
Generator output root (STCOS_OUTPUT_ROOT): $regenerated_root
Check output root: $check_root
Manuscript-era Data directory modified: no
Later Reviewer3 rerun modified: no
EOF

echo "Prepared: $check_root"
