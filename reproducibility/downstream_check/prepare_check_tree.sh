#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 REPOSITORY REGENERATED_ROOT CHECK_ROOT" >&2
  exit 2
fi

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

mkdir -p "$check_root"/{Data,Deconvolution,Segmentation,SVG_detection,CT_SVG,coexpression,logs,metadata}
mkdir -p "$check_root/metadata/runtime_seeds"
mkdir -p "$check_root/coexpression"/{CSSN,hdWGCNA_Agent,smoothie_Agent,spMOCA_Agent,stIHC_Agent,binSpect_Agent}

for scenario in S1 S2 S3 S4 S5 S6 S7; do
  scenario_dir="$check_root/Data/${scenario}_gene"
  source_dir="$regenerated_root/data/${scenario}_gene/simulated"
  if [[ ! -d "$source_dir" ]]; then
    echo "Missing regenerated inputs: $source_dir" >&2
    exit 1
  fi
  mkdir -p "$scenario_dir"
  ln -sfn "$source_dir" "$scenario_dir/simulated"
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
  find data -type f -name '*.csv' -print0 | sort -z | xargs -0 sha256sum
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
Regenerated input root: $regenerated_root
Check output root: $check_root
Manuscript-era Data directory modified: no
Later Reviewer3 rerun modified: no
EOF

echo "Prepared: $check_root"
