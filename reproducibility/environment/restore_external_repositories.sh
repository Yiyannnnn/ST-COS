#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
destination=${1:?Usage: $0 DESTINATION_ROOT}

restore_repository() {
  local name=$1
  local url=$2
  local commit=$3
  local patch=$4
  local patch_sha=$5
  local target="$destination/Deconvolution/$name"

  if [[ -e "$target" ]]; then
    printf 'Target already exists: %s\n' "$target" >&2
    exit 2
  fi

  printf '%s  %s\n' "$patch_sha" "$patch" | sha256sum -c -
  git clone "$url" "$target"
  git -C "$target" checkout --detach "$commit"
  git -C "$target" apply "$patch"
  git -C "$target" diff --check
}

mkdir -p "$destination/Deconvolution"

# Exact upstream commits and compatibility patches used by the runners.
restore_repository \
  POLARIS_repo_20260223 \
  https://github.com/JiawenChenn/POLARIS.git \
  fcb2ac440c5f540e42a9a6183c85c78061103fe6 \
  "$script_dir/external_patches/POLARIS_fcb2ac440c5f.patch" \
  4e00bf8c5c8cd1d753bc2bd025d376182b429e54184b47f95478f6d18cc8a93b

restore_repository \
  stdgcn \
  https://github.com/luoyuanlab/stdgcn.git \
  295b08c740bf9979098e59add5a0e511afa2cfb2 \
  "$script_dir/external_patches/STdGCN_295b08c740bf.patch" \
  bac3b7a36bc87b85bd349d47fe77fad6e723ec378500de07e2c14ac66f02442d
