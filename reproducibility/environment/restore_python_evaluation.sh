#!/usr/bin/env bash
set -euo pipefail

# Do not let packages in ~/.local satisfy or perturb the frozen environment.
export PYTHONNOUSERSITE=1

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${1:-python3}
requirements=${script_dir}/python-evaluation-requirements.txt
pytorch_index=https://download.pytorch.org/whl/cpu
pyg_wheels=https://data.pyg.org/whl/torch-1.13.0+cpu.html

if ! command -v "${python_bin}" >/dev/null 2>&1; then
  printf 'Python executable not found: %s\n' "${python_bin}" >&2
  exit 2
fi
if ! command -v R >/dev/null 2>&1; then
  printf 'R must be available on PATH because the environment includes rpy2.\n' >&2
  exit 2
fi

python_version=$("${python_bin}" -c \
  'import sys; print(".".join(map(str, sys.version_info[:2])))')
if [[ "${python_version}" != "3.8" ]]; then
  printf 'Python 3.8 is required; found %s.\n' "${python_version}" >&2
  exit 2
fi

"${python_bin}" -m pip install \
  'pip==24.2' 'setuptools==75.1.0' 'wheel==0.44.0'

# Install Torch first so the exact source-built PyG extensions can inspect it.
"${python_bin}" -m pip install \
  --extra-index-url "${pytorch_index}" \
  'torch==1.13.1+cpu' 'torchvision==0.14.1+cpu'

"${python_bin}" -m pip install --no-build-isolation \
  --find-links "${pyg_wheels}" \
  'torch_cluster==1.6.3' \
  'torch_scatter==2.1.2' \
  'torch_sparse==0.6.18' \
  'torch-spline-conv==1.2.2+pt113cpu'

# The requirements file carries the same wheel indexes and the immutable
# STAGATE-pyG source revision. Already installed Torch/PyG pins are reused.
"${python_bin}" -m pip install --no-build-isolation -r "${requirements}"

"${python_bin}" -m pip check
