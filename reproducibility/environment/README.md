# Restorable software environments

`SOFTWARE_VERSIONS.tsv` is the complete provenance inventory. The files in
this directory are the corresponding installation specifications for the
separate environments used by the reproducibility workflows.

## ST-COS and S1-S7 (R 4.4.0)

From the package root, install `renv`, restore the package and scenario
dependencies, and then install ST-COS from the included source:

```r
install.packages("renv")
renv::restore(lockfile = "reproducibility/environment/renv.lock")
```

```bash
R CMD INSTALL .
```

This lock covers ST-COS, its tests, and the S1-S7 generation scripts. System
libraries required by packages such as `sf` must be installed by the host
operating system.

## SpaceX (R 4.1.0)

Restore the separate legacy environment with R 4.1.0:

```r
install.packages("renv")
renv::restore(lockfile = "reproducibility/environment/renv-spacex.lock")
```

The SpaceX lock retains its GitHub source record. It is separate because that
method was run under R 4.1.0 and should not be mixed with the primary R 4.4.0
library.

## Python environments

The evaluation workflow used Python 3.8.20, while the POLARIS workflow used
Python 3.12.0. Create an empty Python 3.8 environment, make R 4.4 available on
`PATH`, and restore the evaluation environment with one command:

```bash
bash reproducibility/environment/restore_python_evaluation.sh /path/to/python3.8
```

Restore the separate POLARIS environment with Python 3.12:

```bash
python -m pip install -r reproducibility/environment/python-polaris-requirements.txt
```

The evaluation requirements file embeds both the official PyTorch CPU wheel
index and the PyTorch Geometric wheel page for PyTorch 1.13 CPU builds. It also
pins STAGATE-pyG to the exact Git commit used in the analysis because that
package is not available from PyPI. The restoration script installs Torch
first, then builds the exact `torch_cluster`, `torch_scatter`, and
`torch_sparse` versions used in the analysis against that Torch installation,
and finally installs the remaining pins. It disables the Python user site so
packages in `~/.local` cannot silently satisfy or conflict with the recorded
environment. No additional wheel-index arguments or manual repository checkout
is required.

Use Python 3.8.20 for `python-evaluation-requirements.txt` and Python 3.12.0 for
`python-polaris-requirements.txt`. External method repositories and pretrained
weights are analysis inputs and are not vendored into these software
specifications.

## External method repositories

POLARIS and STdGCN use exact upstream commits plus the small compatibility
patches pinned in `restore_external_repositories.sh`. Restore both beneath a
resource root with one command:

```bash
bash reproducibility/environment/restore_external_repositories.sh \
  /path/to/STCOS_RESOURCE_ROOT
```

The script verifies each patch checksum before cloning, checks out the pinned
commit in detached-HEAD state, applies the recorded patch, and runs
`git diff --check`. STAGATE-pyG is restored from its pinned Git commit by the
Python evaluation requirements and therefore does not require a second source
checkout.

## Slurm launcher settings

`downstream_check/run_method_task.sh` locates each environment through these
variables:

| Variable | Default | Points to |
| --- | --- | --- |
| `STCOS_USER_R_LIBRARY` | `R_LIBS_USER` | R 4.4 library with the downstream method packages |
| `STCOS_RESOURCE_ROOT` | parent of the repository | Directory holding `Deconvolution/POLARIS_repo_20260223` and `Deconvolution/stdgcn` from `restore_external_repositories.sh` |
| `STCOS_PYTHON_EVAL` | `$STCOS_RESOURCE_ROOT/env/STCOS_eval/bin/python` | Python 3.8 evaluation environment |
| `STCOS_PYTHON_POLARIS` | `$STCOS_RESOURCE_ROOT/env/AD_drug/bin/python` | Python 3.12 POLARIS environment |
| `STCOS_R_SPACEX` | `$STCOS_RESOURCE_ROOT/env/r410/bin/Rscript` | R 4.1.0 `Rscript` with the SpaceX library |
| `STCOS_STDGCN_ROOT` | `$STCOS_RESOURCE_ROOT/Deconvolution/stdgcn` | Patched STdGCN checkout |

Before submitting, `submit_all_methods.sh` checks that `STCOS_USER_R_LIBRARY`
contains every package `RUNNER_DEPENDENCIES.tsv` lists for the R 4.4 methods.
It stops if any are missing, so jobs never fail later on a missing package.
