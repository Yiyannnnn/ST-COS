# ST-COS reproducibility

This directory contains the minimal code and environment records needed to
regenerate S1-S7, rerun the included downstream methods, and run the four
reference-based simulators used in the comparison. Scenario scripts call the
exported `STCOS` package functions directly; they do not replace package
functions with a second simulator implementation.

## Layout

| Path | Contents |
| --- | --- |
| `analysis/Data/` | S1-S7 generators for Figures 2-3; one seeded replicate per call |
| `analysis/{Deconvolution,Segmentation,SVG_detection,CT_SVG,coexpression}/` | Downstream method runners and task-level evaluation scripts |
| `analysis/Reference_based/` | Standalone ST-COS, SRTsim, scDesign3, and scCube reference-based runners |
| `analysis/evaluation_rules.R` | Shared truth and metric definitions |
| `downstream_check/` | Slurm launchers for every included downstream method and S1-S7 replicate |
| `environment/` | R lock files, pinned Python requirements, external-repository patches, and restoration scripts |
| `data_external/` | Local location for prepared public reference datasets; data are not redistributed |
| `METHOD_SEEDS.tsv` | Method-specific base seeds |
| `RUNNER_DEPENDENCIES.tsv` | Direct package versions and external source commits used by each runner |
| `SOFTWARE_VERSIONS.tsv` | Complete software inventory for the recorded environments |
| `EXTERNAL_INPUTS.txt` | Public dataset identifiers, URLs, expected filenames, and input schemas |

## Restore the package environment

From the package root, use R 4.4.0:

```r
install.packages("renv")
renv::restore(lockfile = "reproducibility/environment/renv.lock")
```

Then install the included package source:

```bash
R CMD INSTALL .
```

The SpaceX and Python environments are separate. Their exact restoration
commands are in `environment/README.md`.

## Generate S1-S7

Each scenario script accepts a replicate number from 1 to 5. The script derives
and records its simulation seed from the scenario and replicate number.

```bash
export STCOS_OUTPUT_ROOT=/path/to/generated_scenarios
Rscript reproducibility/analysis/Data/S1_gene/simulated/S1_gene_simulate.R 1
```

Use the analogous script under `S2_gene` through `S7_gene` for the remaining
scenarios and repeat with replicate numbers 1-5.

## Run downstream methods

On a Slurm system, submit all included methods after generating S1-S7:

```bash
bash reproducibility/downstream_check/submit_all_methods.sh \
  /path/to/generated_scenarios /path/to/downstream_checks
```

`METHOD_SEEDS.tsv` defines the base seed for each method. The launcher combines
that base seed with the scenario and replicate indices, and records the value in
each job log. `RUNNER_DEPENDENCIES.tsv` maps every method to its runner and
software source.

## Run reference-based simulations

The shared prepared-input format is a gene-by-observation count CSV whose first
column contains gene identifiers, plus a coordinate CSV whose first column
contains matching observation identifiers and whose numeric columns are named
`row` and `col`. Run, for example:

```bash
Rscript reproducibility/analysis/Reference_based/stcos.R \
  counts.csv coords.csv output/stcos 20260922
```

The corresponding runners are `srtsim.R`, `scdesign3.R`, and `sccube.py`.
scCube additionally requires a pretrained model and its ordered training-gene
list. `EXTERNAL_INPUTS.txt` declares the public DLPFC Visium, hypothalamus
MERFISH, and olfactory-bulb Stereo-seq resources, expected filenames, and
dataset-specific prepared formats.

## Scope

The repository supplies executable generation, method, and task-evaluation
code. Large public reference datasets, pretrained models, and generated result
directories are inputs or outputs and are not redistributed in the package.
`SOFTWARE_VERSIONS.tsv` is a provenance inventory; the files under
`environment/` are the installable environment specifications.
