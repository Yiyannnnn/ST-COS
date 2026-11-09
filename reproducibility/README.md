# ST-COS reproducibility code

This directory contains the manuscript analysis and figure-generation code,
the audited historical simulation helper, and the maintained package workflow.
These sources are kept separate so that the reported manuscript results are
not conflated with corrections made to the public package after review. See
[`ANALYSIS_PROVENANCE.md`](ANALYSIS_PROVENANCE.md) for the audit of the
historical code and outputs.

## Start here

From the package root:

```bash
export R_LIBS_USER="${HOME}/Rlibs"
R CMD INSTALL .
cd reproducibility
Rscript validate_core.R
Rscript audit_legacy_s2.R
Rscript check_environment.R
```

`validate_core.R` directly checks the issues raised during review: seed
propagation across genes, positive-definite repair by the declared uniform
eigenvalue shift, the Cholesky orientation, and binomial thinning. It writes its
results to `results/core_validation/`.

`audit_legacy_s2.R` uses the helper preserved under `historical/` to reconstruct
the recoverable historical S2 path with seed 123, 5% binomial thinning, and
abundance weighting, then verifies the archived SHA-256 hash. This audit is
separate from the maintained package implementation.

`MANIFEST.tsv` maps each manuscript result to its generating and plotting
scripts. `SEEDS.tsv` records the scenario seeds, and `SOFTWARE.tsv` records the
validated software environment. The complete analysis code is under
`analysis/`.

## External data

Public reference datasets are not duplicated in this repository. Place them
under `data_external/` using the layout in
[`data_external/README.md`](data_external/README.md). The simulated scenario
scripts do not require these data.

## Main workflows

### Simulated scenarios S1-S7

The scenario scripts are supplied to expose the complete simulation and
downstream workflow. From `reproducibility/`, users can generate a new set of
five seeded replicates with:

```bash
bash run_scenarios.sh
```

Each generator can also be run separately from `reproducibility/analysis`
with a replicate number from 1 to 5. For example:

```bash
Rscript Data/S1_gene/simulated/S1_gene_simulate.R 1
Rscript Data/S7_gene/simulated/S7_gene_STCOS_simulation.R 1
```

The corresponding seeds are listed in `../SEEDS.tsv`. Downstream method
scripts are grouped by task in `Deconvolution/`, `Segmentation/`,
`SVG_detection/`, `CT_SVG/`, and `coexpression/`. Summary plotting scripts
write both per-replicate tables and publication figures.

### DLPFC reference-based analysis

For each of the 12 sections listed in `../SEEDS.tsv`:

```bash
Rscript Reference_based/DLPFC_run_methods_by_slide.R 151673 all
python3 Reference_based/DLPFC_scCube_generate.py 151673
Rscript Reference_based/DLPFC_metrics_by_slide.R 151673
Rscript Reference_based/DLPFC_gene_corr_by_slide.R 151673
```

After all sections finish:

```bash
Rscript Reference_based/DLPFC_plot_4summary.R
Rscript Reference_based/DLPFC_figure4_statistics.R
Rscript Reference_based/DLPFC_strong_pair_recovery.R
```

### MERFISH cross-platform analysis

For each field `04`, `09`, `14`, `19`, and `24`:

```bash
Rscript Reference_based/merfish_run_methods_by_slide.R 04
```

Then run the MERFISH comparison scripts listed in `MANIFEST.tsv`.

### Revision experiments

```bash
# Latent versus realized count-scale correlation
Rscript Independent_simulator_anchor/simulate_count_correlation_attenuation.R

# Ectopic white-matter SVG counterfactual
Rscript Independent_simulator_anchor/DLPFC_ectopic_WM_SVG_demo.R

# Parameter sensitivity for one section
Rscript Reference_based/DLPFC_parameter_sensitivity.R all 151673

# Combine the multi-section sensitivity results
Rscript Reference_based/combine_DLPFC_multislide_sensitivity.R
```

## Compute-intensive methods

Several downstream packages require separate environments or scheduler jobs.
The exact launch scripts are included, notably
`Independent_simulator_anchor/submit_S7SRT_heavy_downstream.sh`,
`Reference_based/submit_DLPFC_jobs.sh`, and the method-specific files in each
task directory. Run these only after `check_environment.R` reports the needed
packages and external programs.

## Model and code provenance

The package API in the repository root is the authoritative maintained
implementation for new simulations.
In particular:

- `generate_gene_pattern()` sets the seed once per call; deterministic helper
  functions do not reset it inside the gene loop.
- `makespd()` applies the minimum uniform shift needed to raise the smallest
  eigenvalue to `epsilon`, then renormalizes the repaired matrix to unit
  diagonal.
- `generate_batch_count()` performs binomial thinning, not zero-inflated
  mixture sampling.
- `generate_spot_count()` implements the abundance-weighted spot-level model
  described in the revised manuscript.

The S1-S7 generators under `analysis/Data/` retain the abundance-weighted
spot-level formulation described in the revised manuscript. Their seeds in
`SEEDS.tsv` are provided for deterministic new runs and are not presented as
recovered seeds for every historical output. The manuscript analyses remain
linked to their archived inputs and the historical helper documented in
`ANALYSIS_PROVENANCE.md`; correcting the maintained package does not replace or
silently regenerate those reported results.

## Reproducibility record

From `reproducibility/`, run:

```bash
Rscript check_environment.R
```

This writes R and package versions to `results/session/`. Record the Git commit
alongside archived result files. Do not commit raw controlled data, package
caches, or third-party model weights.
