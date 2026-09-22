# ST-COS

**ST-COS** is an R package for simulating spatial transcriptomics data with
configurable spatial expression patterns, cell-type compositions,
negative-binomial marginals, and editable gene-gene dependence.

ST-COS supports two workflows:

- **Reference-free simulation:** design coordinates, cell-type composition,
  spatial gene patterns, and cell-type-specific latent correlation matrices.
- **Reference-based simulation:** estimate marginal, spatial, and global
  spot-level dependence parameters from an observed spatial transcriptomics
  dataset and generate new profiles from the fitted model.
- **Image-guided design:** upload a histology image or coordinate file and
  configure a reference-free simulation interactively in the hosted Shiny
  application.

<img width="966" alt="ST-COS simulation and downstream analysis workflow" src="man/figures/figure1.png" />

## Installation

Install the development version from GitHub:

```r
install.packages("devtools") # Skip if devtools is already installed
devtools::install_github("Yiyannnnn/ST-COS")
```

Then load the package:

```r
library(STCOS)
```

## Image-guided application

Open the hosted ST-COS application from an interactive R session:

```r
launch_stcos_app()
```

The application supports histology-image upload and segmentation, uploaded
coordinates, interactive region selection, gene-pattern assignment,
cell-type composition, and gene-dependence controls. On a headless system,
retrieve the URL without trying to open a browser:

```r
launch_stcos_app(browser = FALSE)
```

The hosted application is also available directly at
[yiyanz.shinyapps.io/st-cos](https://yiyanz.shinyapps.io/st-cos/).

## Quick start

The example below creates two cell types, four genes, cell-type-specific
spatial patterns, and one correlated gene pair per cell type.

```r
library(STCOS)

set.seed(1)

# 1. Create spatial locations.
coord <- generate_coordinates(
  x_len = 10,
  y_len = 10,
  pattern = "grid",
  spot_distance = 1
)
rownames(coord) <- paste0("spot_", seq_len(nrow(coord)))

# 2. Assign a spatial mean pattern to each gene and cell type.
genes <- paste0("gene_", 1:4)
gene_pattern <- data.frame(
  Excitatory = c("hotspot", "hotspot", "no_pattern", "no_pattern"),
  Inhibitory = c("no_pattern", "no_pattern", "streak", "streak"),
  row.names = genes
)
mean_param <- generate_gene_pattern(coord, gene_pattern, seed = 1)

# 3. Define the genes participating in each dependence block.
gene_relationship <- data.frame(
  Excitatory = c(1, 1, NA, NA),
  Inhibitory = c(NA, NA, 1, 1),
  row.names = genes
)

# 4. Specify cell-type-specific latent correlation matrices.
rho <- array(
  0,
  dim = c(4, 4, 2),
  dimnames = list(genes, genes, c("Excitatory", "Inhibitory"))
)
diag(rho[, , "Excitatory"]) <- 1
diag(rho[, , "Inhibitory"]) <- 1
rho[1, 2, "Excitatory"] <- rho[2, 1, "Excitatory"] <- 0.6
rho[3, 4, "Inhibitory"] <- rho[4, 3, "Inhibitory"] <- 0.6

# 5. Generate one spot-level expression profile for each cell type.
cell_type_profiles <- generate_true_count(
  coord = coord,
  mean_param = mean_param,
  gene_relationship = gene_relationship,
  cov_gene = rho,
  theta_g = list(
    Excitatory = rep(10, 4),
    Inhibitory = rep(10, 4)
  )
)

# 6. Mix the cell-type profiles using spatially varying abundances.
p_excitatory <- plogis((coord$x_coord - 5) / 1.5)
celltype_proportion <- cbind(
  Excitatory = p_excitatory,
  Inhibitory = 1 - p_excitatory
)
rownames(celltype_proportion) <- rownames(coord)

spot_expression <- generate_spot_count(
  true_count = cell_type_profiles,
  nGenes = length(genes),
  cell_number = 8,
  celltype_proportion = celltype_proportion
)

dim(spot_expression)
spot_expression[1:4, ]
```

## Reference-based simulation

For an observed count matrix with genes in rows and spots in columns, estimate
the reference-based parameters with `get_real_param()`:

```r
# expr: genes x spots count matrix
# spatial_coord: data frame with columns named "row" and "col"
fitted <- get_real_param(
  expr = expr,
  coord = spatial_coord,
  gene_names = rownames(expr),
  n_cores = 4
)

simulated <- generate_true_count(
  coord = spatial_coord,
  mean_param = fitted$mean_param,
  gene_relationship = fitted$gene_relationship,
  cov_gene = fitted$cov_gene,
  theta_g = fitted$theta_g
)
```

The reference-based workflow estimates one global spot-level dependence
matrix. Cell-type-specific dependence matrices must be supplied externally or
specified through the reference-free workflow.

## Input conventions

| Object | Expected structure |
| --- | --- |
| Reference counts | Genes in rows and spots in columns |
| Reference coordinates | One row per spot, with columns `row` and `col` |
| Reference-free coordinates | One row per spot, with columns `x_coord` and `y_coord` |
| Mean parameters | Named list of spot-by-gene matrices, one per cell type |
| Cell-type proportions | Spot-by-cell-type matrix; every row must sum to one |
| Dependence parameters | Gene-by-gene matrix or array containing one matrix per cell type |

## Model interpretation

- The supplied dependence matrix controls **latent Gaussian pairwise
  correlation**. Realized Pearson and Spearman correlations among simulated
  counts may be attenuated by negative-binomial discretization, sparsity, and
  spot-level mixing.
- `generate_spot_count()` constructs a continuous abundance-weighted mixture
  of cell-type-specific spot profiles. It does not simulate or sum independent
  individual cells within a spot.
- `eigen_floor = 1e-6` is the default numerical floor used when repairing a
  latent correlation matrix before Cholesky factorization.

## Documentation and support

- [Full tutorial](https://yiyannnnn.github.io/STCOS_tutorial/)
- [Tutorial source](https://github.com/Yiyannnnn/STCOS_tutorial)
- [Report a problem or request a feature](https://github.com/Yiyannnnn/ST-COS/issues)
- In R, run `help(package = "STCOS")` or `?generate_true_count` for function
  documentation.

## License

ST-COS is released under the [MIT License](LICENSE).
