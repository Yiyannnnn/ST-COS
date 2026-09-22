# ST-COS

ST-COS: Simulator of Spatial Transcriptomic Data Preserving Gene-Gene Co-expression Patterns

ST-COS combines spatial mean patterns, cell-type composition fields,
negative-binomial marginals, and editable latent Gaussian gene-gene
correlation matrices. The specified matrix controls latent pairwise
correlation; realized Pearson and Spearman correlations among counts can be
attenuated by discretization, sparsity, and spot aggregation.

For spot-based data, `generate_spot_count()` constructs the continuous
abundance-weighted mixture described in the manuscript. It does not simulate
independent individual cells within a spot. Reference-based estimation returns
one global spot-level dependence matrix; cell-type-specific matrices must be
provided externally or specified in the reference-free workflow.

<img width="966" height="818" alt="image" src="https://github.com/user-attachments/assets/a2eeff20-c4cf-4f8f-ba6f-e94c725fc86e" />



## Installation

You can install the development version of STCOS from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Yiyannnnn/ST-COS")
```

## Example

Please refer to https://yiyannnnn.github.io/STCOS_tutorial/Main.html


## License

This project is licensed under the MIT License - see the LICENSE file for details
