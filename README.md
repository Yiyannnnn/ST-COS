# STCOS

ST-COS: Simulator of Spatial Transcriptomic Data Preserving Gene-Gene Co-expression Patterns

## Installation

You can install the development version of STCOS from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Yiyannnnn/STCOS")
```

## Example

Basic usage:

```r
library(STCOS)

# Generate coordinates
coords <- generate_coordinates(x_len = 100, 
                             y_len = 100, 
                             pattern = "hex",
                             spot_distance = 10)

# Generate patterns
# Add your example code here
```

## License

This project is licensed under the MIT License - see the LICENSE file for details
