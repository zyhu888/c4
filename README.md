# C⁴: Covariate Connectivity Combined Clustering
R package for Covariate Connectivity Combined Clustering (C⁴), integrating covariate similarity and network structure for community detection.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

C⁴ (**Covariate Connectivity Combined Clustering**) is an adaptive spectral clustering algorithm that fuses network connectivity and node covariates into a unified matrix representation for community detection.  

This R package provides functions to:
- Perform spectral clustering on adjacency matrix (Ng, Jordan, & Weiss, 2002). 
- Implement proposed C⁴, an adaptive method that fuses connectivity and covariates for data-driven community detection.
- Perform CASC（covariate-assisted spectral clustering, Binkiewicz, Vogelstein, & Rohe, 2017).

---

## Installation

You can install the development version of **C4** from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install this package from GitHub
devtools::install_github("zyhu888/c4")

# Load the package
library(c4)
