c4: Covariate Connectivity Combined Clustering
================

[![License:
MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Covariate Connectivity Combined Clustering ($\text{C}^4$)

The package implements the Covariate Connectivity Combined Clustering
($\text{C}^4$) method, which integrates both network topology and
node-level covariates to improve community detection in weighted
networks.

This R package provides functions to:

- simulate weighted networks under stochastic block models with
  user-defined block structures and weight distributions;

- perform $\text{C}^4$ for joint network-covariate community detection
  (Hu et al. 2025);

- apply covariate-assisted spectral clustering as a competing method
  (Binkiewicz, Vogelstein, and Rohe 2017).

## Installation

You can install $\texttt{c4}$ from github with:

``` r
# install.packages("devtools")
devtools::install_github("zyhu888/c4")
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-binkiewicz2017covariate" class="csl-entry">

Binkiewicz, N., J. T. Vogelstein, and K. Rohe. 2017. “Covariate-Assisted
Spectral Clustering.” *Biometrika* 104 (2): 361–77.

</div>

<div id="ref-hu2025covariate" class="csl-entry">

Hu, Z., W. Li, J. Yan, and P. Zhang. 2025. “Covariate Connectivity
Combined Clustering for Weighted Networks.” University of Connecticut.

</div>

</div>
