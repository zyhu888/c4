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
  (Hu et al. 2026);

- apply covariate-assisted spectral clustering as a competing method
  (Binkiewicz et al. 2017).

## Installation

You can install $\texttt{c4}$ from github with:

``` r
# install.packages("devtools")
devtools::install_github("zyhu888/c4")
```

## Citation

If you use `c4` in your research, please cite:

``` bibtex
@article{hu2026covariate,
  title   = {Covariate Connectivity Combined Clustering for Weighted Networks},
  author  = {Hu, Zeyu and Li, Wenrui and Yan, Jun and Zhang, Panpan},
  journal = {Physica A: Statistical Mechanics and its Applications},
  volume  = 704,
  pages   = 132177,
  year    = 2026
}
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-hu2026covariate" class="csl-entry">

Hu, Zeyu, Wenrui Li, Jun Yan, and Panpan Zhang. 2026. “[Covariate
Connectivity Combined Clustering for Weighted
Networks](https://doi.org/10.1016/j.physa.2026.132177).” *Physica A:
Statistical Mechanics and its Applications* 704: 132177.

</div>

<div id="ref-binkiewicz2017covariate" class="csl-entry">

Binkiewicz, N., J. T. Vogelstein, and K. Rohe. 2017. “Covariate-Assisted
Spectral Clustering.” *Biometrika* 104 (2): 361–377.

</div>

</div>
