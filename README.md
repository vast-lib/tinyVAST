
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tinyVAST

> Multivariate spatio-temporal models using dynamic structural equations

<!-- badges: start -->
[![R-CMD-check](https://github.com/James-Thorson-NOAA/tinyVAST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/James-Thorson-NOAA/tinyVAST/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/vast-lib/tinyVAST/branch/main/graph/badge.svg)](https://app.codecov.io/gh/vast-lib/tinyVAST?branch=main)
[![Documentation](https://img.shields.io/badge/documentation-tinyVAST-orange.svg?colorB=E91E63)](https://vast-lib.github.io/tinyVAST/)
[![DOI](https://zenodo.org/badge/704919773.svg)](https://doi.org/10.5281/zenodo.15001856)
[![](https://www.r-pkg.org/badges/version/tinyVAST)](https://cran.r-project.org/package=tinyVAST)
[![](https://cranlogs.r-pkg.org/badges/tinyVAST)](https://cran.r-project.org/package=tinyVAST)
[![](https://cranlogs.r-pkg.org/badges/grand-total/tinyVAST)](https://cran.r-project.org/package=tinyVAST)
<!-- badges: end -->

tinyVAST is an R package that fits multivariate spatio-temporal models while using Gaussian Markov random fields to represent nonseparable interactions among variables over time. See a preprint:

Thorson, J. T., Anderson, S. C., Goddard, P., & Rooper, C. N. (2024). tinyVAST: R package with an expressive interface to specify lagged and simultaneous effects in multivariate spatio-temporal models (arXiv:2401.10193). arXiv. http://arxiv.org/abs/2401.10193

## Table of contents

- [Installation](#installation)
- [Citation](#citation)
- [Related software](#related-software)

## Installation

tinyVAST can be installed from GitHub:

``` r
library(devtools)
install_github("vast-lib/tinyVAST", dependencies = TRUE)
```

## Citation

To cite tinyVAST in publications use:

``` r
citation("tinyVAST")
```

Thorson, J. T., Anderson, S. C., Goddard, P., & Rooper, C. N. (2025).
tinyVAST: R package with an expressive interface to specify lagged and
simultaneous effects in multivariate spatio-temporal models. 
Global Ecology and Biogeography.  <https://doi.org/10.1111/geb.70035>

## Related software

tinyVAST is builds upon many packages. This includes
[VAST](https://github.com/James-Thorson-NOAA/VAST) R package:

Thorson, J.T. 2019. Guidance for decisions using the Vector
Autoregressive Spatio-Temporal (VAST) package in stock, ecosystem,
habitat and climate assessments. Fisheries Research 210: 143–161.
<https://doi.org/10.1016/j.fishres.2018.10.013>.

and [sdmTMB](https://github.com/pbs-assess/sdmTMB) R package:

Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett. 2022. sdmTMB:
an R package for fast, flexible, and user-friendly generalized linear
mixed effects models with spatial and spatiotemporal random fields.
bioRxiv 2022.03.24.485545; doi:
<https://doi.org/10.1101/2022.03.24.485545>

and the [glmmTMB](https://github.com/glmmTMB/glmmTMB) R package:

Brooks, M.E., Kristensen, K., van Benthem, K.J., Magnusson, A., Berg,
C.W., Nielsen, A., Skaug, H.J., Maechler, M., and Bolker, B.M. 2017.
glmmTMB balances speed and flexibility among packages for zero-inflated
generalized linear mixed modeling. The R Journal 9(2): 378–400.
<https://doi.org/10.32614/rj-2017-066>.

[INLA](https://www.r-inla.org/) and
[inlabru](https://sites.google.com/inlabru.org/inlabru) can fit many of
the same models as sdmTMB (and many more) in an approximate Bayesian
inference framework.

[mgcv](https://cran.r-project.org/package=mgcv) can fit similar
SPDE-based Gaussian random field models with code included in [Miller et
al. (2019)](https://doi.org/10.1007/s13253-019-00377-z).
