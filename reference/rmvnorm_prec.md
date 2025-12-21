# Multivariate Normal Random Deviates using Sparse Precision

This function provides a random number generator for the multivariate
normal distribution with mean equal to `mu` and sparse precision matrix
`prec`.

## Usage

``` r
rmvnorm_prec(prec, n = 1, mu = rep(0, nrow(prec)))
```

## Arguments

- prec:

  sparse precision (inverse-covariance) matrix.

- n:

  number of observations.

- mu:

  mean vector.

## Value

a matrix with dimension `length(mu)` by `n`, containing realized draws
from the specified mean and precision
