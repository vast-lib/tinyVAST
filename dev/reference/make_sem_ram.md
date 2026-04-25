# Make a RAM (Reticular Action Model) from a SEM (structural equation model)

`make_sem_ram` converts SEM arrow notation to `ram` describing SEM
parameters

## Usage

``` r
make_sem_ram(sem, variables, quiet = FALSE, covs = variables)
```

## Arguments

- sem:

  structural equation model structure, passed to either
  [`specifyModel`](https://rdrr.io/pkg/sem/man/specifyModel.html) or
  [`specifyEquations`](https://rdrr.io/pkg/sem/man/specifyModel.html)
  and then parsed to control the set of path coefficients and
  variance-covariance parameters

- variables:

  A character vector listing the set of variables

- quiet:

  if `FALSE`, the default, then the number of input lines is reported
  and a message is printed suggesting that `specifyEquations` or `cfa`
  be used.

- covs:

  optional: a character vector of one or more elements, with each
  element giving a string of variable names, separated by commas.
  Variances and covariances among all variables in each such string are
  added to the model. For confirmatory factor analysis models specified
  via `cfa`, `covs` defaults to all of the factors in the model, thus
  specifying all variances and covariances among these factors.
  *Warning*: `covs="x1, x2"` and `covs=c("x1", "x2")` are *not*
  equivalent: `covs="x1, x2"` specifies the variance of `x1`, the
  variance of `x2`, *and* their covariance, while `covs=c("x1", "x2")`
  specifies the variance of `x1` and the variance of `x2` *but not*
  their covariance.

## Value

An S3-class `"sem_ram"` containing:

- `model`:

  Output from
  [`specifyEquations`](https://rdrr.io/pkg/sem/man/specifyModel.html) or
  [`specifyModel`](https://rdrr.io/pkg/sem/man/specifyModel.html) that
  defines paths and parameters

- `ram`:

  reticular action module (RAM) describing dependencies
