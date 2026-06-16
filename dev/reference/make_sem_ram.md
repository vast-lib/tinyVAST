# Make a RAM (Reticular Action Model) from a SEM (structural equation model)

`make_sem_ram` converts SEM arrow notation to `ram` describing SEM
parameters

## Usage

``` r
make_sem_ram(sem, variables, quiet = FALSE, covs = variables)
```

## Arguments

- sem:

  structural equation model structure, passed to either `specifyModel`
  or `specifyEquations` and then parsed to control the set of path
  coefficients and variance-covariance parameters

- variables:

  A character vector listing the set of variables

- quiet:

  Boolean indicating whether to print messages to terminal

- covs:

  A character vector listing variables for which to estimate a standard
  deviation

## Value

An S3-class `"sem_ram"` containing:

- `model`:

  Output from
  [`specifyEquations`](https://rdrr.io/pkg/sem/man/specifyModel.html) or
  [`specifyModel`](https://rdrr.io/pkg/sem/man/specifyModel.html) that
  defines paths and parameters

- `ram`:

  reticular action module (RAM) describing dependencies
