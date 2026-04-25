# summarize tinyVAST

summarize parameters from a fitted tinyVAST

## Usage

``` r
# S3 method for class 'tinyVAST'
summary(
  object,
  what = c("space_term", "time_term", "spacetime_term", "fixed"),
  predictor = c("one", "two"),
  ...
)
```

## Arguments

- object:

  Output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md)

- what:

  What component to summarize, whether `space_term`, `spacetime_term`,
  or `fixed` for the fixed effects included in the GAM formula

- predictor:

  whether to get the 1st or 2nd linear predictor (the latter is only
  applicable in delta models)

- ...:

  Not used

## Value

A data-frame containing the estimate (and standard errors, two-sided
Wald-test z-value, and associated p-value if the standard errors are
available) for model parameters, including the fixed-effects specified
via `formula`, or the path coefficients for the spatial SEM specified
via `space_term`, the dynamic SEM specified via `time_term`, or the
spatial dynamic SEM specified via `spacetime_term`

## Details

`tinyVAST` includes three components:

- Space-variable interaction:

  a separable Gaussian Markov random field (GMRF) constructed from a
  structural equation model (SEM) and a spatial variable

- Space-variable-time interaction:

  a separable GMRF constructed from a a dynamic SEM (a nonseparable
  time-variable interaction) and a spatial variable

- Additive variation:

  a generalized additive model (GAM), representing exogenous covariates

Each of these are summarized and interpreted differently, and
`summary.tinyVAST` facilitates this.

Regarding the DSEM componennt, tinyVAST includes an "arrow and lag"
notation, which specifies the set of path coefficients and exogenous
variance parameters to be estimated. Function `tinyVAST` then estimates
the maximum likelihood value for those coefficients and parameters by
maximizing the log-marginal likelihood.

However, many users will want to associate individual parameters and
standard errors with the path coefficients that were specified using the
"arrow and lag" notation. This task is complicated in models where some
path coefficients or variance parameters are specified to share a single
value a priori, or were assigned a name of NA and hence assumed to have
a fixed value a priori (such that these coefficients or parameters have
an assigned value but no standard error). The `summary` function
therefore compiles the MLE for coefficients (including duplicating
values for any path coefficients that assigned the same value) and
standard error estimates, and outputs those in a table that associates
them with the user-supplied path and parameter names. It also outputs
the z-score and a p-value arising from a two-sided Wald test (i.e.
comparing the estimate divided by standard error against a standard
normal distribution).
