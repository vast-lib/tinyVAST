# Sanity check of a tinyVAST model

Sanity check of a tinyVAST model

## Usage

``` r
sanity(object, big_sd_log10 = 2, gradient_thresh = 0.001, silent = FALSE)
```

## Arguments

- object:

  Fitted model from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md).

- big_sd_log10:

  Value to check size of standard errors against. A value of 2 would
  indicate that standard errors greater than `10^2` (i.e., 100) should
  be flagged.

- gradient_thresh:

  Gradient threshold to issue warning.

- silent:

  Logical: suppress messages? Useful to set to `TRUE` if running large
  numbers of models and just interested in returning sanity list
  objects.

## Value

An invisible named list of checks.
