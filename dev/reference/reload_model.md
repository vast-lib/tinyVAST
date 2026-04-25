# Reload a previously fitted model

`reload_model` allows a user to save a fitted model, reload it in a new
R terminal, and then relink the DLLs so that it functions as expected.

## Usage

``` r
reload_model(x, check_gradient = TRUE)
```

## Arguments

- x:

  Output from
  [`tinyVAST`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md),
  potentially with DLLs not linked

- check_gradient:

  Whether to check the gradients of the reloaded model

## Value

Output from
[`tinyVAST`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md)
with DLLs relinked
