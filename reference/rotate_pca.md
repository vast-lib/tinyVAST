# Rotate factors to match Principal-Components Analysis

Rotate lower-triangle loadings matrix to order factors from largest to
smallest variance.

## Usage

``` r
rotate_pca(
  L_tf,
  x_sf = matrix(0, nrow = 0, ncol = ncol(L_tf)),
  order = c("none", "increasing", "decreasing")
)
```

## Arguments

- L_tf:

  Loadings matrix with dimension \\T \times F\\.

- x_sf:

  Spatial response with dimensions \\S \times F\\.

- order:

  Options for resolving label-switching via reflecting each factor to
  achieve a given order across dimension \\T\\.

## Value

List containing the rotated loadings `L_tf`, the inverse-rotated
response matrix `x_sf`, and the rotation `H`
