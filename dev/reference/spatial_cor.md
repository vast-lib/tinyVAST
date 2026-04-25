# Approximate spatial correlation

Extract the approximated spatial correlation between one coordinate and
other coordinates using a sparse precision and SPDE mesh

## Usage

``` r
spatial_cor(Q, mesh, coord, pred)
```

## Arguments

- Q:

  sparse precision matrix

- mesh:

  SPDE mesh

- coord:

  vector of length-2 with spatial coordinates for focal point

- pred:

  matrix with two columns and multiple rows, with location for points to
  predict correlation

## Value

A vector with length `nrow(pred)` giving the spatial correlation
