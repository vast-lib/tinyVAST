# Construct projection matrix for stream network

Make sparse matrix to project from stream-network nodes to user-supplied
points

## Usage

``` r
sfnetwork_evaluator(stream, loc, tolerance = 0.01)
```

## Arguments

- stream:

  sfnetworks object representing stream network

- loc:

  sf object representing points to which are being projected

- tolerance:

  error-check tolerance

## Value

the sparse interpolation matrix, with rows for each row of `data`
supplied during fitting and columns for each spatial random effect.
