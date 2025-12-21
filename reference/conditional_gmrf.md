# Conditional simulation from a GMRF

Generates samples from a Gaussian Markov random field (GMRF) conditional
upon fixed values for some elements.

## Usage

``` r
conditional_gmrf(
  Q,
  observed_idx,
  x_obs,
  n_sims = 1,
  what = c("simulate", "predict")
)
```

## Arguments

- Q:

  precision for a zero-centered GMRF.

- observed_idx:

  integer vector listing rows of `Q` corresponding to fixed measurements

- x_obs:

  numeric vector with fixed values for indices `observed_idx`

- n_sims:

  integer listing number of simulated values

- what:

  Whether to simulate from the conditional GMRF, or predict the mean and
  precision

## Value

A matrix with `n_sims` columns and a row for every row of `Q` not in
`observed_idx`, with simulations for those rows
