# Sample from predictive distribution of a variable

`sample_variable` samples from the joint distribution of random (and
optionally fixed) effects to approximate the predictive distribution for
a variable.

## Usage

``` r
sample_variable(
  object,
  newdata = NULL,
  variable_name = "mu_i",
  n_samples = 100,
  sample_fixed = TRUE,
  seed = 123456
)
```

## Arguments

- object:

  output from `\code{tinyVAST()}`

- newdata:

  data frame of new data, used to sample model components for
  predictions e.g., `mu_g`

- variable_name:

  name of variable available in report using `Obj$report()` or
  parameters using `Obj$env$parList()`

- n_samples:

  number of samples from the joint predictive distribution for fixed and
  random effects. Default is 100, which is slow.

- sample_fixed:

  whether to sample fixed and random effects, `sample_fixed=TRUE` as by
  default, or just sample random effects, `sample_fixed=FALSE`

- seed:

  integer used to set random-number seed when sampling variables, as
  passed to `set.seed(.)`

## Value

A matrix with a row for each `data` supplied during fitting, and
`n_samples` columns, where each column in a vector of samples for a
requested quantity given sampled uncertainty in fixed and/or random
effects

## Details

Using `sample_fixed=TRUE` (the default) in `sample_variable` propagates
variance in both fixed and random effects, while using
`sample_fixed=FALSE` does not. Sampling fixed effects will sometimes
cause numerical under- or overflow (i.e., output values of `NA`) in
cases when variance parameters are estimated imprecisely. In these
cases, the multivariate normal approximation being used is a poor
representation of the tail probabilities, and results in some samples
with implausibly high (or negative) variances, such that the associated
random effects then have implausibly high magnitude.

## Examples

``` r
 set.seed(101)
 x = runif(n = 100, min = 0, max = 2*pi)
 y = 1 + sin(x) + 0.1 * rnorm(100)

 # Do fit with getJointPrecision=TRUE
 fit = tinyVAST( formula = y ~ s(x),
                 data = data.frame(x=x,y=y) )

 # samples from distribution for the mean
 # excluding fixed effects due to CRAN checks
 samples = sample_variable(fit, sample_fixed = FALSE)
#> # Obtaining samples from predictive distribution for variable mu_i
#>   Finished sample 10 of 100
#>   Finished sample 20 of 100
#>   Finished sample 30 of 100
#>   Finished sample 40 of 100
#>   Finished sample 50 of 100
#>   Finished sample 60 of 100
#>   Finished sample 70 of 100
#>   Finished sample 80 of 100
#>   Finished sample 90 of 100
#>   Finished sample 100 of 100
```
