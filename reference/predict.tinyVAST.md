# Predict using vector autoregressive spatio-temporal model

Predicts values given new covariates using a tinyVAST model

## Usage

``` r
# S3 method for class 'tinyVAST'
predict(
  object,
  newdata,
  remove_origdata = FALSE,
  what = c("mu_g", "p_g", "p1_g", "palpha1_g", "pgamma1_g", "pepsilon1_g", "pomega1_g",
    "pdelta1_g", "pxi1_g", "p2_g", "palpha2_g", "pgamma2_g", "pepsilon2_g", "pomega2_g",
    "pdelta2_g", "pxi2_g"),
  se.fit = FALSE,
  bias.correct = FALSE,
  ...
)
```

## Arguments

- object:

  Output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/reference/tinyVAST.md).

- newdata:

  New data-frame of independent variables used to predict the response.

- remove_origdata:

  Whether to eliminate the original data from the TMB object, thereby
  speeding up the TMB object construction. However, this also eliminates
  information about random-effect variance, and is not appropriate when
  requesting predictive standard errors or epsilon bias-correction.

- what:

  What REPORTed object to output, where `mu_g` is the inverse-linked
  transformed predictor including both linear components, `p_g` is the
  sum of the first and second linear predictors (which only makes sense
  to inspect when using the Poisson-linked delta model), `p1_g` is the
  first linear predictor, `palpha_g` is the first predictor from fixed
  covariates in `formula`, `pgamma_g` is the first predictor from random
  covariates in `formula` (e.g., splines), `pomega_g` is the first
  predictor from spatial variation, `pepsilon_g` is the first predictor
  from spatio-temporal variation, `pxi_g` is the first predictor from
  spatially varying coefficients, `p2_g` is the second linear predictor,
  `palpha2_g` is the second predictor from fixed covariates in
  `formula`, `pgamma2_g` is the second predictor from random covariates
  in `formula` (e.g., splines), `pomega2_g` is the second predictor from
  spatial variation, `pepsilon2_g` is the second predictor from
  spatio-temporal variation, and `pxi2_g` is the second predictor from
  spatially varying coefficients.

- se.fit:

  Calculate standard errors?

- bias.correct:

  whether to epsilon bias-correct the predicted value

- ...:

  Not used.

## Value

Either a vector with the prediction for each row of `newdata`, or a
named list with the prediction and standard error (when
`se.fit = TRUE`).
