# Calculate deviance explained

`deviance_explained` fits a null model, calculates the deviance relative
to a saturated model for both the original and the null model, and uses
these to calculate the proportion of deviance explained.

This implementation conditions upon the maximum likelihood estimate of
fixed effects and the empirical Bayes ("plug-in") prediction of random
effects. It can be described as "conditional deviance explained". A
state-space model that estimates measurement error variance approaching
zero (i.e., collapses to a process-error-only model) will have a
conditional deviance explained that approaches 1.0

For several families (tweedie, negbin1, negbin2, and student), the null
model is fitted using the MLE for an overdispersion parameter from the
full model. This is done because, e.g., the negbin1 and negbin2 only
belong to the exponential family when the overdispersion parameter is
fixed, and the deviance relative to a saturated model is only defined
for the exponential family.

## Usage

``` r
deviance_explained(x, null_formula, null_delta_formula = ~1)
```

## Arguments

- x:

  output from `\code{tinyVAST()}`

- null_formula:

  formula for the null model. If missing, it uses
  `null_formula = response ~ 1`. For multivariate models, it might make
  sense to use `null_formula = response ~ category`

- null_delta_formula:

  formula for the null model for the delta component. If missing, it
  uses `null_formula = response ~ 1`. For multivariate models, it might
  make sense to use `null_delta_formula = response ~ category`

## Value

the proportion of conditional deviance explained.
