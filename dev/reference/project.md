# Project tinyVAST to future times (EXPERIMENTAL)

Projects a fitted model forward in time.

## Usage

``` r
project(
  object,
  extra_times,
  newdata,
  what = "mu_g",
  future_var = TRUE,
  past_var = FALSE,
  parm_var = FALSE
)
```

## Arguments

- object:

  fitted model from `tinyVAST(.)`

- extra_times:

  a vector of extra times, matching values in `newdata`

- newdata:

  data frame including new values for `time_variable`

- what:

  What REPORTed object to output, where `mu_g` is the inverse-linked
  transformed predictor including both linear components, `p_g` is the
  sum of the first and second linear predictors (which only makes sense
  to inspect when using the Poisson-linked delta model), `p1_g` is the
  first linear predictor, `palpha1_g` is the first predictor from fixed
  covariates in `formula`, `pgamma1_g` is the first predictor from
  random covariates in `formula` (e.g., splines), `pomega1_g` is the
  first predictor from spatial variation, `pepsilon1_g` is the first
  predictor from spatio-temporal variation, `pxi1_g` is the first
  predictor from spatially varying coefficients, `p2_g` is the second
  linear predictor, `palpha2_g` is the second predictor from fixed
  covariates in `formula`, `pgamma2_g` is the second predictor from
  random covariates in `formula` (e.g., splines), `pomega2_g` is the
  second predictor from spatial variation, `pepsilon2_g` is the second
  predictor from spatio-temporal variation, and `pxi2_g` is the second
  predictor from spatially varying coefficients.

- future_var:

  logical indicating whether to simulate future process errors from
  GMRFs, or just compute the predictive mean

- past_var:

  logical indicating whether to re-simulate past process errors from
  predictive distribution of random effects, thus changing the boundary
  condition of the forecast

- parm_var:

  logical indicating whether to re-sample fixed effects from their
  predictive distribution, thus changing the GMRF for future process
  errors

## Value

A vector of values corresponding to rows in `newdata`

## Examples

``` r
# Convert to long-form
set.seed(123)
n_obs = 100
rho = 0.9
sigma_x = 0.2
sigma_y = 0.1
x = rnorm(n_obs, mean=0, sd = sigma_x)
for(i in 2:length(x)) x[i] = rho * x[i-1] + x[i]
y = x + rnorm( length(x), mean = 0, sd = sigma_y )
data = data.frame( "val" = y, "var" = "y", "time" = seq_along(y) )

# Define AR2 time_term
time_term = "
  y -> y, 1, rho1
  y -> y, 2, rho2
  y <-> y, 0, sd
"

# fit model
mytiny = tinyVAST(
  time_term = time_term,
  data = data,
  times = unique(data$t),
  variables = "y",
  formula = val ~ 1,
  control = tinyVASTcontrol( getJointPrecision = TRUE )
)

# Deterministic projection
extra_times = length(x) + 1:100
n_sims = 10
newdata = data.frame( "time" = c(seq_along(x),extra_times), "var" = "y" )
Y = project(
  mytiny,
  newdata = newdata,
  extra_times = extra_times,
  future_var = FALSE
)
#> Error: object does not inherit from class sdmTMB
plot( x = seq_along(Y),
      y = Y,
      type = "l", lty = "solid", col = "black" )
#> Error: object 'Y' not found

# Stochastic projection with future process errors
if (FALSE) { # \dontrun{
extra_times = length(x) + 1:100
n_sims = 10
newdata = data.frame( "time" = c(seq_along(x),extra_times), "var" = "y" )
Y = NULL
for(i in seq_len(n_sims) ){
  tmp = project(
    mytiny,
    newdata = newdata,
    extra_times = extra_times,
    future_var = TRUE,
    past_var = TRUE,
    parm_var = TRUE
  )
  Y = cbind(Y, tmp)
}
matplot( x = row(Y),
         y = Y,
         type = "l", lty = "solid", col = "black" )
} # }
```
