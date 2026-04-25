# Control parameters for tinyVAST

Control parameters for tinyVAST

## Usage

``` r
tinyVASTcontrol(
  opt_loops = 1,
  newton_loops = 0,
  eval.max = 10000,
  iter.max = 10000,
  getsd = TRUE,
  silent = getOption("tinyVAST.silent", TRUE),
  trace = getOption("tinyVAST.trace", 0),
  verbose = getOption("tinyVAST.verbose", FALSE),
  profile = c(),
  tmb_par = NULL,
  tmb_map = NULL,
  tmb_random = NULL,
  gmrf_parameterization = c("separable", "projection"),
  reml = FALSE,
  getJointPrecision = FALSE,
  calculate_deviance_explained = TRUE,
  run_model = TRUE,
  suppress_user_warnings = FALSE,
  get_rsr = FALSE,
  extra_reporting = FALSE,
  use_anisotropy = FALSE,
  sar_adjacency = "queen",
  barrier_stiffness = 0.01
)
```

## Arguments

- opt_loops:

  Integer number of times to call nonlinear optimizer.

- newton_loops:

  Integer number of Newton steps to do after running
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- eval.max:

  Maximum number of evaluations of the objective function allowed.
  Passed to `control` in
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- iter.max:

  Maximum number of iterations allowed. Passed to `control` in
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- getsd:

  Boolean indicating whether to call
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html)

- silent:

  Disable terminal output for inner optimizer?

- trace:

  Parameter values are printed every `trace` iteration for the outer
  optimizer. Passed to `control` in
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- verbose:

  Output additional messages about model steps during fitting?

- profile:

  Character-vector passed to
  [TMB::MakeADFun](https://rdrr.io/pkg/TMB/man/MakeADFun.html) and see
  description there. Fixed effects that are highly correlated with
  random effects can often be estimated faster (i.e., with fewer
  iterations) by adding them to `profile`. The most common use-case is
  `profile = c("alpha_j","alpha2_j")`. However, doing so will have a
  small impact on model estimates and predictions.

- tmb_par:

  named list of parameters used as starting values. Elements that have
  names that match those constructed internally then replace the
  internally constructed starting values. Those that match must have
  identical shape to `tinyVAST(...)$internal$parlist`

- tmb_map:

  input passed to
  [TMB::MakeADFun](https://rdrr.io/pkg/TMB/man/MakeADFun.html) as
  argument `map`, over-writing the version
  `tinyVAST(...)$tmb_inputs$tmb_map` and allowing detailed control over
  estimated parameters (advanced feature)

- tmb_random:

  input passed to
  [TMB::MakeADFun](https://rdrr.io/pkg/TMB/man/MakeADFun.html) as
  argument `random`, over-writing the version
  `tinyVAST(...)$tmb_inputs$tmb_random` and allowing detailed control
  over parameters treated as random effects. `tmb_random = NULL` uses
  the default (internal) construction, and use `tmb_random = c("empty")`
  (or some other name that doesn't match actual parameters) to use
  penalized likelihood (presumably while also modifying `tmb_map` and
  `tmb_par`)

- gmrf_parameterization:

  Parameterization to use for the Gaussian Markov random field, where
  the default `separable` constructs a full-rank and separable precision
  matrix, and the alternative `projection` constructs a full-rank and
  IID precision for variables over time, and then projects this using
  the inverse-cholesky of the precision, where this projection allows
  for rank-deficient covariance.

- reml:

  Logical: use REML (restricted maximum likelihood) estimation rather
  than maximum likelihood? Internally, this adds the fixed effects to
  the list of random effects to integrate over.

- getJointPrecision:

  whether to get the joint precision matrix. Passed to
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).

- calculate_deviance_explained:

  whether to calculate proportion of deviance explained. See
  [`deviance_explained()`](https://vast-lib.github.io/tinyVAST/dev/reference/deviance_explained.md)

- run_model:

  whether to run the model of export TMB objects prior to compilation
  (useful for debugging)

- suppress_user_warnings:

  whether to suppress warnings from package author regarding dangerous
  or non-standard options

- get_rsr:

  Experimental option, whether to report restricted spatial regression
  (RSR) adjusted estimator for covariate responses

- extra_reporting:

  Whether to report a much larger set of quantities via
  `obj$env$report()`

- use_anisotropy:

  Whether to estimate two parameters representing geometric anisotropy

- sar_adjacency:

  Whether to use queen or rook adjacency when defining a Simultaneous
  Autoregressive spatial precision from a sfc_GEOMETRY (default is
  queen)

- barrier_stiffness:

  The ratio of local stiffness (the scale of diffusion rate and
  resulting decorrelation distance) for barriers relative to normal
  areas in the SPDE method when using `add_mesh_covariates`. The default
  `barrier_stiffness = 0.01` is the value from Bakka et al. 2019.

## Value

An object (list) of class `tinyVASTcontrol`, containing either default
or updated values supplied by the user for model settings

## References

Bakka, H., Vanhatalo, J., Illian, J., Simpson, D., Rue, H. (2019).
Non-stationary Gaussian models with physical barriers. Spatial
Statistics, 29, 268-288.
[doi:10.1016/j.spasta.2019.01.002](https://doi.org/10.1016/j.spasta.2019.01.002)
