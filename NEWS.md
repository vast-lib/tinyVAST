# tinyVAST 0.8.0

* Adding nbinom1 and nbinom2 families

# tinyVAST 0.7.1

* Fixed bug (wrong output) when using `predict(fit, what="mu_g")` and
  a Poisson-linked delta model
* Fixed bug (cryptic error message) when using `integrate_output`
* Add `cAIC` (but disabling EDF calculation for now)

# tinyVAST 0.7.0

* Adding option for spatially-varying-coefficient (SVC) models
* Add error-check for when `data` has a factor with extra levels, which
  conflicted with the logic of adding all `origdata` levels to `newdata`
  when calling `predict`, and hence caused an uniformative error previously

# tinyVAST 0.6.0

* Change `integrate_output` interface by splitting `W_gz` and `V_gz`
  into four vectors `area`, `type`, `covariate`, and `weighting_index`
  to simplify documentations and improve naming 
* Fix bug where cloglog and logit links were not previously implemented 
  for use in `predict` and `integrate_output`

# tinyVAST 0.5.0

* Adding vignette showing how to fit multiple data types in an SDM
* Adding `deviance_explained` and calculating this by default

# tinyVAST 0.4.0

* Adding code for simulation residuals, and examples to vignettes

# tinyVAST 0.3.0

* Adding `sdmTMB` as dependency, and importing family options from it
* Adding vignette for joint analysis of condition and density

# tinyVAST 0.2.0

* Add option to specify covariance in SEM and DSEM notation

# tinyVAST 0.1.0

* Initial public alpha release
