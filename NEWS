# tinyVAST 1.2.0

* Adding a new spatial domain using `sf::st_make_grid`
* Adding option for geometric anisotropy when using `sf::st_make_grid`

# tinyVAST 1.1.1

* Adding option for geometric anisotropy when using the SPDE method
* For some reason, the CPP edits also address an error message
  during `devtools::check_win_devel` "array subscript 'const __m128i[0]' is
  partly outside array bounds of 'unsigned char [12]'"

# tinyVAST 1.1.0

* Adding `term_covariance` to calculate the covariance among variables for SEM term,
  or covariance among variables-and-lags for DSEM terms
* Adding Restricted Spatial Regression estimator for covariates `alphaprime_j` and `alphaprime2_j`
  to `fit$rep` output.
* Adding methods to allow use of `cv` to calculate crossvalidation skill
* Add `bias.correct` option to predict (but still no flag for SEs for anything except p_i)
* Switch `sample_variable` from using `obj$env$MC` to `obj$env$spHess(random=TRUE)`
  which seems more stable as dependency
* Add functionality for `te` and `ti` splines, although they remain poorly tested
* Add error check to `sfnetwork_mesh` to detect if the stream network is not ordered
  as a tree
* Improve stream network vignette to use matrix notation for joint precision, and
  modify `simulate_sfnetwork` to use that
* Change `tinyVAST.cpp` to use matrix notation constructor and fix bug in previous
  constructor where the covariance between first and second nodes was not right
* Expand `test-sfnetworks.R` integrated test to confirm that matrix-notation
  precision constructor is identical to the inverse of Ornstein-Uhlenbeck covariance
  as intended.

# tinyVAST 1.0.1

* Modify examples for `simulate.tinyVAST` and `sample_variable` to try to avoid
  terminal output giving error in valgrind check
* Add `ivector_minus_one` function to satisfy clang-UBSAN
* Swap `GMRF(Q).Quadform(x)` to `x.matrix().transpose() * (Q * x.matrix())` to
  avoid calculating log-determinant of `Q` in the smoothers penalty, to avoid
  valgrind errors
* Add tinyVASTcontrol option that suppresses nlminb warning messages by default,
  which are typically not informative to casual users

# tinyVAST 1.0.0

* Adding nbinom1 and nbinom2 families
* Simplify argument names, by changing `sem` to `space_term`, `dsem` to `spacetime_term`
  and `spatial_graph` to `spatial_domain`, and eliminating `delta_` in the names
  for arguments to `delta_options`
* Add `time_term` to allow time-variable interaction (e.g., AR1 intercepts)
* Adding overview and model-description vignettes
* Add `simulate` S3 generic

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
