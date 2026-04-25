# Package index

## Fitting and diagnostics

Core tools for model fitting and diagnostics.

- [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/reference/tinyVAST.md)
  : Fit vector autoregressive spatio-temporal model
- [`tinyVASTcontrol()`](https://vast-lib.github.io/tinyVAST/reference/tinyVASTcontrol.md)
  : Control parameters for tinyVAST
- [`residuals(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/residuals.tinyVAST.md)
  : Calculate deviance or response residuals for tinyVAST
- [`logLik(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/logLik.tinyVAST.md)
  : Extract the (marginal) log-likelihood of a tinyVAST model
- [`cAIC()`](https://vast-lib.github.io/tinyVAST/reference/cAIC.md) :
  Calculate conditional AIC
- [`deviance_explained()`](https://vast-lib.github.io/tinyVAST/reference/deviance_explained.md)
  : Calculate deviance explained
- [`delta_lognormal()`](https://vast-lib.github.io/tinyVAST/reference/families.md)
  [`delta_gamma()`](https://vast-lib.github.io/tinyVAST/reference/families.md)
  : Additional families
- [`simulate(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/simulate.tinyVAST.md)
  : Simulate new data from a fitted model
- [`print(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/print.tinyVAST.md)
  : print summary of tinyVAST model
- [`sanity()`](https://vast-lib.github.io/tinyVAST/reference/sanity.md)
  : Sanity check of a tinyVAST model

## SEM and DSEM constructors

Tools to specify interactions among variables and over time.

- [`make_sem_ram()`](https://vast-lib.github.io/tinyVAST/reference/make_sem_ram.md)
  : Make a RAM (Reticular Action Model) from a SEM (structural equation
  model)
- [`make_dsem_ram()`](https://vast-lib.github.io/tinyVAST/reference/make_dsem_ram.md)
  : Make a RAM (Reticular Action Model)
- [`make_eof_ram()`](https://vast-lib.github.io/tinyVAST/reference/make_eof_ram.md)
  : Make a RAM (Reticular Action Model)

## Predicting

Core tools for model predictions.

- [`predict(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/predict.tinyVAST.md)
  : Predict using vector autoregressive spatio-temporal model
- [`integrate_output()`](https://vast-lib.github.io/tinyVAST/reference/integrate_output.md)
  : Integration for target variable
- [`sample_variable()`](https://vast-lib.github.io/tinyVAST/reference/sample_variable.md)
  : Sample from predictive distribution of a variable
- [`project()`](https://vast-lib.github.io/tinyVAST/reference/project.md)
  : Project tinyVAST to future times (EXPERIMENTAL)

## Interpret output

Tools for interpreting output.

- [`summary(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/summary.tinyVAST.md)
  : summarize tinyVAST
- [`reload_model()`](https://vast-lib.github.io/tinyVAST/reference/reload_model.md)
  : Reload a previously fitted model
- [`term_covariance()`](https://vast-lib.github.io/tinyVAST/reference/term_covariance.md)
  : Extract covariance
- [`vcov(`*`<tinyVAST>`*`)`](https://vast-lib.github.io/tinyVAST/reference/vcov.tinyVAST.md)
  : Extract Variance-Covariance Matrix
- [`rotate_pca()`](https://vast-lib.github.io/tinyVAST/reference/rotate_pca.md)
  : Rotate factors to match Principal-Components Analysis

## Stream network utilities

Tools to work with stream networks.

- [`sfnetwork_evaluator()`](https://vast-lib.github.io/tinyVAST/reference/sfnetwork_evaluator.md)
  : Construct projection matrix for stream network
- [`sfnetwork_mesh()`](https://vast-lib.github.io/tinyVAST/reference/sfnetwork_mesh.md)
  : Make mesh for stream network
- [`simulate_sfnetwork()`](https://vast-lib.github.io/tinyVAST/reference/simulate_sfnetwork.md)
  : Simulate GMRF for stream network

## Data sets

Data sets used for illustration and testing.

- [`sea_ice`](https://vast-lib.github.io/tinyVAST/reference/sea_ice.md)
  : Arctic September sea ice concentrations
- [`salmon_returns`](https://vast-lib.github.io/tinyVAST/reference/salmon_returns.md)
  : North Pacific salmon returns
- [`condition_and_density`](https://vast-lib.github.io/tinyVAST/reference/condition_and_density.md)
  : Condition and density example
- [`bering_sea_pollock_ages`](https://vast-lib.github.io/tinyVAST/reference/bering_sea_pollock_ages.md)
  : Survey catch-rates at age for Alaska pollock in the Eastern and
  Northern Bering Sea
- [`bering_sea`](https://vast-lib.github.io/tinyVAST/reference/bering_sea.md)
  : Survey domain for the eastern and northern Bering Sea surveys
- [`bering_sea_pollock_vast`](https://vast-lib.github.io/tinyVAST/reference/bering_sea_pollock_vast.md)
  : Estimated proportion-at-age for Alaska pollock using VAST
- [`red_snapper`](https://vast-lib.github.io/tinyVAST/reference/red_snapper.md)
  : Presence/absence, count, and biomass data for red snapper
- [`red_snapper_shapefile`](https://vast-lib.github.io/tinyVAST/reference/red_snapper_shapefile.md)
  : Shapefile for red snapper analysis
- [`alaska_sponge_coral_fish`](https://vast-lib.github.io/tinyVAST/reference/alaska_sponge_coral_fish.md)
  : Data to analyze sponge-coral-fish associations
- [`bering_sea_capelin_forecasts`](https://vast-lib.github.io/tinyVAST/reference/bering_sea_capelin_forecasts.md)
  : Data to demonstrate probabilistic forecasting
- [`red_grouper_diet`](https://vast-lib.github.io/tinyVAST/reference/red_grouper_diet.md)
  : Data to demonstrate model-based diet proportions
- [`atlantic_yellowtail`](https://vast-lib.github.io/tinyVAST/reference/atlantic_yellowtail.md)
  : Northwest Atlantic yellowtail
