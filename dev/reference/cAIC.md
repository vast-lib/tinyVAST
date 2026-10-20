# Calculate conditional AIC

Calculates the conditional Akaike Information criterion (cAIC).

## Usage

``` r
cAIC(object, method = c("stable", "fast"))
```

## Arguments

- object:

  Output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md).

- method:

  whether to compute effective degrees of freedom using the stable
  method `diag(solve(Hess, Hess_new))` or a faster and experimental
  method `-rowSums(sparseinv::Takahashi_Davis(-Hess) * t(Hess_new))`,
  where the latter avoids constructing a dense matrix prior to
  extracting its diagonal

## Value

cAIC value

## Details

cAIC is designed to optimize the expected out-of-sample predictive
performance for new data that share the same random effects as the
in-sample (fitted) data, e.g., spatial interpolation. In this sense, it
should be a fast approximation to optimizing the model structure based
on k-fold cross-validation.

By contrast, [`AIC()`](https://rdrr.io/r/stats/AIC.html) calculates the
marginal Akaike Information Criterion, which is designed to optimize
expected predictive performance for new data that have new random
effects, e.g., extrapolation, or inference about generative parameters.

Both cAIC and EDF are calculated using Eq. 6 of Zheng, Cadigan, and
Thorson (2024).

For models that include profiled fixed effects, these profiles are
turned off.

## References

Zheng, N., Cadigan, N., & Thorson, J. T. (2026). A unifying theory for
the evaluation of conditional Akaike information for mixed-effects
models. Journal of Statistical Computation and Simulation.
[doi:10.1080/00949655.2026.2736799](https://doi.org/10.1080/00949655.2026.2736799)

## Examples

``` r
data( red_snapper )
red_snapper = droplevels(subset(red_snapper, Data_type=="Biomass_KG"))

# Define mesh
mesh = fmesher::fm_mesh_2d( red_snapper[,c('Lon','Lat')],
                           cutoff = 1 )

# define formula with a catchability covariate for gear
formula = Response_variable ~ factor(Year) + offset(log(AreaSwept_km2))

# make variable column
red_snapper$var = "logdens"
# fit using tinyVAST
fit = tinyVAST( data = red_snapper,
                formula = formula,
                space_term = "logdens <-> logdens, sd_space",
                space_columns = c("Lon",'Lat'),
                spatial_domain = mesh,
                family = tweedie(link="log"),
                variable_column = "var",
                control = tinyVASTcontrol( getsd = FALSE,
                                           profile = "alpha_j" ) )

cAIC(fit) # conditional AIC
#> Error in UseMethod("cAIC", object): no applicable method for 'cAIC' applied to an object of class "tinyVAST"
AIC(fit) # marginal AIC
#> [1] 43073.52
```
