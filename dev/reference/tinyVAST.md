# Fit vector autoregressive spatio-temporal model

Fits a vector autoregressive spatio-temporal (VAST) model using a
minimal feature-set and a widely used interface.

## Usage

``` r
tinyVAST(
  formula,
  data,
  time_term = NULL,
  space_term = NULL,
  spacetime_term = NULL,
  family = gaussian(),
  delta_options = list(formula = ~1),
  spatial_varying = NULL,
  weights = NULL,
  spatial_domain = NULL,
  development = list(),
  control = tinyVASTcontrol(),
  space_columns = c("x", "y"),
  time_column = "time",
  times = NULL,
  variable_column = "var",
  variables = NULL,
  distribution_column = "dist"
)
```

## Arguments

- formula:

  Formula with response on left-hand-side and predictors on
  right-hand-side, parsed by `mgcv` and hence allowing `s(.)` for
  splines or `offset(.)` for an offset.

- data:

  Data-frame of predictor, response, and offset variables. Also includes
  variables that specify space, time, variables, and the distribution
  for samples, as identified by arguments `variable_column`,
  `time_column`, `space_columns`, and `distribution_column`.

- time_term:

  Specification for time-series structural equation model structure for
  constructing a time-variable interaction that defines a time-varying
  intercept for each variable (i.e., applies uniformly across space).
  `time_term=NULL` disables the space-variable interaction; see
  [`make_dsem_ram()`](https://vast-lib.github.io/tinyVAST/dev/reference/make_dsem_ram.md)
  for notation.

- space_term:

  Specification for structural equation model structure for constructing
  a space-variable interaction. `space_term=NULL` disables the
  space-variable interaction; see
  [`make_sem_ram()`](https://vast-lib.github.io/tinyVAST/dev/reference/make_sem_ram.md)
  for notation.

- spacetime_term:

  Specification for time-series structural equation model structure
  including lagged or simultaneous effects for constructing a
  time-variable interaction, which is then combined in a separable
  process with the spatial correlation to form a space-time-variable
  interaction (i.e., the interaction occurs locally at each site).
  `spacetime_term=NULL` disables the space-variable interaction; see
  [`make_dsem_ram()`](https://vast-lib.github.io/tinyVAST/dev/reference/make_dsem_ram.md)
  or
  [`make_eof_ram()`](https://vast-lib.github.io/tinyVAST/dev/reference/make_eof_ram.md).

- family:

  A function returning a class `family`, including
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.html),
  [`tweedie()`](https://sdmTMB.github.io/sdmTMB/reference/families.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`Gamma()`](https://rdrr.io/r/stats/family.html),
  [`student()`](https://sdmTMB.github.io/sdmTMB/reference/families.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`nbinom1()`](https://sdmTMB.github.io/sdmTMB/reference/families.html),
  or
  [`nbinom2()`](https://sdmTMB.github.io/sdmTMB/reference/families.html).
  Alternatively, can be a named list of these functions, with names that
  match levels of `data$distribution_column` to allow different families
  by row of data. Delta model families are possible, and see `Families`
  for delta-model options. For binomial family options, see 'Binomial
  families' in the Details section below.

- delta_options:

  a named list with slots for `formula`, `space_term`, `spacetime_term`,
  `time_term`, and `spatial_varying`. These specify options for the
  second linear predictor of a delta model, and are only used (or
  estimable) when a `delta family` is used for some samples.

- spatial_varying:

  a formula specifying spatially varying coefficients (SVC). Note that
  using formulas in R, `spatial_varying = ~ X` automatically adds an
  intercept to implicitly read as `spatial_varying = ~ 1 + X`, so
  tinyVAST then estimates an SVC for an intercept in addition to
  covariate `X`. Therefore, if you only want an SVC for a single
  covariate, use `spatial_varying = ~ 0 + X` to suppress the default
  behavior of formulas in R.

- weights:

  A numeric vector representing optional likelihood weights for the data
  likelihood. Weights do not have to sum to one and are not internally
  modified. Thee weights argument needs to be a vector and not a name of
  the variable in the data frame.

- spatial_domain:

  Object that represents spatial relationships, either using
  [`fmesher::fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
  to apply the SPDE method,
  [`igraph::make_empty_graph()`](https://r.igraph.org/reference/make_empty_graph.html)
  for independent time-series,
  [`igraph::make_graph()`](https://r.igraph.org/reference/make_graph.html)
  to apply a simultaneous autoregressive (SAR) process to a
  user-supplied graph,
  [`sfnetwork_mesh()`](https://vast-lib.github.io/tinyVAST/dev/reference/sfnetwork_mesh.md)
  for stream networks, or class `sf` or `sfc` containing polygons e.g
  constructed using
  [sf::st_make_grid](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  to apply a SAR to an areal model with adjacency based on the geometry
  of the object, or `NULL` to specify a single site. If using `igraph`
  then the graph must have vertex names `V(graph)$name` that match
  levels of `data[,'space_columns']`.

- development:

  Specify options that are under active development. Please do not use
  these features without coordinating with the package authors.

- control:

  Output from
  [`tinyVASTcontrol()`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVASTcontrol.md),
  used to define user settings.

- space_columns:

  A string or character vector that indicates the column(s) of `data`
  indicating the location of each sample. When `spatial_domain` is an
  `igraph` object, `space_columns` is a string with with levels matching
  the names of vertices of that object. When `spatial_domain` is an
  `fmesher` or `sfnetwork` object, space_columns is a character vector
  indicating columns of `data` with coordinates for each sample.

- time_column:

  A character string indicating the column of `data` listing the
  time-interval for each sample, from the set of times in argument
  `times`.

- times:

  A integer vector listing the set of times in order. If `times=NULL`,
  then it is filled in as the vector of integers from the minimum to
  maximum value of `data$time`. Alternatively, it could be the minimum
  value of `data$time` through future years, such that the model can
  forecast those future years.

- variable_column:

  A character string indicating the column of `data` listing the
  variable for each sample, from the set of times in argument
  `variables`.

- variables:

  A character vector listing the set of variables. if `variables=NULL`,
  then it is filled in as the unique values from
  `data$variable_columns`.

- distribution_column:

  A character string indicating the column of `data` listing the
  distribution for each sample, from the set of names in argument
  `family`. if `variables=NULL`, then it is filled in as the unique
  values from `data$variables`.

## Value

An object (list) of class `tinyVAST`. Elements include:

- data:

  Data-frame supplied during model fitting

- spatial_domain:

  the spatial domain supplied during fitting

- formula:

  the formula specified during model fitting

- obj:

  The TMB object from
  [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- opt:

  The output from [`nlminb`](https://rdrr.io/r/stats/nlminb.html)

- opt:

  The report from `obj$report()`

- sdrep:

  The output from
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html)

- tmb_inputs:

  The list of inputs passed to
  [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- call:

  A record of the function call

- run_time:

  Total time to run model

- interal:

  Objects useful for package function, i.e., all arguments passed during
  the call

- deviance_explained:

  output from
  [`deviance_explained`](https://vast-lib.github.io/tinyVAST/dev/reference/deviance_explained.md)

## Details

`tinyVAST` includes several basic inputs that specify the model
structure:

- `formula` specifies covariates and splines in a Generalized Additive
  Model;

- `time_term` specifies interactions among variables and over time that
  are constant across space, constructing the time-variable interaction.

- `space_term` specifies interactions among variables and over time that
  occur based on the variable values at each location, constructing the
  space-variable interaction.

- `spacetime_term` specifies interactions among variables and over time,
  constructing the space-time-variable interaction.

These inputs require defining the *domain* of the model. This includes:

- `spatial_domain` specifies spatial domain, with determines spatial
  correlations

- `times` specifies the temporal domain, i.e., sequence of time-steps

- `variables` specifies the set of variables, i.e., the variables that
  will be modeled

The default `spacetime_term=NULL` and `space_term=NULL` turns off all
multivariate and temporal indexing, such that `spatial_domain` is then
ignored, and the model collapses to a generalized additive model using
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html). To specify a univariate
spatial model, the user must specify `spatial_domain` and either
`space_term=""` or `spacetime_term=""`, where the latter two are then
parsed to include a single exogenous variance for the single variable

|                                                                                                                                     |                                                                                                                                                                                                                                    |
|-------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Model type**                                                                                                                      | **How to specify**                                                                                                                                                                                                                 |
| Generalized additive model                                                                                                          | specify `spatial_domain=NULL` `space_term=""` and `spacetime_term=""`, and then use `formula` to specify splines and covariates                                                                                                    |
| Dynamic structural equation model (including vector autoregressive, dynamic factor analysis, ARIMA, and structural equation models) | specify `spatial_domain=NULL` and use `spacetime_term` to specify interactions among variables and over time                                                                                                                       |
| Univariate spatio-temporal model, or multiple independence spatio-temporal variables                                                | specify `spatial_domain` and `spacetime_term=""`, where the latter is then parsed to include a single exogenous variance for the single variable                                                                                   |
| Multivariate spatial model including interactions                                                                                   | specify `spatial_domain` and use `space_term` to specify spatial interactions                                                                                                                                                      |
| Vector autoregressive spatio-temporal model (i.e., lag-1 interactions among variables)                                              | specify `spatial_domain` and use `spacetime_term=""` to specify interactions among variables and over time, where spatio-temporal variables are constructed via the separable interaction of `spacetime_term` and `spatial_domain` |

**Model building notes**

- `binomial familes`: A binomial family can be specified in only one
  way: the response is the observed proportion (proportion = successes /
  trials), and the 'weights' argument is used to specify the Binomial
  size (trials, N) parameter (`proportion ~ ..., weights = N`).

- `factor models`: If a factor model is desired, the factor(s) must be
  named and included in the `variables`. The factor is then modeled for
  `space_term`, `time_term`, and `spacetime_term` and it's variance must
  be fixed a priori for any term where it is not being used.

## See also

Details section of
[`make_dsem_ram()`](https://vast-lib.github.io/tinyVAST/dev/reference/make_dsem_ram.md)
for a summary of the math involved with constructing the DSEM, and
[doi:10.1111/2041-210X.14289](https://doi.org/10.1111/2041-210X.14289)
for more background on math and inference

[doi:10.48550/arXiv.2401.10193](https://doi.org/10.48550/arXiv.2401.10193)
for more details on how GAM, SEM, and DSEM components are combined from
a statistical and software-user perspective

[`summary.tinyVAST()`](https://vast-lib.github.io/tinyVAST/dev/reference/summary.tinyVAST.md)
to visualize parameter estimates related to SEM and DSEM model
components

## Examples

``` r
# Simulate a seperable two-dimensional AR1 spatial process
n_x = n_y = 25
n_w = 10
R_xx = exp(-0.4 * abs(outer(1:n_x, 1:n_x, FUN="-")) )
R_yy = exp(-0.4 * abs(outer(1:n_y, 1:n_y, FUN="-")) )
z = mvtnorm::rmvnorm(1, sigma=kronecker(R_xx,R_yy) )

# Simulate nuissance parameter z from oscillatory (day-night) process
w = sample(1:n_w, replace=TRUE, size=length(z))
Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), w=w, z=as.vector(z) + cos(w/n_w*2*pi))
Data$n = Data$z + rnorm(nrow(Data), sd=1)

# Add columns for multivariate and/or temporal dimensions
Data$var = "n"

# make SPDE mesh for spatial term
mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], n=100 )

# fit model with cyclic confounder as GAM term
out = tinyVAST( data = Data,
                formula = n ~ s(w),
                spatial_domain = mesh,
                space_term = "n <-> n, sd_n" )

# Run crossvalidation (too slow for CRAN)
# \donttest{
CV = cv::cv( out, k = 4 )
#> R RNG seed set to 299103
#> Error in eval(call, parent.frame()): object 'mesh' not found
CV
#> Error: object 'CV' not found
# }
```
