# Extract covariance

Extract the covariance resulting from a specified path structure and
estimated parameters for a SEM or DSEM term in tinyVAST

## Usage

``` r
term_covariance(
  object,
  what = c("space_term", "time_term", "spacetime_term"),
  pred = c("one", "two"),
  n_times = NULL
)
```

## Arguments

- object:

  Output from
  [`tinyVAST`](https://vast-lib.github.io/tinyVAST/reference/tinyVAST.md)

- what:

  Which SEM or DSEM term to extract

- pred:

  Extract the term `what` for which linear predictor

- n_times:

  The number of times to include when calculating covariance for a DSEM
  component, i.e., `time_term` or `spacetime_term`. If missing, the
  default is to use the one more than the maximum specified lag (e.g.,
  `n_times=2` by default when the maximum `lag=1`)

## Value

The covariance matrix among variables

## Details

tinyVAST constructs the covariance from specified path structure and
estimated parameters

## Examples

``` r
# Extract covariance for spatial factor analysis (too slow for CRAN)
# \donttest{
# Simulate settings
set.seed(101)
theta_xy = 0.4
n_x = n_y = 10
n_c = 3           # Number of species
n_f = 1           # Number of factors
rho = 0.8
resid_sd = 0.5

# Simulate GMRFs
R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
R_ss = kronecker(X=R_s, Y=R_s)
delta_fs = mvtnorm::rmvnorm(n_c, sigma=R_ss )

# Simulate loadings for two factors
L_cf = matrix( rnorm(n_c^2), nrow=n_c )
L_cf[,seq(from=n_f+1, to=n_c)] = 0
L_cf = L_cf + resid_sd * diag(n_c)

# Simulate correlated densities
d_cs = L_cf %*% delta_fs

# Shape into longform data-frame and add error
Data = data.frame( expand.grid(species=1:n_c, x=1:n_x, y=1:n_y),
                   "var"="logn", "z"=exp(as.vector(d_cs)) )
Data$n = rnorm( n=nrow(Data), mean=Data$z, sd=1 )

# make mesh
mesh = fmesher::fm_mesh_2d( Data[,c('x','y')] )

# Specify factor model with two factors and additional independent variance with shared SD
sem = "
  # Loadings matrix
  f1 -> 1, l1
  f1 -> 2, l2
  f1 -> 3, l3

  # Factor variance = 1
  f1 <-> f1, NA, 1

  # Shared residual variance
  1 <-> 1, sd, 1
  2 <-> 2, sd, 1
  3 <-> 3, sd, 1
"

# fit model
out = tinyVAST( space_term = sem,
           data = Data,
           formula = n ~ 0 + factor(species),
           spatial_domain = mesh,
           variables = c( "f1", 1:n_c ),
           space_columns = c("x","y"),
           variable_column = "species",
           time_column = "time",
           distribution_column = "dist" )

# Extract covariance among species and factors, where
# estimated covariance is obtained by ignoring factors
V = term_covariance( out, what = "space_term", pred = "one" )
# }
```
