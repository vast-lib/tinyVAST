# Spatial factor analysis

``` r

library(tinyVAST)
library(fmesher)
set.seed(101)
options("tinyVAST.verbose" = FALSE)
```

`tinyVAST` is an R package for fitting vector autoregressive
spatio-temporal (VAST) models. We here explore the capacity to specify a
spatial factor analysis, where the spatial pattern for multiple
variables is described via their estimated association with a small
number of spatial latent variables (Thorson et al. 2015).

## Spatial factor analysis

We first explore the ability to specify two latent variables for five
manifest variables. To start we simulate two spatial latent variables,
project via a simulated loadings matrix, and then simulate a Tweedie
response for each manifest variable:

``` r

# Simulate settings
theta_xy = 0.4
n_x = n_y = 10
n_c = 5
rho = 0.8
resid_sd = 0.5

# Simulate GMRFs
R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
R_ss = kronecker(X=R_s, Y=R_s)
delta_fs = mvtnorm::rmvnorm(n_c, sigma=R_ss )

#
L_cf = matrix( rnorm(n_c^2), nrow=n_c )
L_cf[,3:5] = 0
L_cf = L_cf + resid_sd * diag(n_c)

#
d_cs = L_cf %*% delta_fs
```

Where we can inspect the simulated loadings matrix

``` r

dimnames(L_cf) = list( paste0("Var ", 1:nrow(L_cf)),
                       paste0("Factor ", 1:ncol(L_cf)) )
knitr::kable( L_cf,
              digits=2, caption="True loadings")
```

|       | Factor 1 | Factor 2 | Factor 3 | Factor 4 | Factor 5 |
|:------|---------:|---------:|---------:|---------:|---------:|
| Var 1 |    -0.77 |     0.26 |      0.0 |      0.0 |      0.0 |
| Var 2 |    -0.66 |     1.03 |      0.0 |      0.0 |      0.0 |
| Var 3 |    -0.30 |    -0.54 |      0.5 |      0.0 |      0.0 |
| Var 4 |    -0.57 |     0.34 |      0.0 |      0.5 |      0.0 |
| Var 5 |    -0.03 |    -0.24 |      0.0 |      0.0 |      0.5 |

True loadings {.table}

We then specify the model as expected by *tinyVAST*:

``` r

# Shape into longform data-frame and add error
Data = data.frame( expand.grid(species=1:n_c, x=1:n_x, y=1:n_y),
                   "var"="logn", "z"=exp(as.vector(d_cs)) )
Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$z, phi=0.5, power=1.5 )
mean(Data$n==0)
#> [1] 0.03

# make mesh
mesh = fm_mesh_2d( Data[,c('x','y')] )

#
sem = "
  f1 -> 1, l1
  f1 -> 2, l2
  f1 -> 3, l3
  f1 -> 4, l4
  f1 -> 5, l5
  f2 -> 2, l6
  f2 -> 3, l7
  f2 -> 4, l8
  f2 -> 5, l9
  f1 <-> f1, NA, 1
  f2 <-> f2, NA, 1
  1 <-> 1, NA, 0
  2 <-> 2, NA, 0
  3 <-> 3, NA, 0
  4 <-> 4, NA, 0
  5 <-> 5, NA, 0
"

# fit model
out = tinyVAST( space_term = sem,
           data = Data,
           formula = n ~ 0 + factor(species),
           spatial_domain = mesh,
           family = tweedie(),
           variables = c( "f1", "f2", 1:n_c ),
           space_columns = c("x","y"),
           variable_column = "species",
           time_column = "time",
           distribution_column = "dist",
           control = tinyVASTcontrol(gmrf="proj") )
out
#> Call: 
#> tinyVAST(formula = n ~ 0 + factor(species), data = Data, space_term = sem, 
#>     family = tweedie(), spatial_domain = mesh, control = tinyVASTcontrol(gmrf = "proj"), 
#>     space_columns = c("x", "y"), time_column = "time", variable_column = "species", 
#>     variables = c("f1", "f2", 1:n_c), distribution_column = "dist")
#> 
#> Run time: 
#> Time difference of 1.509853 secs
#> 
#> Family: 
#> $obs
#> 
#> Family: tweedie 
#> Link function: log 
#> 
#> 
#> 
#> 
#> sdreport(.) result
#>              Estimate Std. Error
#> alpha_j    0.07570783 0.31850774
#> alpha_j   -0.02014888 0.39763412
#> alpha_j    0.22318277 0.21847161
#> alpha_j    0.14728087 0.27057852
#> alpha_j   -0.26514652 0.14638317
#> theta_z    0.68016356 0.11510781
#> theta_z    0.68285927 0.15773444
#> theta_z    0.31701846 0.10358197
#> theta_z    0.52123769 0.10914344
#> theta_z    0.14819781 0.09200614
#> theta_z    0.51873790 0.13709521
#> theta_z   -0.31998633 0.10049164
#> theta_z    0.23601680 0.10640963
#> theta_z   -0.21613586 0.09652980
#> log_sigma -0.52205331 0.06761637
#> log_sigma  0.21851154 0.13313019
#> log_kappa -0.26761196 0.21030826
#> Maximum gradient component: 0.002233904 
#> 
#> Proportion conditional deviance explained: 
#> [1] 0.5312466
#> 
#> space_term: 
#>    heads to from parameter start   Estimate  Std_Error   z_value      p_value
#> 1      1  1   f1         1  <NA>  0.6801636 0.11510781  5.908926 3.443444e-09
#> 2      1  2   f1         2  <NA>  0.6828593 0.15773444  4.329170 1.496721e-05
#> 3      1  3   f1         3  <NA>  0.3170185 0.10358197  3.060556 2.209261e-03
#> 4      1  4   f1         4  <NA>  0.5212377 0.10914344  4.775713 1.790719e-06
#> 5      1  5   f1         5  <NA>  0.1481978 0.09200614  1.610738 1.072368e-01
#> 6      1  2   f2         6  <NA>  0.5187379 0.13709521  3.783778 1.544654e-04
#> 7      1  3   f2         7  <NA> -0.3199863 0.10049164 -3.184209 1.451504e-03
#> 8      1  4   f2         8  <NA>  0.2360168 0.10640963  2.218002 2.655468e-02
#> 9      1  5   f2         9  <NA> -0.2161359 0.09652980 -2.239058 2.515211e-02
#> 10     2 f1   f1         0     1  1.0000000         NA        NA           NA
#> 11     2 f2   f2         0     1  1.0000000         NA        NA           NA
#> 12     2  1    1         0     0  0.0000000         NA        NA           NA
#> 13     2  2    2         0     0  0.0000000         NA        NA           NA
#> 14     2  3    3         0     0  0.0000000         NA        NA           NA
#> 15     2  4    4         0     0  0.0000000         NA        NA           NA
#> 16     2  5    5         0     0  0.0000000         NA        NA           NA
#> 
#> Fixed terms: 
#>                     Estimate Std_Error     z_value    p_value
#> factor(species)1  0.07570783 0.3185077  0.23769543 0.81211733
#> factor(species)2 -0.02014888 0.3976341 -0.05067191 0.95958696
#> factor(species)3  0.22318277 0.2184716  1.02156419 0.30698721
#> factor(species)4  0.14728087 0.2705785  0.54431841 0.58622238
#> factor(species)5 -0.26514652 0.1463832 -1.81131833 0.07009159
```

We can compare the true loadings (rotated to optimize comparison):

``` r

Lrot_cf = rotate_pca( L_cf )$L_tf
dimnames(Lrot_cf) = list( paste0("Var ", 1:nrow(Lrot_cf)),
                       paste0("Factor ", 1:ncol(Lrot_cf)) )
knitr::kable( Lrot_cf,
              digits=2, caption="Rotated true loadings")
```

|       | Factor 1 | Factor 2 | Factor 3 | Factor 4 | Factor 5 |
|:------|---------:|---------:|---------:|---------:|---------:|
| Var 1 |    -0.71 |     0.34 |     0.00 |    -0.11 |    -0.19 |
| Var 2 |    -1.20 |    -0.18 |     0.00 |    -0.18 |     0.12 |
| Var 3 |     0.22 |     0.73 |    -0.17 |    -0.09 |     0.10 |
| Var 4 |    -0.70 |     0.25 |     0.06 |     0.37 |     0.03 |
| Var 5 |     0.17 |     0.23 |     0.47 |    -0.08 |     0.03 |

Rotated true loadings {.table}

with the estimated loadings

``` r

# Extract and rotate estimated loadings
Lhat_cf = matrix( 0, nrow=n_c, ncol=2 )
Lhat_cf[lower.tri(Lhat_cf,diag=TRUE)] = as.list(out$sdrep, what="Estimate")$theta_z
Lhat_cf = rotate_pca( L_tf=Lhat_cf, order="decreasing" )$L_tf
#> Warning in sqrt(Eigen$values): NaNs produced
```

Where we can compared the estimated and true loadings matrices:

``` r

dimnames(Lhat_cf) = list( paste0("Var ", 1:nrow(Lhat_cf)),
                       paste0("Factor ", 1:ncol(Lhat_cf)) )
knitr::kable( Lhat_cf,
              digits=2, caption="Rotated estimated loadings" )
```

|       | Factor 1 | Factor 2 |
|:------|---------:|---------:|
| Var 1 |     0.64 |    -0.23 |
| Var 2 |     0.82 |     0.26 |
| Var 3 |     0.19 |    -0.41 |
| Var 4 |     0.57 |     0.05 |
| Var 5 |     0.07 |    -0.25 |

Rotated estimated loadings {.table}

Or we can specify the model while ensuring that residual spatial
variation is also captured:

``` r

#
sem = "
  f1 -> 1, l1
  f1 -> 2, l2
  f1 -> 3, l3
  f1 -> 4, l4
  f1 -> 5, l5
  f2 -> 2, l6
  f2 -> 3, l7
  f2 -> 4, l8
  f2 -> 5, l9
  f1 <-> f1, NA, 1
  f2 <-> f2, NA, 1
  1 <-> 1, sd_resid
  2 <-> 2, sd_resid
  3 <-> 3, sd_resid
  4 <-> 4, sd_resid
  5 <-> 5, sd_resid
"

# fit model
out = tinyVAST( space_term = sem,
           data = Data,
           formula = n ~ 0 + factor(species),
           spatial_domain = mesh,
           family = list( "obs"=tweedie() ),
           variables = c( "f1", "f2", 1:n_c ),
           space_columns = c("x","y"),
           variable_column = "species",
           time_column = "time",
           distribution_column = "dist",
           control = tinyVASTcontrol(gmrf="proj") )

# Extract and rotate estimated loadings
Lhat_cf = matrix( 0, nrow=n_c, ncol=2 )
Lhat_cf[lower.tri(Lhat_cf,diag=TRUE)] = as.list(out$sdrep, what="Estimate")$theta_z
#> Warning in Lhat_cf[lower.tri(Lhat_cf, diag = TRUE)] = as.list(out$sdrep, :
#> number of items to replace is not a multiple of replacement length
Lhat_cf = rotate_pca( L_tf=Lhat_cf, order="decreasing" )$L_tf
```

Where we can again compared the estimated and true loadings matrices:

``` r

dimnames(Lhat_cf) = list( paste0("Var ", 1:nrow(Lhat_cf)),
                       paste0("Factor ", 1:ncol(Lhat_cf)) )
knitr::kable( Lhat_cf,
              digits=2, caption="Rotated estimated loadings with full rank" )
```

|       | Factor 1 | Factor 2 |
|:------|---------:|---------:|
| Var 1 |     0.69 |    -0.18 |
| Var 2 |     0.74 |     0.23 |
| Var 3 |     0.07 |    -0.42 |
| Var 4 |     0.47 |    -0.02 |
| Var 5 |     0.07 |    -0.07 |

Rotated estimated loadings with full rank {.table}

Runtime for this vignette: 7.91 secs

## Works cited

Thorson, Mark D. Scheuerell, Andrew O. Shelton, Kevin E. See, Hans J.
Skaug, and Kasper Kristensen. 2015. “Spatial Factor Analysis: A New Tool
for Estimating Joint Species Distributions and Correlations in Species
Range.” *Methods in Ecology and Evolution* 6 (6): 627–37.
<https://doi.org/10.1111/2041-210X.12359>.
