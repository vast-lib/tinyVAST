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

True loadings

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
#> Time difference of 2.192787 secs
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
#> alpha_j    0.07569060 0.31851467
#> alpha_j   -0.02016919 0.39764587
#> alpha_j    0.22317503 0.21847588
#> alpha_j    0.14727256 0.27058070
#> alpha_j   -0.26515247 0.14638935
#> theta_z    0.68017152 0.11511151
#> theta_z    0.68287021 0.15773931
#> theta_z    0.31702161 0.10358374
#> theta_z    0.52124050 0.10914511
#> theta_z    0.14820460 0.09200849
#> theta_z    0.51875214 0.13709851
#> theta_z   -0.31999054 0.10049337
#> theta_z    0.23600834 0.10641093
#> theta_z   -0.21614664 0.09653182
#> log_sigma -0.52203929 0.06761656
#> log_sigma  0.21850709 0.13313013
#> log_kappa -0.26762174 0.21031509
#> Maximum gradient component: 0.002544271 
#> 
#> Proportion conditional deviance explained: 
#> [1] 0.5312474
#> 
#> space_term: 
#>    heads to from parameter start   Estimate  Std_Error   z_value      p_value
#> 1      1  1   f1         1  <NA>  0.6801715 0.11511151  5.908805 3.445975e-09
#> 2      1  2   f1         2  <NA>  0.6828702 0.15773931  4.329106 1.497157e-05
#> 3      1  3   f1         3  <NA>  0.3170216 0.10358374  3.060535 2.209422e-03
#> 4      1  4   f1         4  <NA>  0.5212405 0.10914511  4.775665 1.791141e-06
#> 5      1  5   f1         5  <NA>  0.1482046 0.09200849  1.610771 1.072297e-01
#> 6      1  2   f2         6  <NA>  0.5187521 0.13709851  3.783791 1.544575e-04
#> 7      1  3   f2         7  <NA> -0.3199905 0.10049337 -3.184196 1.451569e-03
#> 8      1  4   f2         8  <NA>  0.2360083 0.10641093  2.217896 2.656195e-02
#> 9      1  5   f2         9  <NA> -0.2161466 0.09653182 -2.239123 2.514790e-02
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
#> factor(species)1  0.07569060 0.3185147  0.23763614 0.81216332
#> factor(species)2 -0.02016919 0.3976459 -0.05072149 0.95954745
#> factor(species)3  0.22317503 0.2184759  1.02150879 0.30701345
#> factor(species)4  0.14727256 0.2705807  0.54428332 0.58624652
#> factor(species)5 -0.26515247 0.1463893 -1.81128250 0.07009713
#> 
#> Sanity check: 
#> 
#> **Possible issues detected! Check output of sanity().**
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

Rotated true loadings

with the estimated loadings

``` r
# Extract and rotate estimated loadings
Lhat_cf = matrix( 0, nrow=n_c, ncol=2 )
Lhat_cf[lower.tri(Lhat_cf,diag=TRUE)] = as.list(out$sdrep, what="Estimate")$theta_z
Lhat_cf = rotate_pca( L_tf=Lhat_cf, order="decreasing" )$L_tf
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

Rotated estimated loadings

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
#> Warning in sqrt(Eigen$values): NaNs produced
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

Rotated estimated loadings with full rank

Runtime for this vignette: 9.28 secs

## Works cited

Thorson, Mark D. Scheuerell, Andrew O. Shelton, Kevin E. See, Hans J.
Skaug, and Kasper Kristensen. 2015. “Spatial Factor Analysis: A New Tool
for Estimating Joint Species Distributions and Correlations in Species
Range.” *Methods in Ecology and Evolution* 6 (6): 627–37.
<https://doi.org/10.1111/2041-210X.12359>.
