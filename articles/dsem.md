# Dynamic structural equation models

``` r

library(tinyVAST)
set.seed(101)
options("tinyVAST.verbose" = FALSE)
```

`tinyVAST` includes features to fit a dynamic structural equation model
(Thorson et al. 2024). We here show this using a bivariate vector
autoregressive model for wolf and moose abundance on Isle Royale.

``` r

data(isle_royale, package="dsem")

# Convert to long-form
data = expand.grid( "time"=isle_royale[,1], "var"=colnames(isle_royale[,2:3]) )
data$logn = unlist(log(isle_royale[2:3]))

# Define cross-lagged DSEM
dsem = "
  # Link, lag, param_name
  wolves -> wolves, 1, arW
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
  moose -> moose, 1, arM
  #wolves -> moose, 0, corr
  wolves <-> moose, 0, corr
"

# fit model
mytiny = tinyVAST( spacetime_term = dsem,
                 data = data,
                 times = isle_royale[,1],
                 variables = colnames(isle_royale[,2:3]),
                 formula = logn ~ 0 + var )
mytiny
#> Call: 
#> tinyVAST(formula = logn ~ 0 + var, data = data, spacetime_term = dsem, 
#>     times = isle_royale[, 1], variables = colnames(isle_royale[, 
#>         2:3]))
#> 
#> Run time: 
#> Time difference of 0.6704443 secs
#> 
#> Family: 
#> $obs
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> 
#> 
#> 
#> sdreport(.) result
#>               Estimate   Std. Error
#> alpha_j     3.32526214 2.483494e-01
#> alpha_j     6.44165368 2.116033e-01
#> beta_z      0.89304279 8.420631e-02
#> beta_z      0.01420916 1.279150e-01
#> beta_z     -0.13239998 3.454996e-02
#> beta_z      0.86147624 7.107380e-02
#> beta_z     -0.01285882 5.063467e-02
#> beta_z      0.37727137 3.504315e-02
#> beta_z      0.17062265 1.584742e-02
#> log_sigma -12.50889085 1.839990e+04
#> Maximum gradient component: 4.38807e-05 
#> 
#> Proportion conditional deviance explained: 
#> [1] 1
#> 
#> spacetime_term: 
#>   heads     to   from parameter start lag    Estimate  Std_Error    z_value
#> 1     1 wolves wolves         1  <NA>   1  0.89304279 0.08420631 10.6054140
#> 2     1 wolves  moose         2  <NA>   1  0.01420916 0.12791499  0.1110828
#> 3     1  moose wolves         3  <NA>   1 -0.13239998 0.03454996 -3.8321309
#> 4     1  moose  moose         4  <NA>   1  0.86147624 0.07107380 12.1208692
#> 5     2  moose wolves         5  <NA>   0 -0.01285882 0.05063467 -0.2539528
#> 6     2 wolves wolves         6  <NA>   0  0.37727137 0.03504315 10.7659078
#> 7     2  moose  moose         7  <NA>   0  0.17062265 0.01584742 10.7665883
#>        p_value
#> 1 2.812157e-26
#> 2 9.115507e-01
#> 3 1.270381e-04
#> 4 8.188607e-34
#> 5 7.995320e-01
#> 6 4.986766e-27
#> 7 4.950052e-27
#> 
#> Fixed terms: 
#>           Estimate Std_Error  z_value       p_value
#> varwolves 3.325262 0.2483494 13.38945  6.969522e-41
#> varmoose  6.441654 0.2116033 30.44212 1.522921e-203
#> 
#> Sanity check:

# Deviance explained relative to both intercepts
# Note that a process-error-only estimate with have devexpl -> 1
deviance_explained( mytiny, 
                    null_formula = logn ~ 0 + var )
#> [1] 1

# See summary
knitr::kable( summary(mytiny,"spacetime_term"), digits=3 )
```

| heads | to     | from   | parameter | start | lag | Estimate | Std_Error | z_value | p_value |
|:------|:-------|:-------|:----------|:------|:----|---------:|----------:|--------:|--------:|
| 1     | wolves | wolves | 1         | NA    | 1   |    0.893 |     0.084 |  10.605 |   0.000 |
| 1     | wolves | moose  | 2         | NA    | 1   |    0.014 |     0.128 |   0.111 |   0.912 |
| 1     | moose  | wolves | 3         | NA    | 1   |   -0.132 |     0.035 |  -3.832 |   0.000 |
| 1     | moose  | moose  | 4         | NA    | 1   |    0.861 |     0.071 |  12.121 |   0.000 |
| 2     | moose  | wolves | 5         | NA    | 0   |   -0.013 |     0.051 |  -0.254 |   0.800 |
| 2     | wolves | wolves | 6         | NA    | 0   |    0.377 |     0.035 |  10.766 |   0.000 |
| 2     | moose  | moose  | 7         | NA    | 0   |    0.171 |     0.016 |  10.767 |   0.000 |

And we can specifically inspect the estimated interaction matrix:

|        | wolves |  moose |
|:-------|-------:|-------:|
| wolves |  0.893 | -0.132 |
| moose  |  0.014 |  0.861 |

We can then compare this with package `dsem`

``` r

library(dsem)

# Keep in wide-form
dsem_data = ts( log(isle_royale[,2:3]), start=1959)

# fit without delta0
# Fitting with family = "fixed" seems stable on CRAN checks
mydsem = dsem::dsem( 
  sem = dsem,
  tsdata = dsem_data,
  control = dsem_control(
    getsd = FALSE
  ), 
)
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      7  Fixed
#> 2             mu_j                      2 Random
mydsem
#> $par
#>       beta_z       beta_z       beta_z       beta_z       beta_z       beta_z 
#>  0.895834720  0.007358847 -0.124879928  0.874884847 -0.014394603  0.378522244 
#>       beta_z 
#>  0.172997994 
#> 
#> $objective
#> [1] 7.739638
#> attr(,"logarithm")
#> [1] TRUE
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 41
#> 
#> $evaluations
#> function gradient 
#>       56       42 
#> 
#> $message
#> [1] "relative convergence (4)"

# See summary
knitr::kable( summary(mydsem), digits=3 )
```

| path | lag | name | start | parameter | first | second | direction | Estimate |
|:---|---:|:---|---:|---:|:---|:---|:---|---:|
| wolves -\> wolves | 1 | arW | NA | 1 | wolves | wolves | 1 | 0.896 |
| moose -\> wolves | 1 | MtoW | NA | 2 | moose | wolves | 1 | 0.007 |
| wolves -\> moose | 1 | WtoM | NA | 3 | wolves | moose | 1 | -0.125 |
| moose -\> moose | 1 | arM | NA | 4 | moose | moose | 1 | 0.875 |
| wolves \<-\> moose | 0 | corr | NA | 5 | wolves | moose | 2 | -0.014 |
| wolves \<-\> wolves | 0 | V\[wolves\] | NA | 6 | wolves | wolves | 2 | 0.379 |
| moose \<-\> moose | 0 | V\[moose\] | NA | 7 | moose | moose | 2 | 0.172 |

where we again inspect the estimated interaction matrix:

|        | wolves |  moose |
|:-------|-------:|-------:|
| wolves |  0.896 | -0.125 |
| moose  |  0.007 |  0.875 |

Runtime for this vignette: 2.37 secs

## Works cited

Thorson, James T., Alexander G. Andrews III, Timothy E. Essington, and
Scott I. Large. 2024. “Dynamic Structural Equation Models Synthesize
Ecosystem Dynamics Constrained by Ecological Mechanisms.” *Methods in
Ecology and Evolution* 15 (4): 744–55.
<https://doi.org/10.1111/2041-210X.14289>.
