# Make a RAM (Reticular Action Model)

`make_dsem_ram` converts SEM arrow notation to `ram` describing SEM
parameters

## Usage

``` r
make_dsem_ram(
  dsem,
  times,
  variables,
  covs = NULL,
  quiet = FALSE,
  remove_na = TRUE
)
```

## Arguments

- dsem:

  dynamic structural equation model structure, passed to either
  [`specifyModel`](https://rdrr.io/pkg/sem/man/specifyModel.html) or
  [`specifyEquations`](https://rdrr.io/pkg/sem/man/specifyModel.html)
  and then parsed to control the set of path coefficients and
  variance-covariance parameters

- times:

  A character vector listing the set of times in order

- variables:

  A character vector listing the set of variables

- covs:

  optional: a character vector of one or more elements, with each
  element giving a string of variable names, separated by commas.
  Variances and covariances among all variables in each such string are
  added to the model. For confirmatory factor analysis models specified
  via `cfa`, `covs` defaults to all of the factors in the model, thus
  specifying all variances and covariances among these factors.
  *Warning*: `covs="x1, x2"` and `covs=c("x1", "x2")` are *not*
  equivalent: `covs="x1, x2"` specifies the variance of `x1`, the
  variance of `x2`, *and* their covariance, while `covs=c("x1", "x2")`
  specifies the variance of `x1` and the variance of `x2` *but not*
  their covariance.

- quiet:

  Boolean indicating whether to print messages to terminal

- remove_na:

  Boolean indicating whether to remove NA values from RAM (default) or
  not. `remove_NA=FALSE` might be useful for exploration and diagnostics
  for advanced users

## Value

A reticular action module (RAM) describing dependencies

## Details

**RAM specification using arrow-and-lag notation**

Each line of the RAM specification for `make_dsem_ram` consists of four
(unquoted) entries, separated by commas:

- 1\. Arrow specification::

  This is a simple formula, of the form `A -> B` or, equivalently,
  `B <- A` for a regression coefficient (i.e., a single-headed or
  directional arrow); `A <-> A` for a variance or `A <-> B` for a
  covariance (i.e., a double-headed or bidirectional arrow). Here, `A`
  and `B` are variable names in the model. If a name does not correspond
  to an observed variable, then it is assumed to be a latent variable.
  Spaces can appear freely in an arrow specification, and there can be
  any number of hyphens in the arrows, including zero: Thus, e.g.,
  `A->B`, `A --> B`, and `A>B` are all legitimate and equivalent.

- 2\. Lag (using positive values)::

  An integer specifying whether the linkage is simultaneous (`lag=0`) or
  lagged (e.g., `X -> Y, 1, XtoY` indicates that X in time T affects Y
  in time T+1), where only one-headed arrows can be lagged. Using
  positive values to indicate lags then matches the notational
  convention used in package dynlm.

- 3\. Parameter name::

  The name of the regression coefficient, variance, or covariance
  specified by the arrow. Assigning the same name to two or more arrows
  results in an equality constraint. Specifying the parameter name as
  `NA` produces a fixed parameter.

- 4\. Value::

  start value for a free parameter or value of a fixed parameter. If
  given as `NA` (or simply omitted), the model is provide a default
  starting value.

Lines may end in a comment following \#. The function extends code
copied from package `sem` under licence GPL (\>= 2) with permission from
John Fox.

**Simultaneous autoregressive process for simultaneous and lagged
effects**

This text then specifies linkages in a multivariate time-series model
for variables \\\mathbf X\\ with dimensions \\T \times C\\ for \\T\\
times and \\C\\ variables. `make_dsem_ram` then parses this text to
build a path matrix \\\mathbf P\\ with dimensions \\TC \times TC\\,
where \\\rho\_{k_2,k_1}\\ represents the impact of \\x\_{t_1,c_1}\\ on
\\x\_{t_2,c_2}\\, where \\k_1=T c_1+t_1\\ and \\k_2=T c_2+t_2\\. This
path matrix defines a simultaneous equation

\$\$ \mathrm{vec}(\mathbf X) = \mathbf P \mathrm{vec}(\mathbf X) +
\mathrm{vec}(\mathbf \Delta)\$\$

where \\\mathbf \Delta\\ is a matrix of exogenous errors with covariance
\\\mathbf{V = \Gamma \Gamma}^t\\, where \\\mathbf \Gamma\\ is the
Cholesky of exogenous covariance. This simultaneous autoregressive (SAR)
process then results in \\\mathbf X\\ having covariance:

\$\$ \mathrm{Cov}(\mathbf X) = \mathbf{(I - P)}^{-1} \mathbf{\Gamma
\Gamma}^t \mathbf{((I - P)}^{-1})^t \$\$

Usefully, it is also easy to compute the inverse-covariance (precision)
matrix \\\mathbf{Q = V}^{-1}\\:

\$\$ \mathbf{Q} = (\mathbf{\Gamma}^{-1} \mathbf{(I - P)})^t
\mathbf{\Gamma}^{-1} \mathbf{(I - P)} \$\$

**Example: univariate and first-order autoregressive model**

This simultaneous autoregressive (SAR) process across variables and
times allows the user to specify both simultaneous effects (effects
among variables within year \\T\\) and lagged effects (effects among
variables among years \\T\\). As one example, consider a univariate and
first-order autoregressive process where \\T=4\\. with independent
errors. This is specified by passing
` dsem = X -> X, 1, rho; X <-> X, 0, sigma ` to `make_dsem_ram`. This is
then parsed to a RAM:

|           |        |          |                |           |
|-----------|--------|----------|----------------|-----------|
| **heads** | **to** | **from** | **paarameter** | **start** |
| 1         | 2      | 1        | 1              | NA        |
| 1         | 3      | 2        | 1              | NA        |
| 1         | 4      | 3        | 1              | NA        |
| 2         | 1      | 1        | 2              | NA        |
| 2         | 2      | 2        | 2              | NA        |
| 2         | 3      | 3        | 2              | NA        |
| 2         | 4      | 4        | 2              | NA        |

Rows of this RAM where `heads=1` are then interpreted to construct the
path matrix \\\mathbf P\\:

    \deqn{ \mathbf P = \begin{bmatrix}
        0 & 0 & 0 & 0 \
        \rho & 0 & 0 & 0 \
        0 & \rho & 0 & 0 \
        0 & 0 & \rho & 0\
        \end{bmatrix} }

While rows where `heads=2` are interpreted to construct the Cholesky of
exogenous covariance \\\mathbf \Gamma\\:

    \deqn{ \mathbf \Gamma = \begin{bmatrix}
        \sigma & 0 & 0 & 0 \
        0 & \sigma & 0 & 0 \
        0 & 0 & \sigma & 0 \
        0 & 0 & 0 & \sigma\
        \end{bmatrix} }

with two estimated parameters \\\mathbf \beta = (\rho, \sigma) \\. This
then results in covariance:

    \deqn{ \mathrm{Cov}(\mathbf X) = \sigma^2 \begin{bmatrix}
        1 & \rho^1 & \rho^2 & \rho^3 \
        \rho^1 & 1 & \rho^1 & \rho^2 \
        \rho^2 & \rho^1 & 1 & \rho^1 \
        \rho^3 & \rho^2 & \rho^1 & 1\
        \end{bmatrix} }

Similarly, the arrow-and-lag notation can be used to specify a SAR
representing a conventional structural equation model (SEM),
cross-lagged (a.k.a. vector autoregressive) models (VAR), dynamic factor
analysis (DFA), or many other time-series models.

## Examples

``` r
# Univariate AR1
dsem = "
  X -> X, 1, rho
  X <-> X, 0, sigma
"
make_dsem_ram( dsem=dsem, variables="X", times=1:4 )
#> $model
#>      path      lag name    start parameter first second direction
#> [1,] "X->X"    "1" "rho"   NA    "1"       "X"   "X"    "1"      
#> [2,] "X <-> X" "0" "sigma" NA    "2"       "X"   "X"    "2"      
#> 
#> $ram
#>   heads to from parameter start
#> 1     1  2    1         1  <NA>
#> 2     1  3    2         1  <NA>
#> 3     1  4    3         1  <NA>
#> 5     2  1    1         2  <NA>
#> 6     2  2    2         2  <NA>
#> 7     2  3    3         2  <NA>
#> 8     2  4    4         2  <NA>
#> 
#> attr(,"class")
#> [1] "dsem_ram"

# Univariate AR2
dsem = "
  X -> X, 1, rho1
  X -> X, 2, rho2
  X <-> X, 0, sigma
"
make_dsem_ram( dsem=dsem, variables="X", times=1:4 )
#> $model
#>      path      lag name    start parameter first second direction
#> [1,] "X->X"    "1" "rho1"  NA    "1"       "X"   "X"    "1"      
#> [2,] "X->X"    "2" "rho2"  NA    "2"       "X"   "X"    "1"      
#> [3,] "X <-> X" "0" "sigma" NA    "3"       "X"   "X"    "2"      
#> 
#> $ram
#>    heads to from parameter start
#> 1      1  2    1         1  <NA>
#> 2      1  3    2         1  <NA>
#> 3      1  4    3         1  <NA>
#> 5      1  3    1         2  <NA>
#> 6      1  4    2         2  <NA>
#> 9      2  1    1         3  <NA>
#> 10     2  2    2         3  <NA>
#> 11     2  3    3         3  <NA>
#> 12     2  4    4         3  <NA>
#> 
#> attr(,"class")
#> [1] "dsem_ram"

# Bivariate VAR
dsem = "
  X -> X, 1, XtoX
  X -> Y, 1, XtoY
  Y -> X, 1, YtoX
  Y -> Y, 1, YtoY
  X <-> X, 0, sdX
  Y <-> Y, 0, sdY
"
make_dsem_ram( dsem=dsem, variables=c("X","Y"), times=1:4 )
#> $model
#>      path      lag name   start parameter first second direction
#> [1,] "X->X"    "1" "XtoX" NA    "1"       "X"   "X"    "1"      
#> [2,] "X->Y"    "1" "XtoY" NA    "2"       "X"   "Y"    "1"      
#> [3,] "Y->X"    "1" "YtoX" NA    "3"       "Y"   "X"    "1"      
#> [4,] "Y->Y"    "1" "YtoY" NA    "4"       "Y"   "Y"    "1"      
#> [5,] "X <-> X" "0" "sdX"  NA    "5"       "X"   "X"    "2"      
#> [6,] "Y <-> Y" "0" "sdY"  NA    "6"       "Y"   "Y"    "2"      
#> 
#> $ram
#>    heads to from parameter start
#> 1      1  2    1         1  <NA>
#> 2      1  3    2         1  <NA>
#> 3      1  4    3         1  <NA>
#> 5      1  6    1         2  <NA>
#> 6      1  7    2         2  <NA>
#> 7      1  8    3         2  <NA>
#> 9      1  2    5         3  <NA>
#> 10     1  3    6         3  <NA>
#> 11     1  4    7         3  <NA>
#> 13     1  6    5         4  <NA>
#> 14     1  7    6         4  <NA>
#> 15     1  8    7         4  <NA>
#> 17     2  1    1         5  <NA>
#> 18     2  2    2         5  <NA>
#> 19     2  3    3         5  <NA>
#> 20     2  4    4         5  <NA>
#> 21     2  5    5         6  <NA>
#> 22     2  6    6         6  <NA>
#> 23     2  7    7         6  <NA>
#> 24     2  8    8         6  <NA>
#> 
#> attr(,"class")
#> [1] "dsem_ram"

# Dynamic factor analysis with one factor and two manifest variables
# (specifies a random-walk for the factor, and miniscule residual SD)
dsem = "
  factor -> X, 0, loadings1
  factor -> Y, 0, loadings2
  factor -> factor, 1, NA, 1
  X <-> X, 0, NA, 0           # No additional variance
  Y <-> Y, 0, NA, 0           # No additional variance
"
make_dsem_ram( dsem=dsem, variables=c("X","Y","factor"), times=1:4 )
#> NOTE: adding 1 variances to the model
#> $model
#>                                              parameter first    second  
#> [1,] "factor->X"         "0" "loadings1" NA  "1"       "factor" "X"     
#> [2,] "factor->Y"         "0" "loadings2" NA  "2"       "factor" "Y"     
#> [3,] "factor->factor"    "1" NA          "1" "0"       "factor" "factor"
#> [4,] "X <-> X"           "0" NA          "0" "0"       "X"      "X"     
#> [5,] "Y <-> Y"           "0" NA          "0" "0"       "Y"      "Y"     
#> [6,] "factor <-> factor" "0" "V[factor]" NA  "3"       "factor" "factor"
#>      direction
#> [1,] "1"      
#> [2,] "1"      
#> [3,] "1"      
#> [4,] "2"      
#> [5,] "2"      
#> [6,] "2"      
#> 
#> $ram
#>    heads to from parameter start
#> 1      1  1    9         1  <NA>
#> 2      1  2   10         1  <NA>
#> 3      1  3   11         1  <NA>
#> 4      1  4   12         1  <NA>
#> 5      1  5    9         2  <NA>
#> 6      1  6   10         2  <NA>
#> 7      1  7   11         2  <NA>
#> 8      1  8   12         2  <NA>
#> 9      1 10    9         0     1
#> 10     1 11   10         0     1
#> 11     1 12   11         0     1
#> 13     2  1    1         0     0
#> 14     2  2    2         0     0
#> 15     2  3    3         0     0
#> 16     2  4    4         0     0
#> 17     2  5    5         0     0
#> 18     2  6    6         0     0
#> 19     2  7    7         0     0
#> 20     2  8    8         0     0
#> 21     2  9    9         3  <NA>
#> 22     2 10   10         3  <NA>
#> 23     2 11   11         3  <NA>
#> 24     2 12   12         3  <NA>
#> 
#> attr(,"class")
#> [1] "dsem_ram"

# ARIMA(1,1,0)
dsem = "
  factor -> factor, 1, rho1 # AR1 component
  X -> X, 1, NA, 1          # Integrated component
  factor -> X, 0, NA, 1
  X <-> X, 0, NA, 0         # No additional variance
"
make_dsem_ram( dsem=dsem, variables=c("X","factor"), times=1:4 )
#> NOTE: adding 1 variances to the model
#> $model
#>                                              parameter first    second  
#> [1,] "factor->factor"    "1" "rho1"      NA  "1"       "factor" "factor"
#> [2,] "X->X"              "1" NA          "1" "0"       "X"      "X"     
#> [3,] "factor->X"         "0" NA          "1" "0"       "factor" "X"     
#> [4,] "X <-> X"           "0" NA          "0" "0"       "X"      "X"     
#> [5,] "factor <-> factor" "0" "V[factor]" NA  "2"       "factor" "factor"
#>      direction
#> [1,] "1"      
#> [2,] "1"      
#> [3,] "1"      
#> [4,] "2"      
#> [5,] "2"      
#> 
#> $ram
#>    heads to from parameter start
#> 1      1  6    5         1  <NA>
#> 2      1  7    6         1  <NA>
#> 3      1  8    7         1  <NA>
#> 5      1  2    1         0     1
#> 6      1  3    2         0     1
#> 7      1  4    3         0     1
#> 9      1  1    5         0     1
#> 10     1  2    6         0     1
#> 11     1  3    7         0     1
#> 12     1  4    8         0     1
#> 13     2  1    1         0     0
#> 14     2  2    2         0     0
#> 15     2  3    3         0     0
#> 16     2  4    4         0     0
#> 17     2  5    5         2  <NA>
#> 18     2  6    6         2  <NA>
#> 19     2  7    7         2  <NA>
#> 20     2  8    8         2  <NA>
#> 
#> attr(,"class")
#> [1] "dsem_ram"

# ARIMA(0,0,1)
dsem = "
  factor -> X, 0, NA, 1
  factor -> X, 1, rho1     # MA1 component
  X <-> X, 0, NA, 0        # No additional variance
"
make_dsem_ram( dsem=dsem, variables=c("X","factor"), times=1:4 )
#> NOTE: adding 1 variances to the model
#> $model
#>                                              parameter first    second  
#> [1,] "factor->X"         "0" NA          "1" "0"       "factor" "X"     
#> [2,] "factor->X"         "1" "rho1"      NA  "1"       "factor" "X"     
#> [3,] "X <-> X"           "0" NA          "0" "0"       "X"      "X"     
#> [4,] "factor <-> factor" "0" "V[factor]" NA  "2"       "factor" "factor"
#>      direction
#> [1,] "1"      
#> [2,] "1"      
#> [3,] "2"      
#> [4,] "2"      
#> 
#> $ram
#>    heads to from parameter start
#> 1      1  1    5         0     1
#> 2      1  2    6         0     1
#> 3      1  3    7         0     1
#> 4      1  4    8         0     1
#> 5      1  2    5         1  <NA>
#> 6      1  3    6         1  <NA>
#> 7      1  4    7         1  <NA>
#> 9      2  1    1         0     0
#> 10     2  2    2         0     0
#> 11     2  3    3         0     0
#> 12     2  4    4         0     0
#> 13     2  5    5         2  <NA>
#> 14     2  6    6         2  <NA>
#> 15     2  7    7         2  <NA>
#> 16     2  8    8         2  <NA>
#> 
#> attr(,"class")
#> [1] "dsem_ram"
```
