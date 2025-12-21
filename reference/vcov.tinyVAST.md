# Extract Variance-Covariance Matrix

extract the covariance of fixed effects, or both fixed and random
effects.

## Usage

``` r
# S3 method for class 'tinyVAST'
vcov(object, which = c("fixed", "random", "both"), ...)
```

## Arguments

- object:

  output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/reference/tinyVAST.md)

- which:

  whether to extract the covariance among fixed effects, random effects,
  or both

- ...:

  ignored, for method compatibility

## Value

A square matrix containing the estimated covariances among the parameter
estimates in the model. The dimensions dependend upon the argument
`which`, to determine whether fixed, random effects, or both are
outputted.
