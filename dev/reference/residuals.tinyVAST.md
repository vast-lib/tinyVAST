# Calculate deviance or response residuals for tinyVAST

Calculate residuals

## Usage

``` r
# S3 method for class 'tinyVAST'
residuals(object, type = c("deviance", "response"), ...)
```

## Arguments

- object:

  Output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md)

- type:

  which type of residuals to compute (only option is `"deviance"` or
  `"response"` for now)

- ...:

  Note used

## Value

a vector residuals, associated with each row of `data` supplied during
fitting
