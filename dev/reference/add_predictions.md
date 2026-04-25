# Add predictions to data-list

Given user-provided `newdata`, expand the object `tmb_data` to include
predictions corresponding to those new observations

## Usage

``` r
add_predictions(object, newdata, remove_origdata = FALSE)
```

## Arguments

- object:

  Output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/dev/reference/tinyVAST.md).

- newdata:

  New data-frame of independent variables used to predict the response.

- remove_origdata:

  Whether to remove original-data to allow faster evaluation.
  `remove_origdata=TRUE` eliminates information about the distribution
  for random effects, and cannot be combined with epsilon
  bias-correction. WARNING: feature is experimental and subject to
  change.

## Value

the object `fit$tmb_inputs$tmb_data` representing data used during
fitting, but with updated values for slots associated with predictions,
where this updated object can be recompiled by TMB to provide
predictions
