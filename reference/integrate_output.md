# Integration for target variable

Calculates an estimator for a derived quantity by summing across
multiple predictions. This can be used to approximate an integral when
estimating area-expanded abundance, abundance-weighting a covariate to
calculate distribution shifts, and/or weighting one model variable by
another.

## Usage

``` r
integrate_output(
  object,
  newdata,
  area,
  type = rep(1, nrow(newdata)),
  weighting_index,
  covariate,
  getsd = TRUE,
  bias.correct = TRUE,
  apply.epsilon = FALSE,
  intern = FALSE
)
```

## Arguments

- object:

  Output from
  [`tinyVAST()`](https://vast-lib.github.io/tinyVAST/reference/tinyVAST.md).

- newdata:

  New data-frame of independent variables used to predict the response,
  where a total value is calculated by combining across these individual
  predictions. If these locations are randomly drawn from a specified
  spatial domain, then `integrate_output` applies midpoint integration
  to approximate the total over that area. If locations are drawn
  sysmatically from a domain, then `integrate_output` is applying a
  midpoint approximation to the integral.

- area:

  vector of values used for area-weighted expansion of estimated density
  surface for each row of `newdata` with length of `nrow(newdata)`.

- type:

  Integer-vector indicating what type of expansion to apply to each row
  of `newdata`, with length of `nrow(newdata)`.

  `type=1`

  :   Area-weighting: weight predictor by argument `area`

  `type=2`

  :   Abundance-weighted covariate: weight `covariate` by proportion of
      total in each row of `newdata`

  `type=3`

  :   Abundance-weighted variable: weight predictor by proportion of
      total in a prior row of `newdata`. This option is used to weight a
      prediction for one category based on predicted proportional
      density of another category, e.g., to calculate abundance-weighted
      condition in a bivariate model.

  `type=4`

  :   Abundance-expanded variable: weight predictor by density in a
      prior row of `newdata`. This option is used to weight a prediction
      for one category based on predicted density of another category,
      e.g., to calculate abundance-expanded consumption in a bivariate
      model.

  `type=0`

  :   Exclude from weighting: give weight of zero for a given row of
      `newdata`. Including a row of `newdata` with `type=0` is useful,
      e.g., when calculating abundance at that location, where the
      eventual index uses abundance as weighting term but without
      otherwise using the predicted density in calculating a total
      value.

- weighting_index:

  integer-vector used to indicate a previous row that is used to
  calculate a weighted average that is then applied to the given row of
  `newdata`. Only used for when `type=3`.

- covariate:

  numeric-vector used to provide a covariate that is used in expansion,
  e.g., to provide positional coordinates when calculating the
  abundance-weighted centroid with respect to that coordinate. Only used
  for when `type=2`.

- getsd:

  logical indicating whether to get the standard error, where
  `getsd=FALSE` is faster during initial exploration

- bias.correct:

  logical indicating if bias correction should be applied using standard
  methods in
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html)

- apply.epsilon:

  Apply epsilon bias correction using a manual calculation rather than
  using the conventional method in
  [TMB::sdreport](https://rdrr.io/pkg/TMB/man/sdreport.html)? See
  details for more information.

- intern:

  Do Laplace approximation on C++ side? Passed to
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

## Value

A vector containing the plug-in estimate, standard error, the epsilon
bias-corrected estimate if available, and the standard error for the
bias-corrected estimator. Depending upon settings, one or more of these
will be `NA` values, and the function can be repeatedly called to get
multiple estimators and/or statistics.

## Details

Analysts will often want to calculate some value by combining the
predicted response at multiple locations, and potentially from multiple
variables in a multivariate analysis. This arises in a univariate model,
e.g., when calculating the integral under a predicted density function,
which is approximated using a midpoint or Monte Carlo approximation by
calculating the linear predictors at each location `newdata`, applying
the inverse-link-trainsformation, and calling this predicted response
`mu_g`. Total abundance is then be approximated by multiplying `mu_g` by
the area associated with each midpoint or Monte Carlo approximation
point (supplied by argument `area`), and summing across these
area-expanded values.

In more complicated cases, an analyst can then use `covariate` to
calculate the weighted average of a covariate for each midpoint
location. For example, if the covariate is positional coordinates or
depth/elevation, then `type=2` measures shifts in the average habitat
utilization with respect to that covariate. Alternatively, an analyst
fitting a multivariate model might weight one variable based on another
using `weighting_index`, e.g., to calculate abundance-weighted average
condition, or predator-expanded stomach contents.

In practice, spatial integration in a multivariate model requires two
passes through the rows of `newdata` when calculating a total value. In
the following, we write equations using C++ indexing conventions such
that indexing starts with 0, to match the way that `integrate_output`
expects indices to be supplied. Given inverse-link-transformed predictor
\\ \mu_g \\, function argument `type` as \\ type_g \\ function argument
`area` as \\ a_g \\, function argument `covariate` as \\ x_g \\,
function argument `weighting_index` as `\eqn{ h_g }` function argument
`weighting_index` as `\eqn{ h_g }` the first pass calculates:

\$\$ \nu_g = \mu_g a_g \$\$

where the total value from this first pass is calculated as:

\$\$ \nu^\* = \sum\_{g=0}^{G-1} \nu_g \$\$

The second pass then applies a further weighting, which depends upon \\
type_g \\, and potentially upon \\ x_g \\ and \\ h_g \\.

If \\type_g = 0\\ then \\\phi_g = 0\\

If \\type_g = 1\\ then \\\phi_g = \nu_g\\

If \\type_g = 2\\ then \\\phi_g = x_g \frac{\nu_g}{\nu^\*} \\

If \\type_g = 3\\ then \\\phi_g = \frac{\nu\_{h_g}}{\nu^\*} \mu_g \\

If \\type_g = 4\\ then \\\phi_g = \nu\_{h_g} \mu_g \\

Finally, the total value from this second pass is calculated as:

\$\$ \phi^\* = \sum\_{g=0}^{G-1} \phi_g \$\$

and \\\phi^\*\\ is outputted by `integrate_output`, along with a
standard error and potentially using the epsilon bias-correction
estimator to correct for skewness and retransformation bias.

Standard bias-correction using `bias.correct=TRUE` can be slow, and in
some cases it might be faster to do `apply.epsilon=TRUE` and
`intern=TRUE`. However, that option is somewhat experimental, and a user
might want to confirm that the two options give identical results.
Similarly, using `bias.correct=TRUE` will still calculate the
standard-error, whereas using `apply.epsilon=TRUE` and `intern=TRUE`
will not.
