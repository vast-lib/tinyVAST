# Additional families

Additional families compatible with
[`tinyVAST()`](https://vast-lib.github.io/tinyVAST/reference/tinyVAST.md).

## Usage

``` r
delta_lognormal(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_gamma(link1, link2 = "log", type = c("standard", "poisson-link"))
```

## Arguments

- link1:

  Link for first part of delta/hurdle model.

- link2:

  Link for second part of delta/hurdle model.

- type:

  Delta/hurdle family type. `"standard"` for a classic hurdle model.
  `"poisson-link"` for a Poisson-link delta model (Thorson 2018).

- link:

  Link.

## Value

A list with elements common to standard R family objects including
`family`, `link`, `linkfun`, and `linkinv`. Delta/hurdle model families
also have elements `delta` (logical) and `type` (standard vs.
Poisson-link).

## References

*Poisson-link delta families*:

Thorson, J.T. 2018. Three problems with the conventional delta-model for
biomass sampling data, and a computationally efficient alternative.
Canadian Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
[doi:10.1139/cjfas-2017-0266](https://doi.org/10.1139/cjfas-2017-0266)

*Poisson-link delta families*:

Thorson, J.T. 2018. Three problems with the conventional delta-model for
biomass sampling data, and a computationally efficient alternative.
Canadian Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
[doi:10.1139/cjfas-2017-0266](https://doi.org/10.1139/cjfas-2017-0266)
