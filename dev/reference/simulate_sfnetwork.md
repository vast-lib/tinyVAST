# Simulate GMRF for stream network

Simulate values from a GMRF using a tail-down (flow-unconnected)
exponential model on a stream network

## Usage

``` r
simulate_sfnetwork(sfnetwork_mesh, theta, n = 1, what = c("samples", "Q"))
```

## Arguments

- sfnetwork_mesh:

  Output from
  [`sfnetwork_mesh`](https://vast-lib.github.io/tinyVAST/dev/reference/sfnetwork_mesh.md)

- theta:

  Decorrelation rate

- n:

  number of simulated GMRFs

- what:

  Whether to return the simulated GMRF or its precision matrix

## Value

a matrix of simulated values for a Gaussian Markov random field arising
from a stream-network spatial domain, with row for each spatial random
effect and `n` columns, using the sparse precision matrix defined in
Charsley et al. (2023)

## References

Charsley, A. R., Gruss, A., Thorson, J. T., Rudd, M. B., Crow, S. K.,
David, B., Williams, E. K., & Hoyle, S. D. (2023). Catchment-scale
stream network spatio-temporal models, applied to the freshwater stages
of a diadromous fish species, longfin eel (Anguilla dieffenbachii).
Fisheries Research, 259, 106583.
[doi:10.1016/j.fishres.2022.106583](https://doi.org/10.1016/j.fishres.2022.106583)
