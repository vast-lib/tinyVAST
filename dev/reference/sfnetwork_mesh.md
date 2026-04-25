# Make mesh for stream network

make an object representing spatial information required to specify a
stream-network spatial domain, similar in usage to
`link[fmesher]{fm_mesh_2d}` for a 2-dimensional continuous domain

## Usage

``` r
sfnetwork_mesh(stream)
```

## Arguments

- stream:

  sfnetworks object representing stream network

## Value

An object (list) of class `sfnetwork_mesh`. Elements include:

- N:

  The number of random effects used to represent the network

- table:

  a table containing a description of parent nodes (from), childen nodes
  (to), and the distance separating them

- stream:

  copy of the stream network object passed as argument
