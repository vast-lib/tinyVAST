# Make spatial object for NNGP

Convert an sf areal model object to a class `"nngp_domain"` that is
recognized by tinyVAST for nearest neighbors Gaussian process models

## Usage

``` r
make_nngp_domain(sf_areal, nn)
```

## Arguments

- sf_areal:

  an sf or sfc object containing "POLYGON" and "MULTIPOLYGON" types

- nn:

  Number of nearest neighbors

## Value

a list of NNGP data objects created by tinyVAST:::make_nngp_data, as
well as the provided sf object
