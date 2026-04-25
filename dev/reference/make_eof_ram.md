# Make a RAM (Reticular Action Model)

`make_eof_ram` converts SEM arrow notation to `ram` describing SEM
parameters

## Usage

``` r
make_eof_ram(
  times,
  variables,
  n_eof,
  remove_na = TRUE,
  standard_deviations = "unequal"
)
```

## Arguments

- times:

  A character vector listing the set of times in order

- variables:

  A character vector listing the set of variables

- n_eof:

  Number of EOF modes of variability to estimate

- remove_na:

  Boolean indicating whether to remove NA values from RAM (default) or
  not. `remove_NA=FALSE` might be useful for exploration and diagnostics
  for advanced users

- standard_deviations:

  One of `"equal"`, `"unequal"`, or a numeric vector indicating fixed
  values.

## Value

A reticular action module (RAM) describing dependencies

## Examples

``` r
# Two EOFs for two variables
make_eof_ram( times = 2010:2020, variables = c("pollock","cod"), n_eof=2 )
#> $model
#>      to  from parameter
#> 1  2010 EOF_1         1
#> 2  2011 EOF_1         2
#> 3  2012 EOF_1         3
#> 4  2013 EOF_1         4
#> 5  2014 EOF_1         5
#> 6  2015 EOF_1         6
#> 7  2016 EOF_1         7
#> 8  2017 EOF_1         8
#> 9  2018 EOF_1         9
#> 10 2019 EOF_1        10
#> 11 2020 EOF_1        11
#> 12 2010 EOF_2        NA
#> 13 2011 EOF_2        12
#> 14 2012 EOF_2        13
#> 15 2013 EOF_2        14
#> 16 2014 EOF_2        15
#> 17 2015 EOF_2        16
#> 18 2016 EOF_2        17
#> 19 2017 EOF_2        18
#> 20 2018 EOF_2        19
#> 21 2019 EOF_2        20
#> 22 2020 EOF_2        21
#> 
#> $ram
#>       heads to from parameter start
#>  [1,]     1  3    1         1  0.01
#>  [2,]     1  4    1         2  0.01
#>  [3,]     1  5    1         3  0.01
#>  [4,]     1  6    1         4  0.01
#>  [5,]     1  7    1         5  0.01
#>  [6,]     1  8    1         6  0.01
#>  [7,]     1  9    1         7  0.01
#>  [8,]     1 10    1         8  0.01
#>  [9,]     1 11    1         9  0.01
#> [10,]     1 12    1        10  0.01
#> [11,]     1 13    1        11  0.01
#> [12,]     1  4    2        12  0.01
#> [13,]     1  5    2        13  0.01
#> [14,]     1  6    2        14  0.01
#> [15,]     1  7    2        15  0.01
#> [16,]     1  8    2        16  0.01
#> [17,]     1  9    2        17  0.01
#> [18,]     1 10    2        18  0.01
#> [19,]     1 11    2        19  0.01
#> [20,]     1 12    2        20  0.01
#> [21,]     1 13    2        21  0.01
#> [22,]     1 16   14         1  0.01
#> [23,]     1 17   14         2  0.01
#> [24,]     1 18   14         3  0.01
#> [25,]     1 19   14         4  0.01
#> [26,]     1 20   14         5  0.01
#> [27,]     1 21   14         6  0.01
#> [28,]     1 22   14         7  0.01
#> [29,]     1 23   14         8  0.01
#> [30,]     1 24   14         9  0.01
#> [31,]     1 25   14        10  0.01
#> [32,]     1 26   14        11  0.01
#> [33,]     1 17   15        12  0.01
#> [34,]     1 18   15        13  0.01
#> [35,]     1 19   15        14  0.01
#> [36,]     1 20   15        15  0.01
#> [37,]     1 21   15        16  0.01
#> [38,]     1 22   15        17  0.01
#> [39,]     1 23   15        18  0.01
#> [40,]     1 24   15        19  0.01
#> [41,]     1 25   15        20  0.01
#> [42,]     1 26   15        21  0.01
#> [43,]     2  1    1         0  1.00
#> [44,]     2  2    2         0  1.00
#> [45,]     2  3    3        22    NA
#> [46,]     2  4    4        22    NA
#> [47,]     2  5    5        22    NA
#> [48,]     2  6    6        22    NA
#> [49,]     2  7    7        22    NA
#> [50,]     2  8    8        22    NA
#> [51,]     2  9    9        22    NA
#> [52,]     2 10   10        22    NA
#> [53,]     2 11   11        22    NA
#> [54,]     2 12   12        22    NA
#> [55,]     2 13   13        22    NA
#> [56,]     2 14   14         0  1.00
#> [57,]     2 15   15         0  1.00
#> [58,]     2 16   16        23    NA
#> [59,]     2 17   17        23    NA
#> [60,]     2 18   18        23    NA
#> [61,]     2 19   19        23    NA
#> [62,]     2 20   20        23    NA
#> [63,]     2 21   21        23    NA
#> [64,]     2 22   22        23    NA
#> [65,]     2 23   23        23    NA
#> [66,]     2 24   24        23    NA
#> [67,]     2 25   25        23    NA
#> [68,]     2 26   26        23    NA
#> 
#> $variances
#>        to    from parameter
#> 1   EOF_1   EOF_1         0
#> 2   EOF_2   EOF_2         0
#> 3 pollock pollock        22
#> 4     cod     cod        23
#> 
#> $standard_deviations
#> [1] "unequal"
#> 
#> attr(,"class")
#> [1] "eof_ram"
```
