

# Test compile locally
if( FALSE ){
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)' )
  TMB::compile('tinyVAST.cpp', framework = "TMBad")
}

# Document
if( FALSE ){
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
  devtools::document()
}

# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

######################
# Provide areal domain
######################

library(sf)

boundary = st_polygon(list(cbind(c(0,1,1,0,0),c(0,0,1,1,0))))
sf_areal = st_make_grid( boundary, square = FALSE, n = c(4,4) )

######################
# Simulate data from SAR
######################

library(Matrix)

st_adjacent <- function(m, ...) st_relate(m, m, pattern="F***1****", ...)
grid_A = st_adjacent(sf_areal, sparse=TRUE)
A_ss = as(grid_A, "sparseMatrix")
A_ss = sweep(A_ss, MARGIN=1, FUN = "/", STAT = rowSums(A_ss))
IminusA_ss = Diagonal( n = length(sf_areal) ) - 0.9 * A_ss
Q_ss = t(IminusA_ss) %*% IminusA_ss

omega_i = RTMB:::rgmrf0( n = 1, Q_ss )[,1]
plot( st_sf(sf_areal,omega_i) )
p_i = 3 + omega_i - mean(omega_i)
y_i = rpois( n = length(sf_areal), lambda = exp(p_i) )

##################
# Run via package
##################

library(tinyVAST)

spatial_domain = make_nngp_domain(
  sf_areal = sf_areal,
  nn = 4
)

data = data.frame(
  response = y_i,
  setNames(data.frame(st_coordinates(st_centroid(sf_areal))), c("x","y"))
)

fit = tinyVAST(
  formula = response ~ 1,
  data = data,
  space_term = "",
  spatial_domain = spatial_domain,
  family = poisson(),
  control = tinyVASTcontrol(
    gmrf_parameterization = "projection",
    extra_reporting = TRUE
  )
)


