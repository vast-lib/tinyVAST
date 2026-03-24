

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
######################
# Provide areal domain
######################

library(sf)

boundary = st_polygon(list(cbind(c(0,1,1,0,0),c(0,0,1,1,0))))
spatial_domain = st_make_grid( boundary, square = FALSE, n = c(4,4) )

######################
# Simulate data from SAR
######################

library(Matrix)

st_adjacent <- function(m, ...) st_relate(m, m, pattern="F***1****", ...)
grid_A = st_adjacent(spatial_domain, sparse=TRUE)
A_ss = as(grid_A, "sparseMatrix")
A_ss = sweep(A_ss, MARGIN=1, FUN = "/", STAT = rowSums(A_ss))
IminusA_ss = Diagonal( n = length(spatial_domain) ) - 0.9 * A_ss
Q_ss = t(IminusA_ss) %*% IminusA_ss

omega_i = RTMB:::rgmrf0( n = 1, Q_ss )[,1]
plot( st_sf(spatial_domain,omega_i) )
p_i = 3 + omega_i - mean(omega_i)
y_i = rpois( n = length(spatial_domain), lambda = exp(p_i) )

##################
# Run via package
##################

library(tinyVAST)

data = data.frame(
  response = y_i,
  setNames(data.frame(st_coordinates(st_centroid(spatial_domain))), c("x","y"))
)

fit = tinyVAST(
  formula = response ~ 1,
  data = data,
  space_term = "",
  spatial_domain = spatial_domain,
  family = poisson(),
  control = tinyVASTcontrol(
    nearest_neighbors = 4,
    gmrf_parameterization = "projection",
    extra_reporting = TRUE
  )
)


