

######################
# Provide areal domain
######################

library(sf)

boundary = st_polygon(list(cbind(c(0,1,1,0,0),c(0,0,1,1,0))))
spatial_domain = st_make_grid( boundary, square = FALSE, n = c(4,4) )

######################
# Make NNGP neighborhoods
######################

data = make_nngp_data(
  coords = st_coordinates(st_centroid(spatial_domain)),
  nn = 4
)

