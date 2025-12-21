
library(VAST)
library(sf)

example = load_example( data_set="NWA_yellowtail_seasons" )
atlantic_yellowtail = example$sampling_data[,c("year","season","latitude","longitude","swept","weight")]


#
sf_locs = st_geometry(st_as_sf( 
  atlantic_yellowtail, coords = c("longitude",'latitude') )
)
sf_locs = st_sfc(
  sf_locs, 
  crs = st_crs("EPSG:4326")
)

# Make domain using 10km buffer around samples
atlantic_yellowtail_domain = st_union(st_buffer( 
  sf_locs,
  dist = 10000
))


setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
usethis::use_data( atlantic_yellowtail )
usethis::use_data( atlantic_yellowtail_domain )
