
# Load packages
library(VAST)
library(sf)

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="multimodal_red_snapper" )

# 
red_snapper = example$sampling_data

shape_dir = R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_7)'
red_snapper_shapefile = read_sf( file.path(shape_dir,"red_snapper_extent.shp") )

setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
usethis::use_data( red_snapper )
usethis::use_data( red_snapper_shapefile )
