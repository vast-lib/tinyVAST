
library(VAST)

data(condition_and_density_example)

condition = subset( condition_and_density_example, !is.na(Individual_length_cm) )
condition = condition[,c('Year',"Lat","Lon","Individual_length_cm","Individual_weight_Grams")]

density = subset( condition_and_density_example, is.na(Individual_length_cm) )
density = density[,c('Year',"Lat","Lon","Sample_biomass_KG")]
colnames(density)[4] = "Sample_biomass_KGperHectare"

#
library(sf)
eastern_bering_sea = st_read( R'(C:\Users\James.Thorson\Desktop\Git\FishStatsUtils\inst\region_shapefiles\EBSshelf)' )
eastern_bering_sea = st_geometry(eastern_bering_sea)
eastern_bering_sea = st_transform( eastern_bering_sea, crs = 4326 )

condition_and_density = list(
  "condition" = condition,
  "density" = density,
  "eastern_bering_sea" = eastern_bering_sea
)

setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
usethis::use_data( condition_and_density, overwrite=TRUE )
