
if( packageVersion("tinyVAST") < package_version("1.5.1.9000") ){
  devtools::install_github("vast-lib/tinyVAST@dev", force=TRUE, dep = FALSE)
  #devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force=TRUE)
}

library(tinyVAST)
library(fmesher)
library(sf)
library(units)

data( bering_sea_pollock_ages )

# Get shapefile for survey extent
data( bering_sea )

# Make extrapolation grid based on shapefile
bering_sea = st_transform(
  bering_sea,
  st_crs("+proj=utm +zone=2 +units=km")
)

# Experiment with rescaling ... doesn't play well with st_crs and st_within during tinyVAST(.)
#bering_sea = (bering_sea / 100) # - c(0, 60)

# subset to Years 2017-2023 (to speed up the example)
Data = subset( bering_sea_pollock_ages, Year == 2017 )

# Add Year-_Age interaction
Data$Age = factor( paste0("Age_",Data$Age) )
Data$Year_Age = interaction( Data$Year, Data$Age )

# Project data to UTM
Data = st_as_sf(
  Data,
  coords = c('Lon','Lat'),
  crs = st_crs(4326)
)

# Avoid projection
Data = st_transform(
  Data,
  crs = st_crs("+proj=utm +zone=2 +units=km")
)

# Add UTM coordinates as columns X & Y
Data = cbind( st_drop_geometry(Data), st_coordinates(Data) )

# adds different variances for each age
sem = ""

# Constant AR1 spatio-temporal term across ages
# and adds different variances for each age
dsem = "
  Age_1 -> Age_1, 1, lag1
  Age_2 -> Age_2, 1, lag1
  Age_3 -> Age_3, 1, lag1
  Age_4 -> Age_4, 1, lag1
  Age_5 -> Age_5, 1, lag1
  Age_6 -> Age_6, 1, lag1
  Age_7 -> Age_7, 1, lag1
  Age_8 -> Age_8, 1, lag1
  Age_9 -> Age_9, 1, lag1
  Age_10 -> Age_10, 1, lag1
  Age_11 -> Age_11, 1, lag1
  Age_12 -> Age_12, 1, lag1
  Age_13 -> Age_13, 1, lag1
  Age_14 -> Age_14, 1, lag1
  Age_15 -> Age_15, 1, lag1
"

# Separate intercept for each Year-Age combination
Formula = Abundance_per_hectare ~ 0 + Age

#
control = tinyVASTcontrol(
  getsd = FALSE,
  #profile = c("alpha_j","alpha2_j"),
  trace = 1,
  profile = NULL,
  use_anisotropy = TRUE
)

# Define separate tweedie family for each age
Family = list(
  Age_1 = delta_gamma(type = "poisson-link"),
  Age_2 = delta_gamma(type = "poisson-link"),
  Age_3 = delta_gamma(type = "poisson-link"),
  Age_4 = delta_gamma(type = "poisson-link"),
  Age_5 = delta_gamma(type = "poisson-link"),
  Age_6 = delta_gamma(type = "poisson-link"),
  Age_7 = delta_gamma(type = "poisson-link"),
  Age_8 = delta_gamma(type = "poisson-link"),
  Age_9 = delta_gamma(type = "poisson-link"),
  Age_10 = delta_gamma(type = "poisson-link"),
  Age_11 = delta_gamma(type = "poisson-link"),
  Age_12 = delta_gamma(type = "poisson-link"),
  Age_13 = delta_gamma(type = "poisson-link"),
  Age_14 = delta_gamma(type = "poisson-link"),
  Age_15 = delta_gamma(type = "poisson-link")
)

#########################
# SPDE method
#########################

# Make spatial domain
mesh = fm_mesh_2d(
  loc = Data[,c("X","Y")],
  cutoff = 100
)

# Fit model
myfit = tinyVAST(
  data = Data,
  formula = Formula,
  space_term = sem,
  #spacetime_term = dsem,
  family = Family,
  delta_options = list(
    formula = ~ 0 + Age,
    #spacetime_term = dsem,
    space_term = sem
  ),
  space_column = c("X", "Y"),
  variable_column = "Age",
  #time_column = "Year",
  distribution_column = "Age",
  spatial_domain = mesh,
  control = control
)

# Make extrapolation grid based on shapefile
grid = st_make_grid( bering_sea, n=c(50,50) )
grid = st_intersection( grid, bering_sea )
grid = st_make_valid( grid )
loc_gz = st_coordinates(st_centroid( grid ))

# Get area for extrapolation grid
area_g = set_units(st_area(grid), "hectares") #  / 100^2 # Hectares

# Get abundance
N_jz = expand.grid( Age=myfit$internal$variables, Year=sort(unique(Data$Year)) )
N_jz = cbind( N_jz, "Biomass"=NA, "SE"=NA )

# Make newdata in blocks
areas = newdata = NULL
for( j in seq_len(nrow(N_jz)) ){
  tmp = data.frame( loc_gz, Year=N_jz[j,'Year'], Age=N_jz[j,'Age'])
  newdata = rbind( newdata, cbind(tmp, block = j) )
  areas = c( areas, area_g )
}

# Debug
if( FALSE ){
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R)' )
  source("utility.R")
  source("internal.R")
  object = myfit
  #newdata
  area = areas
  block = rep(1,nrow(newdata))
  type = rep(1,nrow(newdata))
  weighting_index = rep(0, nrow(newdata))
  covariate = rep(0, nrow(newdata))
  getsd = TRUE
  bias.correct = TRUE
  apply.epsilon = FALSE
  intern = FALSE
}

# Area-expansion
index1 = integrate_output(
  myfit,
  block = newdata$block,
  area = areas,
  newdata = newdata,
  apply.epsilon = TRUE,
  bias.correct = FALSE,
  intern = FALSE,
  getsd = FALSE
)


  if( N_jz[j,'Age']==1 ){
    message( "Integrating ", N_jz[j,'Year'], " ", N_jz[j,'Age'], ": ", Sys.time() )
  }
  if( is.na(N_jz[j,'Biomass']) ){
    # Area-expansion
    index1 = integrate_output( myfit,
                    area = areas,
                    newdata = newdata,
                    apply.epsilon = TRUE,
                    bias.correct = FALSE,
                    intern = FALSE,
                    getsd = FALSE )
    #N_jz[j,'Biomass'] = index1[3] / 1e9
    N_jz[j,'Biomass'] = index1[1] / 1e9
  }
}
run_time = Sys.time() - start_time
N_ct = array(
  N_jz$Biomass,
  dim = c(length(myfit$internal$variables),length(unique(Data$Year))),
  dimnames = list(myfit$internal$variables,sort(unique(Data$Year)))
)
N_ct = sweep( N_ct, MARGIN = 1, FUN = "/", STATS = colSums(N_ct) )

