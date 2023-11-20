
## START NEW
root_dir = R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\data-raw)'
## END NEW

library(ncdf4)     # for reading NetCDF
library(lubridate) # for month(.)

url <- 'https://polarwatch.noaa.gov/erddap/griddap/'

#
Date = Sys.Date()
month_to_use = c( "march", "september" )[2]
run_dir = paste0( root_dir, Date, "_", month_to_use, "/" )
  dir.create(run_dir)

################
# Download initial dimensions
################


# Download indices
grid_id <- 'nsidcCDRice_nh_grid'
if( !("grid.nc" %in% list.files(run_dir) ) ){
  grid_urlcall <- paste0(url,grid_id,'.nc?longitude[(5812500.0):1:(-5337500.0)][(-3837500.0):1:(3737500.0)],latitude[(5812500.0):1:(-5337500.0)][(-3837500.0):1:(3737500.0)]')
  grid_nc <- download.file(grid_urlcall, destfile=paste0(run_dir,"grid.nc"), mode='wb')
}

# Read the grid file
gridFid <- nc_open( paste0(run_dir,'grid.nc') )
  ygrid <- ncvar_get(gridFid, varid="ygrid")
  xgrid <- ncvar_get(gridFid, varid="xgrid")
  longitude <- ncvar_get(gridFid, varid="longitude")
  latitude <- ncvar_get(gridFid, varid="latitude")
nc_close(gridFid)

# Find the x and y indices that correspond with your area of interest
# Next we will get the indices of our area of interest. We will use these indices later to
# make subsetted data requests. For this example we are interested in all data north of 75°.
# We will find all coordinate values greater than or equal to 75°N. To subset by longitude
# you can add a longitude query. The range is an extent that covers all data north of 75°
# so we will be requesting the smallest box of data that covers our desired latitude range.
# Note that our returned request will include points south of 75°N to accomplish this.

inds = which(latitude > 57, arr.ind=TRUE)
rowrange <- range(inds[,1])
colrange <- range(inds[,2])

################
# Extract the data
################

# Request the sea ice dataset using the selected indices from the lat/lon grid
# Request data from the monthly science quality dataset for the Northern Hemisphere
# (nsidcCDRiceSQnhmday). There are four different sea ice variables in this dataset.
# Downloading each of them requires adding the name of each variable and the date and
# coordinate constraints. Here we download two of the variables. If you need a refresher
# on the structure of the URL call, go to the ERDDAP "Data Access Form" for a dataset and
# use the 'generate the URL' button.
varnames <- c('seaice_conc_monthly_cdr', 'goddard_merged_seaice_conc_monthly')

spatial_scale = 5 # 5 has about 100 km x 100 km scale

# Generate a data request URL using the indices from the grid
dataid <- 'nsidcCDRiceSQnhmday'
datestring <- '[(1997-01-16T00:00:00Z):1:(2017-12-16T00:00:00Z)]'
coordstring <- paste0('[',colrange[1],':',spatial_scale,':',colrange[2],'][',rowrange[1],':',spatial_scale,':',rowrange[2],']')

if( !("data.nc" %in% list.files(run_dir) ) ){

  for (i in 1:length(varnames)) {
    if (i == 1) {
      urlcall <- paste0(url,dataid,'.nc?',varnames[i],datestring,coordstring)
    }
    else {
      urlcall <- paste0(urlcall,',',varnames[i],datestring,coordstring)
    }
  }

  # Download the netCDF file  (this will take a few minutes, 20 years of data)
  data_nc <- download.file(urlcall, destfile=paste0(run_dir,"data.nc"), mode='wb')
}

# Read the downloaded netCDF file and load the ice data into R variables
dataFid <- nc_open( paste0(run_dir,'data.nc') )
  datatime <- ncvar_get(dataFid, varid="time")
  datatime <- as.Date(as.POSIXlt(datatime,origin='1970-01-01',tz= "GMT"))
  ygrid <- ncvar_get(dataFid, varid="ygrid")
  xgrid <- ncvar_get(dataFid, varid="xgrid")
  #longitude <- ncvar_get(dataFid, varid="longitude")
  #latitude <- ncvar_get(dataFid, varid="latitude")
  seaiceCDR <- ncvar_get(dataFid, varid=varnames[1])
  seaiceGoddard <- ncvar_get(dataFid, varid=varnames[2])
nc_close(dataFid)

# Which month has highest -- MARCH
# Which month has lowest -- SEPTEMBER
SUM = apply( seaiceCDR, MARGIN=3, FUN=sum, na.rm=TRUE)
tapply( SUM, INDEX=month(datatime), FUN=mean ) / 10000

# Subset
seaiceCDR = seaiceCDR[,,which(tolower(months(datatime))==tolower(month_to_use))]
datatime = datatime[which(tolower(months(datatime))==tolower(month_to_use))]

# Replace 2.53 , 2.54, 2.55 with NAs
seaiceCDR = ifelse( seaiceCDR > 1, NA, seaiceCDR )

################
# Extract Lat-Lon for data
################

# Request a grid subset using the same coordinate string used for the data download.
if( !("grid_subset.nc" %in% list.files(run_dir) ) ){
  urlcall <- paste0(url,grid_id,'.nc?longitude',coordstring,',latitude',coordstring)
  grid_subset <- download.file(urlcall, destfile=paste0(run_dir,"grid_subset.nc"), mode='wb')
}

# Read and format the subsetted grid data from the netCDF file
gridSubsetFid <- nc_open( paste0(run_dir,'grid_subset.nc') )
  ygrid <- ncvar_get(gridSubsetFid, varid="ygrid")
  xgrid <- ncvar_get(gridSubsetFid, varid="xgrid")
  longitudeSubset <- ncvar_get(gridSubsetFid, varid="longitude")
  latitudeSubset <- ncvar_get(gridSubsetFid, varid="latitude")
nc_close(gridSubsetFid)

# Choose a date to use for the map and determine the index value for the date

# Make a long-format dataframe for that time period to use with ggplot
sea_ice <- data.frame(
               lon = rep( as.vector(longitudeSubset), length(datatime) ),
               lat = rep( as.vector(latitudeSubset), length(datatime) ),
               year = rep( year(datatime), each=prod(dim(latitudeSubset)) ),
               ice_concentration = as.vector(seaiceCDR)
             )
sea_ice = na.omit(sea_ice)

################
# SImplify output format
################


# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
usethis::use_data( sea_ice, overwrite=TRUE )
