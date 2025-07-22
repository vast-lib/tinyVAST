
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

library(tinyVAST)
library(ggplot2)
library(gridExtra)

# pak::pkg_install("DFO-NOAA-Pacific/surveyjoin")
library(surveyjoin)
library(sf)
#cache_data()

case = c( "pcod_2011", "GOA_pcod", "NorPac_pcod" )[1]

#library(sf)
if( case == "pcod_2011" ){
  #data( pcod, package = "sdmTMB" )
  data( qcs_grid, package = "sdmTMB" )
  data( pcod_2011, package = "sdmTMB" )

  # Rename
  bathy = data.frame( qcs_grid[,c('X','Y')], covar= qcs_grid[,'depth'] )
  pcod_2011 = as.data.frame(pcod_2011)
  data = data.frame( pcod_2011[,c('X','Y','density','year')], covar = pcod_2011[,'depth'] )

  # Piped version
  mesh <- fmesher::fm_mesh_2d(pcod_2011[, c("X", "Y")], cutoff = 5)
  mesh_with_covs <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = "covar",
    coords  = c("X", "Y")
  )
  xy_i = matrix( c( 500, 5685 ), nrow = 1 )
}
if( case == "GOA_pcod" ){
  load_sql_data()
  all_data <- get_data("pacific cod")
  all_data = as.data.frame(all_data)

  # Covariates
  load( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariates for log_tau in SPDE\GOA_bathy.rda)')
  GOA_bathy = terra::unwrap(GOA_bathy)

  # Subset
  data = subset( all_data, survey_name %in% c('Gulf of Alaska') )
  data = data.frame( X = data[,'lon_start'], Y = data[,'lat_start'], year = data[,'year'],
                     density = data[,'catch_weight'],
                     covar = data[,'depth_m'] )
  data[,c("X","Y")] = sf_project( to = st_crs(GOA_bathy), from = st_crs(4326), pts = data[,c("X","Y")] ) / 1000

  # Downsample bathymetry
  GOA_bathy = terra::aggregate(GOA_bathy, fact = 5)

  #
  bathy = as.data.frame(GOA_bathy, xy = TRUE)
  #bathy = sf_project( from = st_crs(GOA_bathy), to = st_crs(4326), pts = bathy )
  #bathy = as.data.frame(bathy)
  colnames(bathy) = c( "X", "Y", "covar" )
  bathy[,c("X","Y")] = bathy[,c("X","Y")] / 1000

  # Make above-land very different
  bathy$covar = ifelse( bathy$covar < 0, bathy$covar * 1000, bathy$covar )

  # Make mesh
  mesh <- fmesher::fm_mesh_2d(data[, c("X", "Y")], cutoff = 25 )  # 0.5 for Lon-Lat
  mesh_with_covs <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = "covar",
    coords  = c("X", "Y")
  )
  terra::plot( GOA_bathy )
  xy_i = data.frame( -88863, 785481  ) / 1000
  #xy_i = sf_project( from = st_crs(GOA_bathy), to = st_crs(4326), pts = xy_i )
}
if( case == "NorPac_pcod" ){
  load_sql_data()
  all_data <- get_data("pacific cod")
  all_data = as.data.frame(all_data)

  # Covariates
  EBS = st_read( R'(C:\Users\James.Thorson\Desktop\Git\VAST\inst\region_shapefiles\EBSshelf)' )
  NBS = st_read( R'(C:\Users\James.Thorson\Desktop\Git\VAST\inst\region_shapefiles\NBS)' )
  GOA = st_read( R'(C:\Users\James.Thorson\Desktop\Git\VAST\inst\region_shapefiles\GOA)' )
  domain = st_union( st_geometry(EBS), st_geometry(NBS) )
  domain = st_union( domain, st_geometry(GOA) )
  #domain = st_transform( domain, crs = )

  # Make indicator for region
  if( FALSE ){
    sf_bathy = st_make_grid( domain, cellsize = 50*1e3 )
    intersects = st_intersects(sf_bathy, domain)
    sf_bathy = sf_bathy[ which(as.integer(intersects)==1) ]
    bathy = st_coordinates(st_centroid(sf_bathy))
    intersects_GOA = st_intersects(sf_bathy, GOA)
    bathy = data.frame(bathy, covar = ifelse(is.na(as.integer(intersects_GOA)), 0, 1) )
    colnames(bathy) = c( "X", "Y", "covar" )
    bathy[,c("X","Y")] = bathy[,c("X","Y")] / 1000
  }
  # Make bathymetry
  if( TRUE ){
    # GOA bathymetry
    load( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariates for log_tau in SPDE\GOA_bathy.rda)')
    GOA_bathy = terra::unwrap(GOA_bathy)
    GOA_bathy = terra::aggregate(GOA_bathy, fact = 5)
    GOA_bathy_df = as.data.frame(GOA_bathy, xy = TRUE)
    colnames(GOA_bathy_df) = c( "X", "Y", "covar" )
    #
    BS_bathy = terra::rast( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariates for log_tau in SPDE\Bathy.grd)' )
    BS_bathy = terra::aggregate(BS_bathy, fact = 5)
    BS_bathy_df = as.data.frame(BS_bathy, xy = TRUE)
    colnames(BS_bathy_df) = c( "X", "Y", "covar" )
    #
    bathy = rbind( GOA_bathy_df, BS_bathy_df)
    bathy[,c("X","Y")] = bathy[,c("X","Y")] / 1e5
    bathy[,"covar"] = bathy[,"covar"] / 1e3
  }

  #
  data = subset( all_data, survey_name %in% c('eastern Bering Sea', 'Gulf of Alaska', 'northern Bering Sea') )
  data = data.frame( X = data[,'lon_start'], Y = data[,'lat_start'], year = data[,'year'],
                     #covar = ifelse(data[,'survey_name'] == "Gulf of Alaska", 1, 0),
                     covar = data[,'depth_m'] / 1e3,
                     density = data[,'catch_weight'] )
  data[,c("X","Y")] = sf_project( to = st_crs(domain), from = st_crs(4326), pts = data[,c("X","Y")] ) / 1e5

  # Make mesh
  mesh <- fmesher::fm_mesh_2d(data[, c("X", "Y")], cutoff = 25 / 1e2 )  # 0.5 for Lon-Lat
  mesh_with_covs <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = "covar",
    coords  = c("X", "Y")
  )
  plot(mesh)
  points( x = mesh$loc[,1], y = mesh$loc[,2], col = viridisLite::viridis(10)[cut(mesh_with_covs$vertex[,1],breaks=10)] )

  #
  xy_i = data.frame( X = c(-224000, -1178000), Y = c(908000, 1279000)  ) / 1e5
}

# Plot covar
ggplot(bathy) +
    geom_tile( aes(X, Y, fill = covar) )  +  # , colour = "grey50"
    coord_fixed()

# Which polygon to plot
polygon_i = RANN::nn2( data = mesh$loc[,1:2], query = xy_i, k = 1 )

# Inputs
#spatial_domain = mesh_with_covs
#kappa_formula = ~ 0

#
#Xkappa_sk = model.matrix(
#  update.formula( old = kappa_formula, new = "~ . + 0" ),
#  mesh_with_covs$vertex_covariates
#)

fit1 <- tinyVAST(
  formula = density ~ 0 + s(covar),    # 0 + factor(year)
  data = data,
  space_term = "",
  #spacetime_term = "",
  #spatial_domain = mesh,
  spatial_domain = mesh_with_covs,
  family = tweedie(link = "log"),
  space_columns = c('X','Y'),
  time_column = "year",
  times = 1982:2024,
  control = tinyVASTcontrol( use_anisotropy = FALSE, run_model = TRUE, trace = 1 )
)

# Fit a Tweedie spatial random field GLMM with a smoother for covar:
fit2 <- tinyVAST(
  formula = density ~ 0 + s(covar),    # 0 + factor(year)
  data = data,
  space_term = "",
  #spacetime_term = "",
  #spatial_domain = mesh,
  spatial_domain = mesh_with_covs,
  family = tweedie(link = "log"),
  development = list(
    kappa_formula = ~ covar
  ),
  space_columns = c('X','Y'),
  time_column = "year",
  times = 1982:2024,
  control = tinyVASTcontrol( use_anisotropy = FALSE, run_model = TRUE, trace = 1 )
)



# Plot Omega
#bathy$omega = predict(fit1, what = "pomega_g", newdata = bathy )
#ggplot(bathy, aes(X, Y, fill = omega)) +
#    geom_tile() +
#    coord_fixed()

#setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)' )
#library(TMB)
#compile("tinyVAST.cpp", framework = "TMBad")
#dyn.load( dynlib("tinyVAST") )
#obj = MakeADFun( data = fit$tmb_data, par = fit$tmb_par, map = fit$tmb_map, random = fit$tmb_random,
#                 DLL = "tinyVAST" )

# Compute initial density
A_gs = fmesher::fm_evaluator( mesh, loc=as.matrix(bathy[c('X','Y')]) )$proj$A
v0_s = rep(0, mesh$n)
v0_s[polygon_i$nn.idx] = 1
bathy$v0_g = as.numeric(A_gs %*% v0_s)
bathy$v0_g = ifelse( bathy$v0_g < 1e-6, NA, bathy$v0_g )

#
#ggplot(bathy, aes(X, Y, fill = v0_g)) +
#    geom_tile()  +
#    coord_fixed()

# Compute diffused density
G1 = fit1$rep$G1
log_kappa = fit1$opt$par['log_kappa']
invD = ( Matrix::Diagonal(n = mesh$n) + exp(-2 * log_kappa) * fit1$tmb_inputs$tmb_data$spatial_list$G0_inv %*% G1 )
v1_s = Matrix::solve( invD, v0_s )
bathy$v1_g = as.numeric(A_gs %*% v1_s)
bathy$v1_g = ifelse( bathy$v1_g < 1e-6, NA, bathy$v1_g )

# Compute diffused density
G1 = fit2$rep$G1
log_kappa = fit2$opt$par['log_kappa']
invD = ( Matrix::Diagonal(n = mesh$n) + exp(-2 * log_kappa) * fit2$tmb_inputs$tmb_data$spatial_list$G0_inv %*% G1 )
v2_s = Matrix::solve( invD, v0_s )
bathy$v2_g = as.numeric(A_gs %*% v2_s)
bathy$v2_g = ifelse( bathy$v2_g < 1e-6, NA, bathy$v2_g )


#
p0 = ggplot(bathy, aes(X, Y, fill = v0_g )) +
    geom_tile()  +
    coord_fixed() + scale_fill_continuous(trans = "log")
p1 = ggplot(bathy, aes(X, Y, fill = v1_g )) +
    geom_tile()  +
    coord_fixed() + scale_fill_continuous(trans = "log")
p2 = ggplot(bathy, aes(X, Y, fill = v2_g )) +
    geom_tile()  +
    coord_fixed() + scale_fill_continuous(trans = "log")
grid.arrange( arrangeGrob(p0, p1, p2, nrow=1) )
ggsave( file = R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\scratch\aniso.png)', width = 5, height = 5,
        plot = grid.arrange(arrangeGrob(p0, p1, nrow=2)) )
