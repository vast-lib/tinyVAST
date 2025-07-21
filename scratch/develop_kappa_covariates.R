
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

library(tinyVAST)
library(ggplot2)
library(gridExtra)

#library(sf)
if( FALSE ){
  data( pcod, package = "sdmTMB" )
  data( qcs_grid, package = "sdmTMB" )
  data( pcod_2011, package = "sdmTMB" )

  # Rename
  bathy = qcs_grid
  data = as.data.frame(pcod_2011)

  # Piped version
  mesh <- fmesher::fm_mesh_2d(pcod[, c("X", "Y")], cutoff = 5)
  mesh_with_covs <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = c("depth_scaled", "depth_scaled2"),
    coords  = c("X", "Y")
  )
  xy_i = matrix( c( 500, 5685 ), nrow = 1 )
}
if( TRUE ){
  # pak::pkg_install("DFO-NOAA-Pacific/surveyjoin")
  library(surveyjoin)
  library(sf)
  #cache_data()
  load_sql_data()
  all_data <- get_data("pacific cod")
  all_data = as.data.frame(all_data)

  # Covariates
  load( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariates for log_tau in SPDE\GOA_bathy.rda)')
  GOA_bathy = terra::unwrap(GOA_bathy)

  # Subset
  #data = subset( all_data, survey_name %in% c('eastern Bering Sea', 'Gulf of Alaska', 'northern Bering Sea') )
  data = subset( all_data, survey_name %in% c('Gulf of Alaska') )
  data = data.frame( X = data[,'lon_start'], Y = data[,'lat_start'], year = data[,'year'],
                     density = data[,'catch_weight'],
                     depth = data[,'depth_m'] )
  data[,c("X","Y")] = sf_project( to = st_crs(GOA_bathy), from = st_crs(4326), pts = data[,c("X","Y")] ) / 1000

  # Downsample bathymetry
  GOA_bathy = terra::aggregate(GOA_bathy, fact = 5)

  #
  bathy = as.data.frame(GOA_bathy, xy = TRUE)
  #bathy = sf_project( from = st_crs(GOA_bathy), to = st_crs(4326), pts = bathy )
  #bathy = as.data.frame(bathy)
  colnames(bathy) = c( "X", "Y", "depth" )
  bathy[,c("X","Y")] = bathy[,c("X","Y")] / 1000

  # Make above-land very different
  bathy$depth = ifelse( bathy$depth < 0, bathy$depth * 1000, bathy$depth )

  # Make mesh
  mesh <- fmesher::fm_mesh_2d(data[, c("X", "Y")], cutoff = 25 )  # 0.5 for Lon-Lat
  mesh_with_covs <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = "depth",
    coords  = c("X", "Y")
  )
  terra::plot( GOA_bathy )
  xy_i = data.frame( -88863, 785481  ) / 1000
  #xy_i = sf_project( from = st_crs(GOA_bathy), to = st_crs(4326), pts = xy_i )
}

# Plot depth
ggplot(bathy, aes(X, Y, fill = depth)) +
    geom_tile()  +
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

# Fit a Tweedie spatial random field GLMM with a smoother for depth:
fit <- tinyVAST(
  formula = density ~ s(depth) + 0 + factor(year),
  data = data,
  space_term = "",
  spacetime_term = "",
  #spatial_domain = mesh,
  spatial_domain = mesh_with_covs,
  family = tweedie(link = "log"),
  development = list(
    kappa_formula = ~ depth
  ),
  space_columns = c('X','Y'),
  time_column = "year",
  times = 1990:2023,
  control = tinyVASTcontrol( use_anisotropy = FALSE, run_model = TRUE, trace = 1 )
)


# Plot Omega
bathy$omega = predict(fit, what = "pomega_g", newdata = bathy )
library(ggplot2)
ggplot(bathy, aes(X, Y, fill = depth, col = omega)) +
    geom_tile() +
    coord_fixed()

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
G1 = fit$rep$G1
log_kappa = fit$opt$par['log_kappa']
invD = ( Matrix::Diagonal(n = mesh$n) + exp(-2 * log_kappa) * fit$tmb_inputs$tmb_data$spatial_list$G0_inv %*% G1 )
v1_s = Matrix::solve( invD, v0_s )
bathy$v1_g = as.numeric(A_gs %*% v1_s)
bathy$v1_g = ifelse( bathy$v1_g < 1e-6, NA, bathy$v1_g )

#
p0 = ggplot(bathy, aes(X, Y, fill = v0_g )) +
    geom_tile()  +
    coord_fixed() + scale_fill_continuous(trans = "log")
p1 = ggplot(bathy, aes(X, Y, fill = v1_g )) +
    geom_tile()  +
    coord_fixed() + scale_fill_continuous(trans = "log")
grid.arrange( arrangeGrob(p0, p1, nrow=2) )
ggsave( file = R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\scratch\aniso.png)', width = 5, height = 5,
        plot = grid.arrange(arrangeGrob(p0, p1, nrow=2)) )
