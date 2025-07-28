
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

library(tinyVAST)
library(ggplot2)
library(gridExtra)

# pak::pkg_install("DFO-NOAA-Pacific/surveyjoin")
library(surveyjoin)
library(sf)
library(rnaturalearth)
#cache_data()

case = c( "pcod_2011", "GOA_pcod", "NorPac_pcod" )[3]

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
  mesh <- fmesher::fm_mesh_2d(pcod_2011[, c("X", "Y")], cutoff = 2)
  mesh_with_covs <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = "covar",
    coords  = c("X", "Y")
  )

  #
  #plot( x = bathy$X, y = bathy$Y, col = viridisLite::viridis(10)[cut(bathy$covar,breaks=10)], cex = 1, pch = 20 )
  #xy_i = bathy[ which.max(bathy$covar), c("X","Y") ]
  xy_i = data.frame( 415, 5783 )
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
  BS = st_union( st_geometry(EBS), st_geometry(NBS) )
  GOA = st_geometry(GOA)
  domain = st_union( BS, st_geometry(GOA) )
  #domain = st_transform( domain, crs = )

  # map
  AK = ne_states(country = "united states of america") # , regions = "alaska" )
  AK = AK[pmatch("Alas", AK$name_en),]
  AK = st_geometry(AK)
  AK = st_transform( AK, crs = st_crs(domain) )
  
  # Make indicator for region
  if( TRUE ){
    # Set up grid ... keeping excess before making covariates
    sf_bathy = st_make_grid( domain, cellsize = 25*1e3 )
    intersects = st_intersects(sf_bathy, domain)
    sf_plot = sf_bathy[ which(as.integer(intersects)==1) ]
    bathy = st_coordinates(st_centroid(sf_bathy))
    # Add GOA
    intersects_GOA = st_intersects(sf_bathy, GOA)
    bathy = data.frame(bathy, GOA = ifelse(is.na(as.integer(intersects_GOA)), 0, 1) )
    # Add EBS
    intersects_BS = st_intersects(sf_bathy, BS)
    bathy = data.frame(bathy, BS = ifelse(is.na(as.integer(intersects_BS)), 0, 1) )
    # Add land
    intersects_AK = st_intersects(sf_bathy, AK)
    bathy = data.frame(bathy, AK = ifelse(is.na(as.integer(intersects_AK)), 0, 1) )
    # Change names
    #colnames(bathy) = c( "X", "Y", "covar" )
    bathy[,c("X","Y")] = bathy[,c("X","Y")] / 1000
  }
  # Make bathymetry
  if( TRUE ){
    # GOA bathymetry
    load( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariate anisotropy SPDE\GOA_bathy.rda)')
    GOA_bathy = terra::unwrap(GOA_bathy)
    GOA_bathy = terra::aggregate(GOA_bathy, fact = 5)
    GOA_bathy_df = as.data.frame(GOA_bathy, xy = TRUE)
    colnames(GOA_bathy_df) = c( "X", "Y", "depth" )
    #
    BS_bathy = terra::rast( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariate anisotropy SPDE\Bathy.grd)' )
    BS_bathy = terra::aggregate(BS_bathy, fact = 5)
    BS_bathy_df = as.data.frame(BS_bathy, xy = TRUE)
    colnames(BS_bathy_df) = c( "X", "Y", "depth" )
    #
    combo_bathy = rbind( GOA_bathy_df, BS_bathy_df)
    nni = RANN::nn2( data = combo_bathy[,c("X","Y")] / 1000, 
                     query = bathy[,c("X","Y")], k = 1 )$nn.idx[,1]
    bathy$depth = combo_bathy[nni,'depth']
  }

  #
  data = subset( all_data, survey_name %in% c('eastern Bering Sea', 'Gulf of Alaska', 'northern Bering Sea') )
  data = data.frame( X = data[,'lon_start'], Y = data[,'lat_start'], year = data[,'year'],
                     #covar = ifelse(data[,'survey_name'] == "Gulf of Alaska", 1, 0),
                     #covar = data[,'depth_m'] / 1e3,
                     density = data[,'catch_weight'] )
  data[,c("X","Y")] = sf_project( to = st_crs(domain), from = st_crs(4326), pts = data[,c("X","Y")] ) / 1e3

  # Make mesh
  mesh = fmesher::fm_mesh_2d(data[, c("X", "Y")], cutoff = 25 )  # 0.5 for Lon-Lat
  mesh_cov <- add_vertex_covariates(
    mesh,
    bathy,
    covariates = c("depth"),
    coords  = c("X", "Y")
  )
  n_tri <- length(mesh$graph$tv[, 1]) 
  posTri <- matrix(0, n_tri, 2)
  for (t in 1:n_tri) {
    temp <- mesh$loc[mesh$graph$tv[t, ], ]
    posTri[t, ] <- colMeans(temp)[c(1, 2)]
  }
  nni = RANN::nn2( data = bathy[,c("X","Y")], query = posTri, k = 1 )$nn.idx[,1]
  mesh_cov$triangle_covariates = data.frame(
    bathy[nni,c("GOA","BS","AK")]
  )
  #mesh_cov$vertex_covariates = data.frame( matrix(nrow=mesh$n,ncol=0) )
  class(mesh_cov) <- c("vertex_coords", class(mesh))
  plot(mesh)
  points( x = mesh$loc[,1], y = mesh$loc[,2], 
          col = viridisLite::viridis(10)[cut(mesh_cov$vertex[,1],breaks=10)] )
  points( x = posTri[,1], y = posTri[,2], 
          col = viridisLite::viridis(10)[cut(mesh_cov$triangle[,3],breaks=10)] )
  
  #
  xy_i = data.frame( X = c(0, -1000000), Y = c(908000, 1279000)  ) / 1e3
}

#
plot( x = bathy$X, y = bathy$Y, col = viridis() )

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

fit0 <- tinyVAST(
  formula = density ~ 1,    # 0 + factor(year)
  data = data,
  space_term = "",
  #spacetime_term = "",
  #spatial_domain = mesh,
  spatial_domain = mesh_cov,
  family = tweedie(link = "log"),
  space_columns = c('X','Y'),
  time_column = "year",
  times = 1982:2024,
  control = tinyVASTcontrol( use_anisotropy = FALSE, run_model = TRUE, trace = 1 )
)

# Fit a Tweedie spatial random field GLMM with a smoother for covar:
fit1 <- tinyVAST(
  formula = density ~ 1,    # 0 + factor(year)
  data = data,
  space_term = "",
  #spacetime_term = "",
  #spatial_domain = mesh,
  spatial_domain = mesh_cov,
  family = tweedie(link = "log"),
  development = list(
    triangle_formula = ~ offset(AK) + GOA + BS,
    vertex_formula = ~ depth
  ),
  space_columns = c('X','Y'),
  time_column = "year",
  times = 1982:2024,
  control = tinyVASTcontrol( use_anisotropy = FALSE, run_model = TRUE, trace = 1 )
)

# Plot correlations
intersects = st_intersects(sf_bathy, domain)
sf_plot = sf_bathy[ which(as.integer(intersects)==1) ]
pred = st_coordinates(st_centroid(sf_plot)) / 1000
r1 = spatial_cor( Q = fit1$rep$Q, mesh = mesh, coord = as.numeric(xy_i[1,]), pred = pred )
r2 = spatial_cor( Q = fit1$rep$Q, mesh = mesh, coord = as.numeric(xy_i[2,]), pred = pred )
sf_plot = st_sf( sf_plot, r1 = r1, r2 = r2 )
plot( sf_plot, border=NA )

# Compare performance
performance = rbind(
  AIC = c( "null" = AIC(fit0), "covar" = AIC(fit1)),
  cAIC = c( cAIC(fit0), cAIC(fit1) ),
  CV = c( cv::cv(fit0)[['CV crit']], cv::cv(fit1)[['CV crit']] )
)


# Plot Omega
#bathy$omega = predict(fit0, what = "pomega_g", newdata = bathy )
#ggplot(bathy, aes(X, Y, fill = omega)) +
#    geom_tile() +
#    coord_fixed()

#setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)' )
#library(TMB)
#compile("tinyVAST.cpp", framework = "TMBad")
#dyn.load( dynlib("tinyVAST") )
#obj = MakeADFun( data = fit$tmb_data, par = fit$tmb_par, map = fit$tmb_map, random = fit$tmb_random,
#                 DLL = "tinyVAST" )

#
trun = \(v, threshold = 0.95){
  v = ifelse( is.na(v), 0, v )
  ord = order(v, decreasing = TRUE )
  cum = cumsum( v[ord] )
  above = ( cum < (max(cum) * threshold) )
  index = which.max( which(above) )
  v = ifelse( seq_along(v) %in% ord[seq_len(index)], v, NA )
  return(v)
}

# Compute initial density
A_gs = fmesher::fm_evaluator( mesh, loc=as.matrix(bathy[c('X','Y')]) )$proj$A
v0_s = rep(0, mesh$n)
v0_s[polygon_i$nn.idx] = 1
bathy$v0_g = as.numeric(A_gs %*% v0_s)
#bathy$v0_g = trun( bathy$v0_g, threshold = 0.95 )
bathy$v0_g = ifelse( bathy$v0_g > 0.1, bathy$v0_g, NA )

#
#ggplot(bathy, aes(X, Y, fill = v0_g)) +
#    geom_tile()  +
#    coord_fixed()

# Compute diffused density   /  correlation
#G1 = fit0$rep$G1
#log_kappa = fit0$opt$par['log_kappa']
#invD = ( Matrix::Diagonal(n = mesh$n) + exp(-2 * log_kappa) * fit0$tmb_inputs$tmb_data$spatial_list$G0_inv %*% G1 )
v1_s = Matrix::solve( fit0$rep$Q_ss, v0_s )
bathy$v1_g = as.numeric(A_gs %*% v1_s / max(v1_s))
#bathy$v1_g = trun( bathy$v1_g, threshold = 0.95 )
bathy$v1_g = ifelse( bathy$v1_g > 0.1, bathy$v1_g, NA )

# Compute diffused density   /  correlation
#G1 = fit1$rep$G1
#log_kappa = fit1$opt$par['log_kappa']
#invD = ( Matrix::Diagonal(n = mesh$n) + exp(-2 * log_kappa) * fit1$tmb_inputs$tmb_data$spatial_list$G0_inv %*% G1 )
v2_s = Matrix::solve( fit1$rep$Q_ss, v0_s )
bathy$v2_g = as.numeric(A_gs %*% v2_s / max(v2_s))
#bathy$v2_g = trun( bathy$v2_g, threshold = 0.95 )
bathy$v2_g = ifelse( bathy$v2_g > 0.1, bathy$v2_g, NA )

#
pbathy = ggplot(bathy, aes(X, Y, fill = covar )) +
    geom_tile()  +
    coord_fixed()
pinit = ggplot(bathy, aes(X, Y, fill = v0_g )) +
    geom_tile()  +
    coord_fixed() + #scale_fill_continuous(trans = "log") +
    scale_x_continuous( limits = xy_i[[1]] + c(-1,1) * 25 ) +
    scale_y_continuous( limits = xy_i[[2]] + c(-1,1) * 25 )
p0 = ggplot(bathy, aes(X, Y, fill = v1_g )) +
    geom_tile()  +
    coord_fixed() + #scale_fill_continuous(trans = "log") +
    scale_x_continuous( limits = xy_i[[1]] + c(-1,1) * 25 ) +
    scale_y_continuous( limits = xy_i[[2]] + c(-1,1) * 25 )
p1 = ggplot(bathy, aes(X, Y, fill = v2_g )) +
    geom_tile()  +
    coord_fixed() + #scale_fill_continuous(trans = "log") +
    scale_x_continuous( limits = xy_i[[1]] + c(-1,1) * 25 ) +
    scale_y_continuous( limits = xy_i[[2]] + c(-1,1) * 25 )
pfull = grid.arrange( arrangeGrob(pinit, pbathy, p0, p1, nrow=2) )
ggsave( file = R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\scratch\aniso.png)', width = 8, height = 8,
        plot = pfull )
