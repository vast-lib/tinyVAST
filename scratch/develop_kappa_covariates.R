
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

library(tinyVAST)

#library(sf)
data( pcod, package = "sdmTMB" )
data( qcs_grid, package = "sdmTMB" )
data( pcod_2011, package = "sdmTMB" )

# Piped version
mesh <- fmesher::fm_mesh_2d(pcod[, c("X", "Y")], cutoff = 5)
mesh_with_covs <- add_vertex_covariates(
  mesh,
  qcs_grid,
  covariates = c("depth_scaled", "depth_scaled2"),
  coords  = c("X", "Y")
)

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
  formula = density ~ s(depth),
  data = as.data.frame(pcod_2011),
  space_term = "",
  spatial_domain = mesh,
  #spatial_domain = mesh_with_covs,
  family = tweedie(link = "log"),
  development = list(
  #  kappa_formula = ~ depth_scaled
  ),
  space_columns = c('X','Y'),
  control = tinyVASTcontrol( use_anisotropy = TRUE, run_model = TRUE )
)

#
qcs_grid$omega = predict(fit, what = "pomega_g", newdata = qcs_grid )

library(ggplot2)
ggplot(qcs_grid, aes(X, Y, fill = omega)) +
    geom_raster()  +
    coord_fixed()

#setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)' )
#library(TMB)
#compile("tinyVAST.cpp", framework = "TMBad")
#dyn.load( dynlib("tinyVAST") )
#obj = MakeADFun( data = fit$tmb_data, par = fit$tmb_par, map = fit$tmb_map, random = fit$tmb_random,
#                 DLL = "tinyVAST" )

library(ggplot2)
ggplot(qcs_grid, aes(X, Y, fill = depth)) +
    geom_raster()  +
    coord_fixed()

#
xy_i = c( 500, 5685 )
polygon_i = RANN::nn2( data = mesh$loc[,1:2], query = matrix(xy_i,ncol=2), k = 1 )

#
A_gs = fmesher::fm_evaluator( mesh, loc=as.matrix(qcs_grid[c('X','Y')]) )$proj$A
v0_s = rep(0, mesh$n)
v0_s[polygon_i$nn.idx] = 1
v0_g = as.numeric(A_gs %*% v0_s)
qcs_grid$v0_g = v0_g

#
ggplot(qcs_grid, aes(X, Y, fill = v0_g)) +
    geom_raster()  +
    coord_fixed()

# Diffusion
G1 = fit$rep$G1
kappa = exp(fit$opt$par['log_kappa'])
invD = ( Matrix::Diagonal(n = mesh$n) + kappa^2 * fit$tmb_inputs$tmb_data$spatial_list$G0_inv %*% G1 )
v1_s = Matrix::solve( invD, v0_s )
v1_g = as.numeric(A_gs %*% v1_s)

#
ggplot(qcs_grid, aes(X, Y, fill = v1_g)) +
    geom_raster()  +
    coord_fixed()

