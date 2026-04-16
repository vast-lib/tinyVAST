

# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

# Install tinyVAST@dev branch
#remotes::install_github( R'(vast-lib/tinyVAST@dev)', force = TRUE )

setwd( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariate anisotropy SPDE)' )

library(tinyVAST)
library(fmesher)
library(ggplot2)
library(Matrix)
library(sf)

# Euclidean and polar coordinates
euclidean_iz = pracma::poisson2disk( n = 1000 ) - 0.5
polar_iz = cbind( d = sqrt(rowSums(euclidean_iz^2)),
                  theta = atan2(euclidean_iz[,2],euclidean_iz[,1]) )

# Simulate preference
data_i = data.frame( "X" = euclidean_iz[,1], "Y" = euclidean_iz[,2],
                     polar_iz, "p" = -1 * abs(polar_iz[,'d'] - 0.3) / 0.2 )

# Make plotting grid
domain = st_polygon( list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0,0)) - 0.5) )
sf_grid = st_make_grid( domain, cellsize = 0.01 )

# Make counter lines
if( FALSE ){
  # Step 1
  grid_sf <- st_sf(id = 1:length(sf_grid), geometry = sf_grid)
  # Step 2: Convert to terra::SpatVector
  grid_vect <- vect(grid_sf)
  # Step 3: Create raster template with same extent/resolution
  r_template <- rast(grid_vect, resolution = 0.01)
  # Step 4: Rasterize based on ID or other attribute
  r <- rasterize(grid_vect, r_template, field = "id")
  # Step 1: Create example raster
  values(r) <- matrix(outer(1:100, 1:100, function(x, y) sin(x/10) + cos(y/10)), 100, 100)
  # Step 2: Extract contour lines (as sf LINESTRINGs)
  contours <- as.contour(r, levels = c(0.5, 1.0, 1.5)) |> st_as_sf()
}

# Show preference
ggplot( data_i ) + geom_point( aes(x=X, y=Y, col=d) )
ggplot( data_i ) + geom_point( aes(x=X, y=Y, col=p) )

# Make mesh
mesh = fm_mesh_2d(
  loc = euclidean_iz,
  cutoff = 0
)
spde = fm_fem( mesh )
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(sf_grid)) )$proj$A
A_is = fm_evaluator( mesh, loc=as.matrix(data_i[,c("X","Y")]) )$proj$A

# Visualize mesh
plot(mesh)
points( euclidean_iz )

#make_stiffness = tinyVAST:::make_stiffness
#source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\internal.R)')
make_stiffness = tinyVAST:::make_stiffness
mesh = add_mesh_covariates(
  mesh = mesh,
  data = data_i,
  covariates = c("d","theta"),
  coords  = c("X", "Y")
)
mesh$triangle_covariates = matrix( 0, nrow = nrow(mesh$graph$tv), ncol = 0 )
M1 = make_stiffness(
  mesh,
  loc = cbind( mesh$loc[,1:2], as.matrix(mesh$vertex) ),
  H = diag( 10 ^ (c(1, 1, -1, 1)) )
  #loc = mesh$loc[,1:2],
  #H = diag( 10 ^ c(2,2) )
)
kappa = 100
M0 = spde$c0

# Make indicator to show diffusion
nn = RANN::nn2( data = mesh$loc[,1:2], query = matrix(c(0.3,0),nrow=1), k = 1 )
v0_s = v = rep(0, mesh$n)
v0_s[nn$nn.idx[1]] = 1

# Plot diffusion
invD = Diagonal( n = mesh$n ) + 1/kappa^2 * solve(M0) %*% M1
v1_s = solve(invD, v0_s)
sf_plot = st_sf( sf_grid, v = as.numeric(A_gs %*% v1_s) )
plot(sf_plot, border = NA)

tau = 0.0001 # for  H = diag( 10 ^ (c(1, 1, -1, 1)) )
#tau = 10 # for H = diag( 10 ^ c(0,0,3,3) )
M2 = M1 %*% solve(M0) %*% M1
Q = tau^2 * ( kappa^4 * M0 + 2*kappa^2*M1 + M2 )
w_s = rmvnorm_prec( prec = Q )[,1]
data_i$w_i = as.numeric(A_is %*% w_s)
data_i$mu_i = data_i$w_i + 5*data_i$p
data_i$y_i = rnorm( n = nrow(data_i), mean = data_i$mu_i, sd = 1 * mean(sqrt(diag(solve(Q)))) )

# Plot correlations
#source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\internal.R)' )
#r0_s = cov2cor(solve(Q))[,nn$nn.idx[1]]
#r0_g = (A_gs %*% r0_s)[,1]
pred = st_coordinates(st_centroid(sf_grid))
r_g = spatial_cor( Q, mesh, c(0.3,0), pred )
sf_plot = st_sf( sf_grid, r = r_g ) # r0 = r0_g,
plot(sf_plot, border = NA)

# Plot field
sf_plot = st_sf( sf_grid, v = as.numeric(A_gs %*% w_s)  )
plot(sf_plot, border = NA)

#
ggplot( data_i ) + geom_point( aes( x=X, y=Y, col = y_i) )
tmp = data_i
tmp$X = cut(data_i$X, breaks = 20)
tmp$Y = cut(data_i$Y, breaks = 20)
ggplot( tmp, aes( x=X, y=Y, fill = y_i) ) + geom_tile( )

# Fit with covariate-anisotropy
fit = tinyVAST(
  data = data_i,
  formula = y_i ~ 1 + poly(d, 2),
  #formula = y_i ~ 1,
  spatial_domain = mesh,
  space_term = "",
  space_columns = c("X","Y"),
  development = list(
    vertex_formula = ~ d
  ),
  control = tinyVASTcontrol( trace = 1, use_anisotropy = FALSE)
)

# Plot diffusion
#invD = Diagonal( n = mesh$n ) + exp(-2 * fit$opt$par['log_kappa']) * solve(M0) %*% fit$rep$G1
#v1_s = solve(invD, v0_s)
cor_g = spatial_cor( fit$rep$Q, mesh, c(0.3,0), pred )
sf_plot = st_sf( sf_grid, v = cor_g )
plot(sf_plot, border = NA)

# Fit without covariate-anisotropy
fit0 = tinyVAST(
  data = data_i,
  formula = y_i ~ 1 + poly(d, 2),
  #formula = y_i ~ 1,
  spatial_domain = mesh,
  space_term = "",
  space_columns = c("X","Y"),
  control = tinyVASTcontrol( trace = 1, use_anisotropy = TRUE)
)

# Plot correlation
cor_g = spatial_cor( fit0$rep$Q, mesh, c(0.3,0), pred )
sf_plot = st_sf( sf_grid, v = cor_g )
plot(sf_plot, border = NA)

# Compare performance
performance = rbind(
  AIC = c( "covar" = AIC(fit), "null" = AIC(fit0)),
  cAIC = c( cAIC(fit), cAIC(fit0) ),
  CV = c( cv::cv(fit)[['CV crit']], cv::cv(fit0)[['CV crit']] )
)
performance

