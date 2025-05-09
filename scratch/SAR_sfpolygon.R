
rmvnorm_prec <-
function( mu, # estimated fixed and random effects
          prec, # estimated joint precision
          n.sims) {

  require(Matrix)
  # Simulate values
  z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  # Q = t(P) * L * t(L) * P
  L = Matrix::Cholesky(prec, super=TRUE)
  # Calcualte t(P) * solve(t(L)) * z0 in two steps
  z = Matrix::solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = Matrix::solve(L, z, system = "Pt") # z = Pt    * z
  return(mu + as.matrix(z))
}

# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force=TRUE )
library(tinyVAST)
#library(sf)
# library(Matrix)

#data(bering_sea)
domain = bering_sea
#domain = st_polygon( x = list( cbind(x=c(-1,1,1,-1,-1),y=c(-1,-1,1,1,-1)) ) )

sf_grid = sf::st_make_grid(
  domain,
  cellsize = 1 * c(1,1),
  square = FALSE
)
#sf_grid = st_intersection( sf_grid, bering_sea )
# inherits( grid, "sfc_GEOMETRY" )

# Settings:
# Queen or Rook
st_rook <- function(m, ...) sf::st_relate(m, m, pattern="F***1****", ...)
st_queen <- function(m, ...) sf::st_relate(m, m, pattern = "F***T****", ...)
grid_A  <- st_queen(sf_grid, sparse=TRUE)
A_ss       <- as(grid_A, "sparseMatrix")

grid_xy = sf::st_coordinates( sf::st_centroid(sf_grid))

plot(domain)
plot(sf_grid, add=TRUE)
text( grid_xy, labels = seq_along(sf_grid) )

tripA = Matrix::mat2triplet(A_ss)
#dx_ss = sparseMatrix(i = tripA$i, j = tripA$j, x = abs(grid_xy[tripA$i,1] - grid_xy[tripA$j,1]) )
#dy_ss = sparseMatrix(i = tripA$i, j = tripA$j, x = abs(grid_xy[tripA$i,2] - grid_xy[tripA$j,2]) )
i_z = tripA$i
j_z = tripA$j
delta_z2 = grid_xy[tripA$i,] - grid_xy[tripA$j,]

# Get projection matrix
samples_i = sf::st_sample(domain, size = 10)
#sf_samples_i = st_as_sf( domain, size = 10 )
s_i = as.integer(sf::st_within( samples_i, sf_grid ))
A_is = Matrix::sparseMatrix( i = seq_along(samples_i),
                     j = s_i,
                     x = 1,
                     dims = c(length(samples_i),length(sf_grid)) )

ln_H_input = c(0.5, -0.5)
  H = matrix( NA, 2, 2 )
  H[1,1] = exp(ln_H_input[1]);
  H[2,1] = ln_H_input[2];
  H[1,2] = ln_H_input[2];
  H[2,2] = (1+ln_H_input[2]*ln_H_input[2]) / exp(ln_H_input[1]);
#theta = 1 * pi / 3
#k = 3
#H = diag( c(k,1/k) ) %*% matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2, byrow=TRUE )
rho = 0.95
#kappa = 1
#rho = exp(-kappa)

d_z = ((delta_z2 %*% H)^2 %*% matrix(1,nrow=2,ncol=1))[,1]   # use sqrt ?
W_ss = Matrix::sparseMatrix(
  i = i_z,
  j = j_z,
  x = exp( -1 * d_z )
)

# Row normalize
W_ss = Matrix::Diagonal( n = nrow(W_ss), x = 1/Matrix::rowSums(W_ss) ) %*% W_ss

#
IminusW_ss = Matrix::Diagonal(n = nrow(W_ss)) - rho * W_ss
Q_ss = Matrix::t(IminusW_ss) %*% IminusW_ss
V_ss = solve(Q_ss)
R_ss = cov2cor(V_ss)

x_st = rmvnorm_prec(
  mu = rep(0,nrow(Q_ss)),
  prec = Q_ss,
  n.sims = 9
)

plot_sf = sf::st_sf( sf_grid, x_st )
plot( plot_sf )

#######################
# Fit in tinyVAST
#######################

#
mu_st = exp( 2 + x_st )
data = expand.grid(
  s = seq_len(nrow(x_st)),
  t = seq_len(ncol(x_st))
)
data$y = rpois( n = prod(dim(x_st)), lambda = mu_st )
data = cbind( data, grid_xy[data$s,] )
data$var = "logd"

#library(checkmate)
#library(TMB)
#source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R)','internal.R'))
#source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R)','fit.R'))

#
fit = tinyVAST(
  formula = y ~ 1,
  #space_term = "",
  spacetime_term = "",
  data = data,
  spatial_domain = sf_grid,
  time_column = "t",
  variable_column = "var",
  distribution_column = "var",
  family = list( "logd" = poisson() ),
  space_columns = c("X","Y"),
  control = tinyVASTcontrol( use_anisotropy = TRUE )
)

predict( fit,
         newdata = data[1:10,] )

# Global for inputs
formula = y ~ 1
#data
time_term = NULL
space_term = ""
spacetime_term = ""
family = list( "logd" = poisson() )
space_columns = c("X","Y")
spatial_domain = sf_grid
time_column = "t"
times = NULL
variable_column = "var"
variables = NULL
distribution_column = "var"
delta_options = list(formula = ~ 1)
spatial_varying = NULL
weights = NULL
control = tinyVASTcontrol()

control$sar_adjacency = "queen"
library(checkmate)

