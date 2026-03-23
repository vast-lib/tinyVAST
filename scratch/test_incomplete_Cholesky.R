
# Install necessary dependencies
devtools::install_github( "vast-lib/tinyVAST@dev" )
devtools::install_github( "kaskr/adcomp/TMB@ichol" )

# Packages
library(tinyVAST)
library(TMB)
library(Matrix)
library(ggplot2)
library(igraph)
library(pracma)

# Set threads
Sys.setenv(OPENBLAS_NUM_THREADS = 1)

# 
sim = c("diffusion", "spde")[2]
sigma = 1

if( sim == "diffusion" ){
  kappa = c( 0.8, 0.9 )
  tau = 0.1
  #
  n_t <- 500
  n_s <- 500
  P_tt = bandSparse( 
    n = n_t,
    k = c(-1,0),
    diagonals = list(rep(1,n_t), rep(-1,n_t)) 
  )
  P_ss = bandSparse(
    n = n_s,
    k = c(-1,0,1),
    diagonals = list( rep(0.5,n_s), rep(-1,n_s), rep(0.5,n_s) )
  )
  I_tt = Diagonal(n_t)
  I_ss = Diagonal(n_s)
  I_kk = Diagonal(n_t*n_s)
  invV_kk = tau^2 * Diagonal( n_t * n_s )
  P_kk = kappa[2] * kronecker(P_tt, I_ss) + kappa[1] * kronecker(I_tt, P_ss)
  Q_jj = (I_kk - t(P_kk) ) %*% invV_kk %*% (I_kk - P_kk)

  #
  sigma = 1
  Data = expand.grid( site = seq_len(n_s), time = seq_len(n_t) )
  Data$z = RTMB:::rgmrf0(n = 1, Q = Q_jj)
  Data$n = rnorm(n = nrow(Data), mean = Data$z, sd = sigma)
  Data$c = rpois(n = nrow(Data), lambda = exp(Data$z) )

  # 
  edges <- cbind(seq_len(n_s-1), seq_len(n_s-1)+1)
  mesh <- graph_from_edgelist(edges, directed = FALSE)
  V(mesh)$name = seq_len(n_s)
  
  #
  space_columns = "site"
  space_term = NULL
  spacetime_term = "response -> response, 1, ar"
}
if( sim == "spde" ){
  n_j = 1e5
  
  # Fixed sample density ... 100 samples per unit area
  #Data = poisson2disk( n=n_j, a = sqrt(n_j / 100), b = sqrt(n_j / 100) )
  Data = data.frame(x = runif(n_j), y = runif(n_j) ) * sqrt( n_j / 100 )
  #Data = expand.grid( x = seq(0,1,length = 10), y = seq(0,1, length = 10) )
  
  # make SPDE mesh for spatial term
  mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], cutoff = 0.02 )
  spde = fmesher::fm_fem( mesh )
  A_js = fmesher::fm_evaluator( mesh, as.matrix(Data) )$proj$A
  
  # Simulate
  range = 0.5              # range = sqrt(8) / kappa
  SD = 1                   # Var = 1 / 4*pi / kappa^2 / tau^2
  kappa = sqrt(8) / range
  tau = 1 / sqrt(4*pi) / kappa / SD
  Q_ss = tau^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  z_s = RTMB:::rgmrf0( n = 1, Q_ss )
  
  # Simulate nuissance parameter z from oscillatory (day-night) process
  Data$z = as.numeric(A_js %*% z_s)
  Data$n = Data$z + rnorm(nrow(Data), sd= sigma)
  Data$c = rpois( n = nrow(Data), lambda = exp(Data$z) )
  Data$time = 1
  
  ggplot(Data, aes(x, y)) +
    stat_summary_hex(
      aes(z = z, fill = after_stat(value)),
      fun = mean,
      bins = 20
    ) +
    coord_fixed() +
    scale_fill_viridis_c(name = "Mean density") +
    theme_minimal()
  
  # Settings
  space_columns = c( "x", "y" )
  space_term = ""
  spacetime_term = NULL
}

# fit with incomplete decomposition
out1 = tinyVAST( 
  data = Data,
  formula = c ~ 1,
  family = poisson(),
  spatial_domain = mesh,
  space_columns = space_columns,
  time_column = "time",
  development = list(
    method = "incomplete", 
    #abstol = 0,
    #tol = 0,
    #maxit = 10, 
    trace = FALSE
  ),
  space_term = space_term,
  spacetime_term = spacetime_term 
)
                
# complete decomposition
out2 = tinyVAST( 
  data = Data,
  formula = c ~ 1, 
  family = poisson(),
  spatial_domain = mesh,
  space_columns = space_columns,
  time_column = "time",
  space_term = space_term,
  spacetime_term = spacetime_term 
)

out1$run
out2$run                
cbind(out1$opt$par, out2$opt$par)
c(out1$opt$obj, out2$opt$obj)
