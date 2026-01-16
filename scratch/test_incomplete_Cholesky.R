
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )
library(tinyVAST)
library(TMB)
library(Matrix)

# 
n_j = 1e4

# Fixed sample density ... 100 samples per unit area
Data = data.frame(x = runif(n_j), y = runif(n_j) ) * sqrt( n_j / 100 )

# make SPDE mesh for spatial term
mesh = fmesher::fm_mesh_2d( Data[,c('x','y')] )
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
Data$n = Data$z + rnorm(nrow(Data), sd=1)

# Add columns for multivariate and/or temporal dimensions
Data$var = "n"

# fit model with cyclic confounder as GAM term
out1 = tinyVAST( data = Data,
                formula = n ~ 1, 
                spatial_domain = mesh,
                development = list(method = "incomplete", trace = FALSE),
                space_term = "n <-> n, sd_n" )
                
out2 = tinyVAST( data = Data,
                formula = n ~ 1, 
                spatial_domain = mesh,
                space_term = "n <-> n, sd_n" )

out1$run
out2$run                
cbind(out1$opt$par, out2$opt$par)
c(out1$opt$obj, out2$opt$obj)
