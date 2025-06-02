
#################
# Simulate bivariate VAR
#################

library(tinyVAST)
library(fmesher)

# Simulate settings
theta_xy = 0.2
n_x = n_y = 10
n_t = 20
B = rbind( c( 0.5, -0.25),
           c(-0.1,  0.50) )

# Simulate GMRFs
R = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
d1 = mvtnorm::rmvnorm(n_t, sigma=0.2*kronecker(R,R) )
d2 = mvtnorm::rmvnorm(n_t, sigma=0.2*kronecker(R,R) )
d = abind::abind( d1, d2, along=3 )

# Project through time and add mean
for( t in seq_len(n_t) ){
  if(t>1) d[t,,] = t(B%*%t(d[t-1,,])) + d[t,,]
}

# Shape into longform data-frame and add error
Data = data.frame( expand.grid(time=1:n_t, x=1:n_x, y=1:n_y, "var"=c("d1","d2")),
                   mu = exp(as.vector(d)))
Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$mu, phi=0.5, power=1.5 )

# make mesh
mesh = fm_mesh_2d( Data[,c('x','y')] )

# Define DSEM
dsem = "
  d1 -> d1, 1, b11
  d2 -> d2, 1, b22
  d2 -> d1, 1, b21
  d1 -> d2, 1, b12
  d1 <-> d1, 0, var1
  d2 <-> d2, 0, var1
"

# fit model
out = tinyVAST( spacetime_term = dsem,
           data = Data,
           formula = n ~ 0 + var,
           spatial_domain = mesh,
           family = tweedie() )
out

#object = out
#extra_times = 21:25
newdata = Data[1:10,]
newdata$t = 25

project(
  out,
  extra_times = 21:25,
  newdata = newdata
)

######################
# Univariate version
######################

library(tinyVAST)
library(fmesher)

# Simulate settings
theta_xy = 0.4
n_x = n_y = 10
n_t = 20
t_obs = 1:8
rho = 0.8
spacetime_sd = 0.5
space_sd = 0.5
gamma = 0

# Simulate GMRFs
R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
R_ss = kronecker(R_s, R_s)
Vspacetime_ss = spacetime_sd^2 * R_ss
Vspace_ss = space_sd^2 * R_ss

# make spacetime AR1 over time
eps_ts = mvtnorm::rmvnorm( n_t, sigma=Vspacetime_ss )
for( t in seq_len(n_t) ){
  if(t>1) eps_ts[t,] = rho*eps_ts[t-1,] + eps_ts[t,]/(1 + rho^2)
}

# make space term
omega_s = mvtnorm::rmvnorm( 1, sigma=Vspace_ss )[1,]

# linear predictor
p_ts = gamma + outer( rep(1,n_t),omega_s ) + eps_ts

# Shape into longform data-frame and add error
DF = data.frame( expand.grid(time=1:n_t, x=1:n_x, y=1:n_y),
                   var = "logn",
                   mu = exp(as.vector(p_ts)) )
DF$n = tweedie::rtweedie( n=nrow(DF), mu=DF$mu, phi=1, power=1.5 )

#
Data = DF[ DF$t %in% t_obs, ]
mean( Data$n==0 )

# make mesh
mesh = fm_mesh_2d( Data[,c('x','y')] )

# fit model
mytinyVAST = tinyVAST(
           space_term = "logn <-> logn, sd_space",
           spacetime_term = "logn -> logn, 1, rho
                             logn <-> logn, 0, sd_spacetime",
           data = Data,
           formula = n ~ 1,
           spatial_domain = mesh,
           family = tweedie() )
mytinyVAST

#
DF$nhat = project(
  mytinyVAST,
  extra_times = (max(t_obs)+1):n_t,
  newdata = DF,
  what = "pepsilon_g"
)

library(sf)
data_wide = reshape( DF[,c('x','y','time','mu')],
                     direction = "wide", idvar = c('x','y'), timevar = "time")
sf_data = st_as_sf( data_wide, coords=c("x","y"))
sf_grid = sf::st_make_grid( sf_data )
sf_plot = st_sf( sf_grid, log(st_drop_geometry(sf_data)) )
plot( sf_plot, max.plot=n_t )

dev.new()
data_wide = reshape( DF[,c('x','y','time','nhat')],
                     direction = "wide", idvar = c('x','y'), timevar = "time")
sf_data = st_as_sf( data_wide, coords=c("x","y"))
sf_grid = sf::st_make_grid( sf_data )
sf_plot = st_sf(sf_grid, (st_drop_geometry(sf_data)) )
plot(sf_plot, max.plot=n_t )


##########################
# Bivariate time-series
##########################

library(tinyVAST)

data(isle_royale, package="dsem")

# Convert to long-form
data = expand.grid( "time"=isle_royale[,1], "var"=colnames(isle_royale[,2:3]) )
data$logn = unlist(log(isle_royale[2:3]))

# Define cross-lagged DSEM
dsem = "
  # Link, lag, param_name
  wolves -> wolves, 1, arW
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
  moose -> moose, 1, arM
  #wolves -> moose, 0, corr
  wolves <-> moose, 0, corr
"

# fit model
mytiny = tinyVAST( time_term = dsem,
                 data = data,
                 times = isle_royale[,1],
                 variables = colnames(isle_royale[,2:3]),
                 formula = logn ~ 0 + var )
mytiny

newdata = expand.grid( "time" = max(isle_royale[,1])+1:20, "var" = colnames(isle_royale)[2:3] )
newdata$pred = project(
  mytiny,
  newdata = newdata,
  extra_times = max(isle_royale[,1])+1:20
)

##########################
# Time-series
##########################

library(tinyVAST)

# Convert to long-form
x = cumsum( rnorm(50, mean=0, sd = 0.2) )
data = data.frame( "val" = x, "var" = "x", "time" = seq_along(x) )

# Define cross-lagged DSEM
dsem = "
  x -> x, 1, rho
  x <-> x, 0, sd
"

# fit model
mytiny = tinyVAST(
  time_term = dsem,
  data = data,
  times = unique(data$t),
  variables = "x",
  formula = val ~ 1
)
mytiny

newdata = data.frame( "time" = 1:60, "var" = "x" )
Y = NULL
for(i in 1:100 ){
  tmp = project(
    mytiny,
    newdata = newdata,
    extra_times = 51:60
  )
  Y = cbind(Y, tmp)
}
plot( x = row(Y),
      y = Y )
