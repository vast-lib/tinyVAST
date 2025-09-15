
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
           #space_term = "logn <-> logn, sd_space",
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
  what = "mu_g",
  future_var = FALSE
)

library(sf)
data_wide = reshape( DF[,c('x','y','time','mu')],
                     direction = "wide", idvar = c('x','y'), timevar = "time")
sf_data = st_as_sf( data_wide, coords=c("x","y"))
sf_grid = sf::st_make_grid( sf_data )
sf_plot = st_sf( sf_grid, (st_drop_geometry(sf_data)) )
plot( sf_plot, max.plot=n_t, logz = TRUE )

dev.new()
data_wide = reshape( DF[,c('x','y','time','nhat')],
                     direction = "wide", idvar = c('x','y'), timevar = "time")
sf_data = st_as_sf( data_wide, coords=c("x","y"))
sf_grid = sf::st_make_grid( sf_data )
sf_plot = st_sf(sf_grid, log(st_drop_geometry(sf_data)) )
plot(sf_plot, max.plot=n_t, logz=FALSE, breaks = seq(log(min(DF$nhat)),log(max(DF$nhat)),length=11) )


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
newdata$logn = project(
  mytiny,
  newdata = newdata,
  extra_times = max(isle_royale[,1])+1:20,
  future_var = FALSE
)

library(ggplot2)
ggplot( rbind(data,newdata) ) +
  geom_point( aes(x=time, y=logn, col = var) ) #+
  #facet_wrap( vars(var) )



##########################
# Time-series
##########################

library(tinyVAST)
set.seed(123)

# Convert to long-form
n_obs = 100
rho = 0.9
sigma_x = 0.2
sigma_y = 0.1
x = rnorm(n_obs, mean=0, sd = sigma_x)
for(i in 2:length(x)) x[i] = rho * x[i-1] + x[i]
y = x + rnorm( length(x), mean = 0, sd = sigma_y )
data = data.frame( "val" = y, "var" = "y", "time" = seq_along(y) )

# Define cross-lagged DSEM
dsem = "
  y -> y, 1, rho
  y <-> y, 0, sd
"

# fit model
mytiny = tinyVAST(
  time_term = dsem,
  data = data,
  times = unique(data$t),
  variables = "y",
  formula = val ~ 1
)
mytiny

extra_times = length(x) + 1:20
n_sims = 3
newdata = data.frame( "time" = c(seq_along(x),extra_times), "var" = "y" )
Y = NULL
for(i in seq_len(n_sims) ){
  tmp = project(
    mytiny,
    newdata = newdata,
    extra_times = extra_times,
    future_var = FALSE
  )
  Y = cbind(Y, tmp)
}
plot( x = row(Y),
      y = Y )


################################
# Run VAST projection demo
################################

# Load package
setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\scratch)' )
library(VAST)

# load data set
example = load_example( data_set="EBS_pollock" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 50,
  Region = example$Region,
  purpose = "index2",
  bias.correct = FALSE )
settings$RhoConfig[] = 4

# Run model
vast = fit_model( settings = settings,
  Lat_i = example$sampling_data[,'Lat'],
  Lon_i = example$sampling_data[,'Lon'],
  t_i = example$sampling_data[,'Year'],
  b_i = example$sampling_data[,'Catch_KG'],
  a_i = example$sampling_data[,'AreaSwept_km2'],
  newtonsteps = 0 )

# Generate projections
Index_tr = numeric(0)
for( rI in 1:100 ){
  Sim = project_model( x = vast,
                        n_proj = 80,
                        new_covariate_data = NULL,
                        historical_uncertainty = "random",
                        seed = rI )
  Index_tr = cbind( Index_tr, strip_units(Sim$Index_ctl[1,,1]) )
  #Index_tr = cbind( Index_tr, Sim$mean_Z_ctm[1,,2] )
}

# Plot 90% interval for the index
Index_tq = t(apply(Index_tr, MARGIN=1, FUN=quantile, probs=c(0.1,0.5,0.9) ))
Index_tq[,2] = apply(Index_tr, MARGIN=1, FUN=median)
matplot( y=Index_tq, x=rownames(Index_tq), log="y", col="black", lwd=c(1,2,1), type="l", lty="solid" )


##############################
# Compare with tinyVAST
##############################

# make extrapolation list
extrap = data.frame(
  vast$extrapolation_list$Data_Extrap[,c("Lon","Lat")],
  AreaSwept_km2 = strip_units(vast$extrapolation_list$Area_km2_x)
)

#
years = min(dat$Year):(max(dat$Year)+80)
newdata = expand.grid( Year = years, Row = seq_len(nrow(extrap)) )
newdata = cbind( newdata, extrap[newdata$Row,] )

#
library(fmesher)
library(tinyVAST)

# load data set
example = VAST::load_example( data_set="EBS_pollock" )
dat = example$sampling_data
dat$var = "pollock"

#
mesh = fm_mesh_2d( dat[,c("Lon","Lat")], cutoff = 1 )

#
space_term = "
  pollock <-> pollock, sd_space
"
time_term = "
  pollock <-> pollock, 0, sd_time
  pollock -> pollock, 1, ar1_time
"
spacetime_term = "
  pollock <-> pollock, 0, sd_spacetime
  pollock -> pollock, 1, ar1_spacetime
"

#
tv = tinyVAST(
  formula = Catch_KG ~ 1 + offset(log(AreaSwept_km2)),
  data = dat,
  space_term = space_term,
  spacetime_term = spacetime_term,
  time_term = time_term,
  family = delta_gamma( type = "poisson-link" ),
  delta_options = list(
    formula = ~ 1,
    space_term = space_term,
    spacetime_term = spacetime_term,
    time_term = time_term
  ),
  spatial_domain = mesh,
  control = tinyVASTcontrol(
    trace = 1,
    use_anisotropy = TRUE
  ),
  space_columns = c("Lon","Lat"),
  time_column = "Year",
  variable_column = "var",
  variables = "pollock",
  times = min(dat$Year):max(dat$Year)
)

Index2_tr = numeric(0)
for( rI in 1:100 ){
  message( Sys.time(), ": Finished ", rI )
  proj = project(
    object = tv,
    newdata = newdata,
    extra_times = (max(dat$Year)+1):max(newdata$Year),
    past_var = TRUE,
    future_var = TRUE,
    parm_var = FALSE
  )
  index_t = tapply( proj, INDEX = newdata$Year, FUN = sum )
  Index2_tr = cbind( Index2_tr, index_t )
}

# Plot 90% interval for the index
Index2_tq = t(apply(Index2_tr, MARGIN=1, FUN=quantile, probs=c(0.1,0.5,0.9) ))
Index2_tq[,2] = apply(Index2_tr, MARGIN=1, FUN=median)
matplot( y=Index2_tq, x=rownames(Index2_tq), log="y", col="black", lwd=c(1,2,1), type="l", lty="solid" )
