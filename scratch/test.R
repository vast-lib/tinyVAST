 #
    library(tinyVAST)
     set.seed(123)  # HAS ERRORS

###############
# tinyVAST
###############

# Simulate a seperable two-dimensional AR1 spatial process
n_x = n_y = 25
n_w = 10
R_xx = exp(-0.4 * abs(outer(1:n_x, 1:n_x, FUN="-")) )
R_yy = exp(-0.4 * abs(outer(1:n_y, 1:n_y, FUN="-")) )
z = mvtnorm::rmvnorm(1, sigma=kronecker(R_xx,R_yy) )

# Simulate nuissance parameter z from oscillatory (day-night) process
w = sample(1:n_w, replace=TRUE, size=length(z))
Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), w=w, z=as.vector(z) + cos(w/n_w*2*pi))
Data$n = Data$z + rnorm(nrow(Data), sd=1)

# Add columns for multivariate and/or temporal dimensions
Data$var = "n"

# make SPDE mesh for spatial term
mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], n=100 )

# fit model with cyclic confounder as GAM term
out = tinyVAST( data = Data,
                formula = n ~ s(w),
                spatial_domain = mesh,
                space_term = "n <-> n, sd_n" )

###############
# cAIC
###############

data( red_snapper )
red_snapper = droplevels(subset(red_snapper, Data_type=="Biomass_KG"))

# Define mesh
mesh = fmesher::fm_mesh_2d( red_snapper[,c('Lon','Lat')],
                           cutoff = 1 )

# define formula with a catchability covariate for gear
formula = Response_variable ~ factor(Year) + offset(log(AreaSwept_km2))

# make variable column
red_snapper$var = "logdens"
# fit using tinyVAST
fit = tinyVAST( data = red_snapper,
                formula = formula,
                sem = "logdens <-> logdens, sd_space",
                space_columns = c("Lon",'Lat'),
                spatial_graph = mesh,
                family = tweedie(link="log"),
                variable_column = "var",
                control = tinyVASTcontrol( getsd = FALSE,
                                           profile = "alpha_j" ) )

cAIC(fit) # conditional AIC
AIC(fit) # marginal AIC

###############
# make_dsem_ram
###############


# Univariate AR1
dsem = "
  X -> X, 1, rho
  X <-> X, 0, sigma
"
make_dsem_ram( dsem=dsem, variables="X", times=1:4 )

# Univariate AR2
dsem = "
  X -> X, 1, rho1
  X -> X, 2, rho2
  X <-> X, 0, sigma
"
make_dsem_ram( dsem=dsem, variables="X", times=1:4 )

# Bivariate VAR
dsem = "
  X -> X, 1, XtoX
  X -> Y, 1, XtoY
  Y -> X, 1, YtoX
  Y -> Y, 1, YtoY
  X <-> X, 0, sdX
  Y <-> Y, 0, sdY
"
make_dsem_ram( dsem=dsem, variables=c("X","Y"), times=1:4 )

# Dynamic factor analysis with one factor and two manifest variables
# (specifies a random-walk for the factor, and miniscule residual SD)
dsem = "
  factor -> X, 0, loadings1
  factor -> Y, 0, loadings2
  factor -> factor, 1, NA, 1
  X <-> X, 0, NA, 0           # No additional variance
  Y <-> Y, 0, NA, 0           # No additional variance
"
make_dsem_ram( dsem=dsem, variables=c("X","Y","factor"), times=1:4 )

# ARIMA(1,1,0)
dsem = "
  factor -> factor, 1, rho1 # AR1 component
  X -> X, 1, NA, 1          # Integrated component
  factor -> X, 0, NA, 1
  X <-> X, 0, NA, 0         # No additional variance
"
make_dsem_ram( dsem=dsem, variables=c("X","factor"), times=1:4 )

# ARIMA(0,0,1)
dsem = "
  factor -> X, 0, NA, 1
  factor -> X, 1, rho1     # MA1 component
  X <-> X, 0, NA, 0        # No additional variance
"
make_dsem_ram( dsem=dsem, variables=c("X","factor"), times=1:4 )

###############
# make_eof_ram
###############

make_eof_ram( times = 2010:2020, variables = c("pollock","cod"), n_eof=2 )

###############
# simulate.tinyVAST
###############

set.seed(101)
x = seq(0, 2*pi, length=100)
y = sin(x) + 0.1*rnorm(length(x))
fit = tinyVAST( data=data.frame(x=x,y=y), formula = y ~ s(x) )
simulate(fit, nsim=100, type="mle-mvn")

if(requireNamespace("DHARMa")){
  # simulate new data conditional on fixed effects
  # and sampling random effects from their predictive distribution
  y_iz = simulate(fit, nsim=500, type="mle-mvn")

  # Visualize using DHARMa
  res = DHARMa::createDHARMa( simulatedResponse = y_iz,
                      observedResponse = y,
                      fittedPredictedResponse = fitted(fit) )
  plot(res)
}


###############
# sample_variable
###############

set.seed(101)
 x = runif(n = 100, min = 0, max = 2*pi)
 y = 1 + sin(x) + 0.1 * rnorm(100)

 # Do fit with getJointPrecision=TRUE
 fit = tinyVAST( formula = y ~ s(x),
                 data = data.frame(x=x,y=y),
                 control = tinyVASTcontrol(getJointPrecision = TRUE) )

 # samples from distribution for the mean
 sample_variable(fit)