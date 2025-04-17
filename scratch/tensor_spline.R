
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force=TRUE, dep = FALSE )
library(tinyVAST)

  set.seed(1)
  dat <- mgcv::gamSim(1, n = 400, scale = 2)
  dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
  dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
  dat$y <- as.numeric(dat$y)

  form = y ~ s(x1) + s(x2) + ti(x1, x2)
  m_m <- mgcv::gam( form, data = dat)
  m_v <- tinyVAST( form, data = dat)

  #
  form = y ~ te(x1, x2)
  m_m <- mgcv::gam( form, data = dat)
  m_v <- tinyVAST( form, data = dat)


# Fix add predictions
if( FALSE ){
  object = m_v
  newdata = object$data
  remove_origdata = FALSE

  add_predictions( m_v, newdata )

}


formula = y ~ t2(x1, x2)
data = dat

library(tinyVAST)
time_term = NULL
space_term = NULL
spacetime_term = NULL
family = gaussian()
space_columns = c("x","y")
spatial_domain = NULL
time_column = "time"
times = NULL
variable_column = "var"
variables = NULL
distribution_column = "dist"
delta_options = list(formula = ~ 1)
spatial_varying = NULL
weights = NULL
control = tinyVASTcontrol()

library(checkmate)
library(Matrix)

##################
# After saving image
##################

setwd( R'(C:\Users\James.Thorson\Desktop\Temporary (can be deleted!))' )
load( "debugging.RData" )

library(TMB)

setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
dyn.unload(dynlib("tinyVAST"))
compile("tinyVAST.cpp" , framework = "TMBad" )
dyn.load(dynlib("tinyVAST"))

# Local fix
#tmb_data$S2dims = numeric(0)
#tmb_data$S2block = numeric(0)

obj = MakeADFun( data = tmb_data,
                 parameters = tmb_par,
                 map = tmb_map,
                 random = tmb_random,
                 DLL = "tinyVAST",
                 profile = control$profile,
                 silent = control$silent )  #

opt = nlminb( obj$par, obj$fn, obj$gr )
sdrep = sdreport(obj)

# Compare
  tmp <- mgcv::gam(formula, data = dat)
