
devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

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

#
tv = tinyVAST(
  formula = Catch_KG ~ 1 + offset(log(AreaSwept_km2)),
  data = dat,
  space_term = space_term,
  family = tweedie(),
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

#####################
# EXPERIMENT
#####################

mesh_cov = add_mesh_covariates(
  mesh = mesh,
  data = dat,
  covariates = c(),
  coords = c("Lon","Lat")
)
mesh_cov$vertex_covariates = data.frame(
  x = scale(mesh$loc[,1]),
  y = scale(mesh$loc[,2]),
  xy = scale(mesh$loc[,1] * mesh$loc[,2])
)

#
tv2 = tinyVAST(
  formula = Catch_KG ~ 1 + offset(log(AreaSwept_km2)),
  data = dat,
  space_term = space_term,
  family = tweedie(),
  spatial_domain = mesh_cov,
  control = tinyVASTcontrol(
    trace = 1,
    use_anisotropy = FALSE
  ),
  development = list(
    vertex_formula = ~ 0 + x + I(x-y)
  ),
  space_columns = c("Lon","Lat"),
  time_column = "Year",
  variable_column = "var",
  variables = "pollock",
  times = min(dat$Year):max(dat$Year)
)
tv2$sdrep
c(AIC(tv), AIC(tv2))
