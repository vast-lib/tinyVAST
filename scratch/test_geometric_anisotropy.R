
# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

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

#
(H = tv$rep$H)

#
#library(RConics)
#adjoint(H)
adjoint = function(A){
  minor = function(A,i,j) det(A[-i, -j, drop = FALSE])
  cofactor = function(A,i,j) (-1)^(i + j) * minor(A, i, j)
  n <- nrow(A)
  t(outer(1:n, 1:n, Vectorize(function(i, j) cofactor(A, i, j))))
}
adjoint(H)

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
    #vertex_formula = ~ 0 + I(x+y) + x
    vertex_formula = ~ 0 + I(x+y) + y
  ),
  space_columns = c("Lon","Lat"),
  time_column = "Year",
  variable_column = "var",
  variables = "pollock",
  times = min(dat$Year):max(dat$Year)
)
tv2$sdrep
c(AIC(tv), AIC(tv2))

adjoint(tv2$rep$H)

#####################
# EXPERIMENT #2 ... not working
#####################


# Extract
loc = mesh$loc[,1:2]
TV = mesh$graph$tv
V0 <- loc[ TV[,1], ]
V1 <- loc[ TV[,2], ]
V2 <- loc[ TV[,3], ]
E0 <- V2 - V1
E1 <- V0 - V2
E2 <- V1 - V0

W_rz = array(0, dim=c(nrow(TV),3), dimnames=list(NULL,c("sum_dx2","sum_dy2","sum_dxdy")) )
for( r in 1:nrow(TV) ){
  tmp = rbind(E0[r,], E1[r,], E2[r,])
  tmp2 = t(tmp) %*% tmp
  W_rz[r,] = c( tmp2[1,1], tmp2[2,2], tmp2[1,2] )
}

# Add to mesh
mesh_cov = add_mesh_covariates(
  mesh = mesh,
  data = dat,
  covariates = c(),
  coords = c("Lon","Lat")
)
mesh_cov$triangle_covariates = cbind(
  mesh_cov$triangle_covariates,
  W_rz / (W_rz[,1] + W_rz[,2])
)

#
tv3 = tinyVAST(
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
    #triangle_formula = ~ 0 + sum_dx2 + sum_dxdy
    #triangle_formula = ~ 0 + sum_dy2 + sum_dxdy
    triangle_formula = ~ 0 + sum_dx2 + sum_dy2 + sum_dxdy
  ),
  space_columns = c("Lon","Lat"),
  time_column = "Year",
  variable_column = "var",
  variables = "pollock",
  times = min(dat$Year):max(dat$Year)
)


