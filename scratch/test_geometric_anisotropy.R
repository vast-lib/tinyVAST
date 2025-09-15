
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
# EXPERIMENT ... is WORKING
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

############################
# Larger domain
############################

#
library(fmesher)
library(tinyVAST)
library(surveyjoin)
library(sf)
library(rnaturalearth)

d <- get_data("sablefish")
dat = as.data.frame(d)[,c('survey_name','lat_start','lon_start','effort','catch_weight')]
dat = na.omit(dat)

mesh = fm_mesh_2d( dat[,c('lon_start','lat_start')], cutoff = 0.5 )

#
nn = RANN::nn2( data = dat[,c('lon_start','lat_start')], query = mesh$loc[,1:2], k = 1 )
survey_s = dat[nn$nn.idx, 'survey_name']

mesh_cov = add_mesh_covariates(
  mesh = mesh,
  data = dat,
  covariates = c(),
  coords = c("lon_start","lat_start")
)
mesh_cov$vertex_covariates = data.frame(
  x = (mesh$loc[,1]),
  y = (mesh$loc[,2])
)
mesh_cov$vertex_covariates = cbind(
  mesh_cov$vertex_covariates,
  x_goa = ifelse( survey_s %in% c("Gulf of Alaska","Aleutian Islands"), mesh$loc[,1], NA ),
  y_goa = ifelse( survey_s %in% c("Gulf of Alaska","Aleutian Islands"), mesh$loc[,2], NA ),
  x_bs = ifelse( survey_s %in% c("Bering Sea Slope","eastern Bering Sea","northern Bering Sea"), mesh$loc[,1], NA ),
  y_bs = ifelse( survey_s %in% c("Bering Sea Slope","eastern Bering Sea","northern Bering Sea"), mesh$loc[,2], NA ),
  x_wc = ifelse( survey_s %in% c("Gulf of Alaska","Aleutian Islands","Bering Sea Slope","eastern Bering Sea","northern Bering Sea"), NA, mesh$loc[,1] ),
  y_wc = ifelse( survey_s %in% c("Gulf of Alaska","Aleutian Islands","Bering Sea Slope","eastern Bering Sea","northern Bering Sea"), NA, mesh$loc[,2] )
)

map_sf = ne_countries( country = c("united states of america","russia","canada") )
#map_sf = st_union(map_sf)
intersected = st_intersects( st_make_valid(st_set_crs(fm_as_sfc(mesh), st_crs(map_sf))), st_make_valid(map_sf) )
mesh_cov$triangle_covariates$barrier_polygon = ifelse( is.na(as.numeric(intersected)), 0, 1 )

fit0 = tinyVAST(
  data = dat,
  formula = catch_weight ~ 1 + offset(log(effort)),
  space_term = "",
  space_columns = c("lon_start","lat_start"),
  spatial_domain = mesh_cov,
  family = tweedie(link = "log"),
  control = tinyVASTcontrol(
    trace=1,
    getsd = FALSE,
    profile = c("alpha_j"),
    use_anisotropy = FALSE
  ),
  development = list(
    #triangle_formula = ~ offset(barrier_polygon)
  )
)

fit1 = tinyVAST(
  data = dat,
  formula = catch_weight ~ 1 + offset(log(effort)),
  space_term = "",
  space_columns = c("lon_start","lat_start"),
  spatial_domain = mesh_cov,
  family = tweedie(link = "log"),
  control = tinyVASTcontrol(
    trace=1,
    getsd = FALSE,
    profile = c("alpha_j"),
    use_anisotropy = TRUE
  ),
  development = list(
    #triangle_formula = ~ offset(barrier_polygon)
  )
)


fit2 = tinyVAST(
  data = dat,
  formula = catch_weight ~ 1 + offset(log(effort)),
  space_term = "",
  space_columns = c("lon_start","lat_start"),
  spatial_domain = mesh_cov,
  family = tweedie(link = "log"),
  control = tinyVASTcontrol(
    trace=1,
    getsd = FALSE,
    profile = c("alpha_j"),
    use_anisotropy = FALSE
  ),
  development = list(
    #triangle_formula = ~ offset(barrier_polygon),
    vertex_formula = ~ 0 + I(x + y) + y
  )
)

fit3 = tinyVAST(
  data = dat,
  formula = catch_weight ~ 1 + offset(log(effort)),
  space_term = "",
  space_columns = c("lon_start","lat_start"),
  spatial_domain = mesh_cov,
  family = tweedie(link = "log"),
  control = tinyVASTcontrol(
    trace=1,
    getsd = FALSE,
    profile = c("alpha_j"),
    use_anisotropy = FALSE
  ),
  development = list(
    #triangle_formula = ~ offset(barrier_polygon),
    vertex_formula = ~
      #I(x + y) + y +
      I(x_goa + y_goa) + y_goa +
      I(x_wc + y_wc) + y_wc +
      #I(x_bs + y_bs) + y_bs +
      0
  )
)

aic = sapply( FUN = AIC, list(fit0,fit1,fit2,fit3) )
aic - min(aic)

mesh_sf = fm_as_sfc( mesh )
domain_sf = st_union( mesh_sf )
grid_sf = st_make_grid( domain_sf, cellsize = 0.5*c(1,1) )
grid_sf = st_intersection(grid_sf, domain_sf)
pred = st_coordinates(st_centroid(grid_sf))

r1 = spatial_cor(
  fit3$rep$Q,
  mesh = mesh_cov,
  coord = c(-122, 36),
  pred = pred
)
r2 = spatial_cor(
  fit3$rep$Q,
  mesh = mesh_cov,
  coord = c(-130, 53),
  pred = pred
)
r3 = spatial_cor(
  fit3$rep$Q,
  mesh = mesh_cov,
  coord = c(-153, 57),
  pred = pred
)
r4 = spatial_cor(
  fit3$rep$Q,
  mesh = mesh_cov,
  coord = c(-171, 61),
  pred = pred
)

#png( "aniso.png", width = 8, height = 8, res=200, units = "in" )
  par(mfrow=c(2,2))
  plot_sf = st_sf( st_geometry(grid_sf), r1 = r1, r2 = r2, r3 = r3, r4 = r4 )
  for(i in 1:4){
    plot(plot_sf[,i], border = NA, reset = FALSE, key.pos = NULL )
    plot(map_sf, add=TRUE, col = "black", border=NA)
  }
#dev.off()
