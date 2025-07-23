
library(sdmTMBextra)
?add_barrier_mesh

library(sdmTMB)
library(dplyr)
library(ggplot2)

# For applied situations on finer scales, you may with to use scale = "large".
# For that, first: remotes::install_github("ropensci/rnaturalearthhires")
map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "canada")

# Crop the polygon for plotting and efficiency:
sf::st_bbox(map_data)
bc_coast <-
  sf::st_crop(map_data, c(xmin = -134, ymin = 46, xmax = -120, ymax = 57))

crs_utm9 <- 3156 # Pick a projection, here UTM9

sf::st_crs(bc_coast) <- 4326 # 'WGS84'
bc_coast <- sf::st_transform(bc_coast, crs_utm9)

# Project our survey data coordinates:
dat <- pcod
survey <- dat %>% dplyr::select(lon, lat, density) %>%
  sf::st_as_sf(crs = 4326, coords = c("lon", "lat")) %>%
  sf::st_transform(crs_utm9)

# Plot our coast and survey data:
ggplot(bc_coast) +
  geom_sf() +
  geom_sf(data = survey, size = 0.5)

# Note that a barrier mesh won't do much here for this
# example data set, but we nonetheless use it as an example.

# Prepare for making the mesh
# First, we will extract the coordinates:
surv_utm_coords <- sf::st_coordinates(survey)

# Then we will scale coordinates to km so the range parameter
# is on a reasonable scale for estimation:
dat$X1000 <- surv_utm_coords[,1] / 1000
dat$Y1000 <- surv_utm_coords[,2] / 1000

mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X1000", "Y1000"),
  n_knots = 200, type = "kmeans")
plot(mesh)

# Add on the barrier mesh component:
bspde <- sdmTMBextra::add_barrier_mesh(
  mesh, bc_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

# In the above, the grey dots are the centre of triangles that are in the
# ocean. The red crosses are centres of triangles that are over land. The
# spatial range will be assumed to be 0.1 (`range_fraction`) over land
# compared to over water.

# We can make a more advanced plot if we want:
mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
mesh_df_land <- bspde$mesh_sf[bspde$barrier_triangles, ]

ggplot(bc_coast) +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green")

# Now, when we fit our model with the new mesh, it will automatically
# include a barrier structure in the spatial correlation:
fit <- sdmTMB(density ~ 1, # s(depth, k = 3),
              data = pcod, mesh = bspde,
  family = tweedie(link = "log"))
fit

#
fit0 <- sdmTMB(density ~ 1, # s(depth, k = 3),
               data = pcod, mesh = mesh,
  family = tweedie(link = "log"))
fit0

sdmTMB_mesh = mesh

######################
# Replicate in tinyVAST
######################


library(tinyVAST)
library(sf)

mesh_with_covs = mesh = sdmTMB_mesh$mesh

# Manually add covariates for now
boundary_sf = st_geometry( bc_coast )
vertex_sf = st_as_sf( as.data.frame(mesh$loc[,1:2]) * 1000, crs = st_crs(boundary_sf), coords = c("V1","V2") )
intersected = st_intersects( vertex_sf, boundary_sf )
which_land = which(lengths(intersected) == 1)
mesh_with_covs$vertex_covariates = data.frame( "land" = ifelse(lengths(intersected)==1, 10, 0) )
class(mesh_with_covs) <- c("vertex_coords", class(mesh_with_covs))

# Run
tv0 <- tinyVAST(
  density ~ 1,  # s(depth, k = 3)
  data = as.data.frame(pcod),
  spatial_domain = mesh,
  space_term = "",
  space_columns = c("X","Y"),
  family = tweedie(link = "log")
)

cbind(
  "tinyVAST" = c(tv0$opt$par, AIC = AIC(tv0), cAIC = tinyVAST::cAIC(tv0)),
  "sdmTMB" = c(fit0$model$par, AIC(fit0), sdmTMB::cAIC(fit0))
)

# Run
tv <- tinyVAST(
  density ~ 1,  # s(depth, k = 3)
  data = as.data.frame(pcod),
  spatial_domain = mesh_with_covs,
  space_term = "",
  development = list(
    kappa_formula = ~ land
  ),
  space_columns = c("X","Y"),
  family = tweedie(link = "log")
)

cbind(
  "tinyVAST" = c(tv$opt$par, AIC = AIC(tv), cAIC = tinyVAST::cAIC(tv)),
  "sdmTMB" = c(fit$model$par, NA, AIC(fit), sdmTMB::cAIC(fit))
)


