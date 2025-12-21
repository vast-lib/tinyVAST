
get_grid <-
function( model,
          exclude_terms = NULL,
          length_out = 50,
          values = NULL,
          ... ){

  n_terms <- length(model[["var.summary"]])
  term_list <- list()
  for (term in 1:n_terms) {
    term_summary <- model[["var.summary"]][[term]]
    term_name <- names(model[["var.summary"]])[term]
    if (term_name %in% names(values)) {
      new_term <- values[[which(names(values) == term_name)]]
      if (is.null(new_term)) {
        new_term <- model[["var.summary"]][[term]][[1]]
      }
    }
    else {
      if (is.numeric(term_summary)) {
        min_value <- min(term_summary)
        max_value <- max(term_summary)
        new_term <- seq(min_value, max_value, length.out = length_out)
      }
      else if (is.factor(term_summary)) {
        new_term <- levels(term_summary)
      }
      else {
        stop("The terms are not numeric or factor.\n")
      }
    }
    term_list <- append(term_list, list(new_term))
    names(term_list)[term] <- term_name
  }
  new_data <- expand.grid(term_list)
  return(new_data)
}


library(mvtweedie)
library(mgcv)
library(tinyVAST)
library(fmesher)

#
library(rnaturalearth)
library(raster)
library(sf)
library(dplyr)
library(ggplot2)

# Load data
data(southeast_alaska_wolf)
groups = c("Black tailed deer","Marine mammal", "Mountain goat", "Beaver")
southeast_alaska_wolf = subset( southeast_alaska_wolf,
                                group %in% groups )

#
southeast_alaska_wolf$group = factor(southeast_alaska_wolf$group)

# Illustrate format
knitr::kable( head(southeast_alaska_wolf), digits=1, row.names=FALSE)

###############
# Using mgcv
###############

Formula = Response ~
  0 + group +
  s(Latitude,Longitude,m=c(1,0.5),bs="ds") +
  s(Latitude,Longitude,by=group,m=c(1,0.5),bs="ds")

mygam = gam(
  formula = Formula,
  data = southeast_alaska_wolf,
  family = tw
)
class(mygam) = c( "mvtweedie", class(mygam) )

#############
# Plot for mgcv

# Predict raster on map
new_data = get_grid( mygam,
                 length_out = 100 )
pred <- predict( mygam,
                 newdata = new_data,
                 se.fit = TRUE)
predicted <- as.data.frame(pred)
pred_wolf <- cbind(new_data, predicted)

pred_wolf$cv.fit = pred_wolf$se.fit / pred_wolf$fit

# Map oceanmap layer
US_high = ne_countries(scale=50, country="united states of america")
st_box = st_polygon( list(cbind( x=c(-140,-125,-124,-140,-140),
                                 y=c(50,50,60,60,50))) )
st_box = st_sfc(st_box, crs=st_crs(US_high) )
wmap = st_intersection( US_high, st_box )
oceanmap = st_difference( st_as_sfc(st_bbox(wmap)), wmap )

sf.isl <- data.frame(island = c("Baranof", "Chichagof", "Admiralty"),
                     lat = c(57.20583, 57.88784, 57.59644),
                     lon = c(-135.1866, -136.0024, -134.5776)) %>%
  st_as_sf(., coords = c("lon", "lat"), crs = 4326)

mask.land = ne_countries(scale=50, country="united states of america", returnclass = 'sf') %>%
  st_set_crs(., 4326) %>%
  st_cast(., "POLYGON") %>%
  st_join(., sf.isl) %>%
  filter(!is.na(island))

# Make figure
my_breaks = c(0.02,0.1,0.25,0.5,0.75)
ggplot(oceanmap) +
  geom_tile(data=pred_wolf, aes(x=Longitude,y=Latitude,fill=fit)) +
  geom_sf() +
  geom_sf(data = mask.land, inherit.aes = FALSE, fill = "darkgrey") +
  coord_sf(xlim=range(pred_wolf$Longitude), ylim=range(pred_wolf$Latitude), expand = FALSE) +
  facet_wrap(vars(group), ncol = 5) +
  scale_fill_gradient2(name = "Proportion", trans = "sqrt", breaks = my_breaks) +
  scale_y_continuous(breaks = seq(55,59,2)) +
  scale_x_continuous(breaks = c(-135, -133.5, -132)) +
  theme(axis.text = element_text(size = 7))

##############
# using tinyVAST
##############

mesh = fm_mesh_2d(
  loc = southeast_alaska_wolf[,c("Longitude","Latitude")],
  cutoff = 0.05,
  refine = TRUE
)

southeast_alaska_wolf$category = sapply(
  X = as.character(southeast_alaska_wolf$group),
  FUN = switch,
  "Beaver" = "beaver",
  "Black tailed deer" = "deer",
  "Marine mammal" = "marine",
  "Mountain goat" = "goat"
)

space_term = "
  beaver <-> beaver, sd_beaver, 1
  deer <-> deer, sd_deer, 1
  marine <-> marine, sd_marine, 1
  goat <-> goat, sd_goat, 1
"

mytv = tinyVAST(
  formula = Response ~ 0 + factor(category),
  data = southeast_alaska_wolf,
  family = tweedie(),
  space_term = space_term,
  space_columns = c("Longitude","Latitude"),
  spatial_varying = ~ 1,
  spatial_domain = mesh,
  variable_column = "category"
)

#############
# Plot for mgcv

new_data = get_grid( mygam,
                 length_out = 100 )
# Predict raster on map
new_data$category = sapply(
  X = as.character(new_data$group),
  FUN = switch,
  "Beaver" = "beaver",
  "Black tailed deer" = "deer",
  "Marine mammal" = "marine",
  "Mountain goat" = "goat"
)

#class(mytv) = "tinyVAST"
class(mytv) = c( "mvtweedie", class(mytv) )

pred <- predict.mvtweedie(
  mytv,
  newdata = new_data,
  category_name = "category",
  se.fit = FALSE
)
#predicted <- as.data.frame(pred)
pred_wolf <- cbind(new_data, fit = pred)
#pred_wolf$cv.fit = pred_wolf$se.fit / pred_wolf$fit

# Make figure
my_breaks = c(0.02,0.1,0.25,0.5,0.75)
ggplot(oceanmap) +
  geom_tile(data=pred_wolf, aes(x=Longitude,y=Latitude,fill=fit)) +
  geom_sf() +
  geom_sf(data = mask.land, inherit.aes = FALSE, fill = "darkgrey") +
  coord_sf(xlim=range(pred_wolf$Longitude), ylim=range(pred_wolf$Latitude), expand = FALSE) +
  facet_wrap(vars(group), ncol = 5) +
  scale_fill_gradient2(name = "Proportion", trans = "sqrt", breaks = my_breaks) +
  scale_y_continuous(breaks = seq(55,59,2)) +
  scale_x_continuous(breaks = c(-135, -133.5, -132)) +
  theme(axis.text = element_text(size = 7))
