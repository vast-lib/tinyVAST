
#' Arctic September sea ice concentrations
#'
#' Data used to demonstrate and test empirical orthogonal function
#' generalized linear latent variable model (EOF-GLLVM)
#'
#' @name sea_ice
#' @docType data
#' @usage data(sea_ice)
#' @keywords data
NULL

#' North Pacific salmon returns
#'
#' Data used to demonstrate and test multivariate second-order autoregressive
#' models using a simultaneous autoregressive (SAR) process across regions.
#' Data are from \doi{10.1002/mcf2.10023}
#'
#' @name salmon_returns
#' @docType data
#' @usage data(salmon_returns)
#' @keywords data
NULL

#' Condition and density example
#'
#' Data used to demonstrate and test a bivariate model for morphometric condition
#' (i.e., residuals in a weight-at-length relationship) and density for fishes, using
#' the same example as was provided as a wiki example for VAST.
#' Data are from \doi{10.3354/meps13213}
#'
#' @name condition_and_density
#' @docType data
#' @usage data(condition_and_density)
#' @keywords data
NULL

#' Survey catch-rates at age for Alaska pollock in the Eastern and Northern Bering Sea
#'
#' Data used to demonstrate and test model-based age expansion, using density=
#' dependence corrected survey catch rates after first=stage expansion
#' from the bottom trawl survey for ages 1-15, conducted by
#' by the Alaska Fisheries Science Center, including annual surveys in the eastern 
#' Bering Sea 1982-2019 and 2021-2023, as well as the northern Bering Sea 
#' in 1982/85/88/91 and 2010/17/18/19/21/22/23.
#'
#' @name bering_sea_pollock_ages
#' @docType data
#' @usage data(bering_sea_pollock_ages)
#' @keywords data
NULL

#' Survey domain for the eastern and northern Bering Sea surveys
#'
#' Shapefile defining the spatial domain for the eastern and northern
#' Bering Sea bottom trawl surveys.
#'
#' @name bering_sea
#' @docType data
#' @usage data(bering_sea)
#' @keywords data
NULL

#' Estimated proportion-at-age for Alaska pollock using VAST 
#'
#' Estimated proporrtion-at-age for Alaska pollock using the 
#' package VAST, for comparison with output using tinyVAST.
#'
#' @name bering_sea_pollock_vast
#' @docType data
#' @usage data(bering_sea_pollock_vast)
#' @keywords data
NULL

#' Presence/absence, count, and biomass data for red snapper
#'
#' Data used to demonstrate and test analysis using multiple data types
#'
#' @name red_snapper
#' @docType data
#' @usage data(red_snapper)
#' @keywords data
NULL

#' Shapefile for red snapper analysis
#'
#' Spatial extent used for red snapper analysis, derived from Chap-7 of \doi{10.1201/9781003410294}
#'
#' @name red_snapper_shapefile
#' @docType data
#' @usage data(red_snapper_shapefile)
#' @keywords data
NULL

#' Data to demonstrate model-based diet proportions
#'
#' Data to estimate predator-expanded stomach contents (PESC), i.e.,
#' a multi-prey single-predator model for stomach contents of red grouper
#' in the Gulf of Mexico that also estimates biomass density for red grouper.
#' Diet proportions are then the product of predicted diet proportions and prey
#' specific consumption, normalized across prey categories. Copied from
#' VAST data set \code{PESC_example_red_grouper}
#'
#' @name red_grouper_diet
#' @docType data
#' @usage data(red_grouper_diet)
#' @keywords data
NULL

#' Data to demonstrate probabilistic forecasting
#'
#' Data to estimate probabilistic forecasting, which bridges uncertainty
#' between short-term (decadal) projections and long-term (end-of-century)
#' forecasts.  Here we use capelin sampling from a survey trawl survey
#' in the eastern and northern Bering Sea, as well as annual
#' sea surface temperature anomalies.  The data set includes data frames
#' for both the fitted data \code{Data} and the projected values \code{New_Data}
#'
#' @name bering_sea_capelin_forecasts
#' @docType data
#' @usage data(bering_sea_capelin_forecasts)
#' @keywords data
NULL

#' Data to analyze sponge-coral-fish associations
#'
#' Data used in
#' \href{https://onlinelibrary.wiley.com/doi/10.1111/geb.70035}{Thorson et al. 2025}
#' and copied from \href{https://zenodo.org/records/15001896}{Zenodo} to allow
#' inclusion as a package vignette.  This includes bottom-trawl samples of rockfish
#' and flatfish from the Alaska Fisheries Science Center in the Gulf of Alaska
#' and Aleutian Islands, and drop-camera samples of rockfish, flatfish, sponges,
#' and corals from the Alaska Fisheries Science Center and
#' funded by the Deep Sea Coral Initiative.
#'
#' @name alaska_sponge_coral_fish
#' @docType data
#' @usage data(alaska_sponge_coral_fish)
#' @keywords data
NULL

