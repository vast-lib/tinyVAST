


test_that("Basic index standardization works", {
  library(fmesher)
  set.seed(101)
  options("tinyVAST.verbose" = FALSE)

  # Simulate settings
  theta_xy = 0.4
  n_x = n_y = 10
  n_t = 5
  rho = 0.8
  spatial_sd = 0.5
  time_sd = 0.5

  # Simulate GMRFs
  R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
  V_ss = spatial_sd^2*kronecker(R_s, R_s)
  d = mvtnorm::rmvnorm(n_t, sigma=V_ss )
  nu_t = rnorm( n_t, mean=0.5, sd = time_sd)

  # Project through time and add mean
  for( t in seq_len(n_t) ){
    if(t>1) d[t,] = rho*d[t-1,] + d[t,]
  }
  d = sweep( d, MARGIN=1, FUN="+", STATS=nu_t)

  # Shape into longform data-frame and add error
  Data = data.frame( expand.grid(time=1:n_t, x=1:n_x, y=1:n_y), "var"="logn", z=exp(as.vector(d)))
  Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$z, phi=0.5, power=1.5 )
  mean(Data$n==0)

  # make mesh
  mesh = fm_mesh_2d( Data[,c('x','y')] )

  # fit with spacetime random-walk using GMRF-projection
  my1 = tinyVAST( spacetime_term = "logn -> logn, 1, NA, 1",
             data = Data,
             formula = n ~ 0 + factor(time),
             spatial_domain = mesh,
             family = delta_gamma(type="poisson-link"),
             control = tinyVASTcontrol(gmrf="proj") )
  # fit model with random walk using standard GMRF
  my2 = tinyVAST( spacetime_term = "logn -> logn, 1, NA, 1",
             data = Data,
             formula = n ~ 0 + factor(time),
             spatial_domain = mesh,
             family = delta_gamma(type="poisson-link"),
             control = tinyVASTcontrol(gmrf="sep") )
  expect_equal( my1$opt, my2$opt, tolerance=0.001 )

  # Predicted sample-weighted total
  index = integrate_output( my1,
                    newdata = subset(Data,time==t) )

  # fit with time & spacetime random-walk using GMRF-projection
  my3 = tinyVAST( spacetime_term = "logn -> logn, 1, NA, 1",
             time_term = "logn -> logn, 1, NA, 0",
             data = Data,
             formula = n ~ 1,
             spatial_domain = mesh,
             family = delta_gamma(type="poisson-link"),
             control = tinyVASTcontrol(gmrf="proj") )
  # fit with time & spacetime random walk using standard GMRF
  my4 = tinyVAST( spacetime_term = "logn -> logn, 1, NA, 1",
             time_term = "logn -> logn, 1, NA, 0",
             data = Data,
             formula = n ~ 1,
             spatial_domain = mesh,
             family = delta_gamma(type="poisson-link"),
             control = tinyVASTcontrol(gmrf="sep") )
  expect_equal( my3$opt, my4$opt, tolerance=0.001 )
})

test_that("Index standardization results are identical in VAST, tinyVAST, and sdmTMB", {
  skip_if_not(require(VAST))
  skip_if_not(require(sdmTMB))
  library(fmesher)
  set.seed(101)
  options("tinyVAST.verbose" = FALSE)

  #
  data(red_snapper)
  data(red_snapper_shapefile)

  #
  Data = subset( red_snapper, Data_type == "Biomass_KG")

  #
  settings = make_settings(
    purpose = "index3",
    Region = "Other",
    n_x = 100,
    use_anisotropy = FALSE
  )
  # Eliminate Omega/Epsilon for 2nd linear predictor because tinyVAST has only one kappa
  settings$FieldConfig[c("Omega","Epsilon"),'Component_2'] = 0
  vast = fit_model(
    settings = settings,
    Lat_i = Data[,'Lat'],
    Lon_i = Data[,'Lon'],
    t_i = Data[,'Year'],
    b_i = Data[,'Response_variable'],
    a_i = rep(1,nrow(Data)),
    grid_dim_km = c(0.1,0.1),
    projargs = "EPSG:4326",
    observations_LL = cbind('Lat'=Data[,'Lat'], 'Lon'=Data[,'Lon'])
  )

  # fit model with random walk using standard GMRF
  # Won't exactly match because using a single log_kappa for both linear predictors'
  tv = tinyVAST(
    spacetime_term = "",
    space_term = "",
    data = droplevels(Data),
    space_columns = c("Lon","Lat"),
    variable_column = "Data_type",
    time_column = "Year",
    formula = Response_variable ~ 0 + factor(Year),
    spatial_domain = vast$spatial_list$Mesh$anisotropic_mesh,
    family = delta_gamma(type="poisson-link"),
    delta_options = list(
      formula = ~ 0 + factor(Year)
    ),
    control = tinyVASTcontrol( trace = 1 )
  )

  # Predicted index using VAST grid
  extrap = vast$extrapolation_list
  index = sapply(
    sort(unique(Data$Year)),
    FUN = \(t){
      integrate_output(
        tv,
        newdata = data.frame( extrap$Data_Extrap,
                              Year = t,
                              Data_type = "Biomass_KG" ),
        area = extrap$Area_km2,
      )
    }
  )
  # cbind( tv$rep$mu_i, vast$Report$D_i )

  # tinyVAST index
  index_vast = rbind(
    "Estimate" = as.list(vast$par$SD, report=TRUE, what="Estimate")$Index_ctl[1,,1],
    "Std. Error" = as.list(vast$par$SD, report=TRUE, what="Std. Error")$Index_ctl[1,,1],
    "Est. (bias.correct)" = as.list(vast$par$SD, report=TRUE, what="Est. (bias.correct)")$Index_ctl[1,,1]
  )
  index_tv = as.matrix(index[1:3,])

  mesh2 = make_mesh( data = droplevels(Data),
                     xy_cols = c("Lon","Lat"),
                     mesh = vast$spatial_list$Mesh$anisotropic_mesh )
  sdm = sdmTMB(
    formula = Response_variable ~ 0 + factor(Year),
    data = droplevels(Data),
    mesh = mesh2,
    family = delta_gamma(type="poisson-link"),
    spatiotemporal = list("iid","off"),
    time = "Year",
    spatial = list("on","off"),
    share_range = TRUE,
    do_index = TRUE
  )
  pred_sdm = predict( sdm,
                      newdata = data.frame(extrap$Data_Extrap),
                      xy_cols = c("Lon","Lat"),
                      return_tmb_object = TRUE,
                      year =  )
  index_sdm = sapply(
    sort(unique(Data$Year)),
    FUN = \(t){
      pred_sdm = predict( sdm,
                          newdata = data.frame(extrap$Data_Extrap, Year = t),
                          xy_cols = c("Lon","Lat"),
                          return_tmb_object = TRUE )
      index_sdm = get_index( pred_sdm,
                             bias_correct = TRUE,
                             area = as.numeric(extrap$Area_km2) )[1:6]
      index_sdm$est
    }
  )

  # index for tinyVAST vs. VAST
  expect_equal( index_tv, index_vast, tolerance = 1 )
  # index for tinyVAST vs. sdmTMB
  expect_equal( index_tv['Est. (bias.correct)',], index_sdm, tolerance = 1 )
  # ML for tinyVAST vs. VAST
  expect_equal( tv$opt$obj, as.numeric(vast$parameter_estimates$obj), tolerance = 0.1 )
  # ML for tinyVAST vs. sdmTMB
  expect_equal( tv$opt$obj, as.numeric(sdm$model$obj), tolerance = 0.1 )
})



