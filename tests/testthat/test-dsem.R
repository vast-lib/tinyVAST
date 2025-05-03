

test_that("dsem and tinyVAST give identical results without lags", {

  library(dsem)
  library(igraph)

  # Simulate data
  set.seed(101)
  n_rows = 100
  X = 1 + rnorm(n_rows)
  Y = 2 + 0.5 * X + rnorm(n_rows)

  # Format wide-form for DSEM and long-form fot tinyVAST
  tsdat = ts( data.frame(X=X, Y=Y) )
  dat = expand.grid( times = time(tsdat), var = c("X","Y") )
  dat$site = "a"
  dat$z = as.vector(tsdat)

  #
  f1 = tinyVAST(
    data = dat,
    formula = z ~ 0 + factor(var),
    variable_column = "var",
    time_term = "X -> Y, 0, beta",
    time_column = "times"
  )

  # dsem build inputs
  f2 = dsem(
    tsdata = tsdat,
    family = c( "normal", "normal" ),
    sem = "X -> Y, 0, beta",
    control = dsem_control( run_model = FALSE, use_REML = FALSE )
  )

  # Refit with same measurement error as tinyVAST
  map = f2$tmb_inputs$map
    map$lnsigma_j = factor( c(NA,NA) )
  pars = f2$tmb_inputs$parameters
    pars$lnsigma_j = rep( f1$internal$parlist$log_sigma, 2 )
  f2 = dsem(
    tsdata = ts( data.frame(X=X, Y=Y) ),
    family = c( "normal", "normal" ),
    sem = "X -> Y, 0, beta",
    control = dsem_control( parameters = pars, map = map, use_REML = FALSE )
  )

  #
  unconnected_graph = make_empty_graph( 1 )
  V(unconnected_graph)$name = "a"
  f3 = tinyVAST(
    data = dat,
    formula = z ~ 0 + factor(var),
    variable_column = "var",
    spacetime_term = "X -> Y, 0, beta",
    time_column = "times",
    space_column = "site",
    spatial_domain = unconnected_graph
  )

  # tinyVAST vs. DSEM
  expect_equal( as.numeric(f1$opt$par[1:5]),
                as.numeric(f2$opt$par[c(4:5,1:3)]),
                tolerance = 0.01 )
  # Different specification for tinyVAST
  expect_equal( as.numeric(f1$opt$par[1:5]),
                as.numeric(f3$opt$par[1:5]),
                tolerance = 0.01 )
})

test_that("dsem and tinyVAST give identical results with lags", {
  # Simulate data
  set.seed(101)
  n_rows = 100
  X = 1 + rnorm(n_rows - 1)
  Y = 2 + 0.5 * X + rnorm(n_rows - 1)

  # Add NAs for manually lagged effect
  X = c(X, NA )
  Y = c(NA, Y )

  # Format wide-form for DSEM and long-form fot tinyVAST
  tsdat = ts( data.frame(X=X, Y=Y) )
  dat = expand.grid( times = time(tsdat), var = c("X","Y") )
  dat$site = "a"
  dat$z = as.vector(tsdat)
  dat = na.omit(dat)

  #
  f1 = tinyVAST(
    data = dat,
    formula = z ~ 0 + factor(var),
    variable_column = "var",
    time_term = "X -> Y, 1, beta",
    time_column = "times",
    control = tinyVASTcontrol( extra_reporting = TRUE )
  )

  # dsem build inputs
  f2 = dsem(
    tsdata = tsdat,
    family = c( "normal", "normal" ),
    sem = "X -> Y, 1, beta",
    control = dsem_control( run_model = TRUE, use_REML = FALSE, nlminb_loops = 0, newton_loops = 0,
                            getsd = FALSE, extra_convergence_checks = FALSE )
  )

  # Refit with same measurement error as tinyVAST
  map = f2$tmb_inputs$map
    map$lnsigma_j = factor( c(NA,NA) )
  pars = f2$tmb_inputs$parameters
    pars$lnsigma_j = rep( f1$internal$parlist$log_sigma, 2 )
  f2 = dsem(
    tsdata = tsdat,
    family = c( "normal", "normal" ),
    sem = "X -> Y, 1, beta",
    control = dsem_control( parameters = pars, map = map, use_REML = FALSE, extra_convergence_checks = FALSE )
  )

  expect_equal( as.numeric(f1$opt$par[1:5]),
                as.numeric(f2$opt$par[c(4:5,1:3)]),
                tolerance = 0.01 )
})
