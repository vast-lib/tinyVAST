

test_that("SVC works", {
  library(fmesher)
  set.seed(101)
  options("tinyVAST.verbose" = FALSE)

  # Simulate settings
  theta_xy = 0.4
  n_x = n_y = 10
  spatial_sd = 0.5

  # Simulate GMRFs
  R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
  V_ss = spatial_sd^2 * kronecker(R_s, R_s)
  xi_s = mvtnorm::rmvnorm(1, sigma=V_ss )[1,]
  logd_s = 1 + xi_s

  # Shape into longform data-frame and add error
  Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), 
                     "var"="logn", 
                     d_s = exp(as.vector(logd_s)) )
  Data$n = tweedie::rtweedie( n = nrow(Data), 
                              mu = Data$d_s, 
                              phi = 0.5, 
                              power = 1.5 )
  mean(Data$n==0)

  # make mesh
  mesh = fm_mesh_2d( Data[,c('x','y')] )

  # fit model with random walk using GMRF-projection
  my1 = tinyVAST( sem = "logn <-> logn, sd",
             data = Data,
             formula = n ~ 1,
             spatial_graph = mesh,
             family = tweedie(),
             control = tinyVASTcontrol() )
  # fit model with random walk using standard GMRF
  my2 = tinyVAST( spatial_varying = ~ 1,
             data = Data,
             formula = n ~ 1,
             spatial_graph = mesh,
             family = tweedie(),
             control = tinyVASTcontrol() )
  expect_equal( my1$opt$objective, my2$opt$objective, tolerance=0.001 )

})
