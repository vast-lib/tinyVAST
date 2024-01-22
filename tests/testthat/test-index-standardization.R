


test_that("Basic spatial factor analysis works", {
  library(fmesher)
  set.seed(101)
  options("tinyVAST.verbose" = FALSE)

  # Simulate settings
  theta_xy = 0.4
  n_x = n_y = 10
  n_t = 5
  rho = 0.8
  spatial_sd = 0.5

  # Simulate GMRFs
  R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
  V_ss = spatial_sd^2*kronecker(R_s, R_s)
  d = mvtnorm::rmvnorm(n_t, sigma=V_ss )

  # Project through time and add mean
  for( t in seq_len(n_t) ){
    if(t>1) d[t,] = rho*d[t-1,] + d[t,]
  }
  #d = d + 0.5

  # Shape into longform data-frame and add error
  Data = data.frame( expand.grid(time=1:n_t, x=1:n_x, y=1:n_y), "var"="logn", z=exp(as.vector(d)))
  Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$z, phi=0.5, power=1.5 )
  mean(Data$n==0)

  # make mesh
  mesh = fm_mesh_2d( Data[,c('x','y')] )

  # fit model with random walk using GMRF-projection
  my1 = tinyVAST( dsem = "logn -> logn, 1, NA, 1",
             data = Data,
             formula = n ~ 0 + factor(time),
             spatial_graph = mesh,
             family = delta_poisson_link_gamma(),
             control = tinyVASTcontrol(gmrf="proj") )
  # fit model with random walk using standard GMRF
  my2 = tinyVAST( dsem = "logn -> logn, 1, NA, 1",
             data = Data,
             formula = n ~ 0 + factor(time),
             spatial_graph = mesh,
             family = delta_poisson_link_gamma(),
             control = tinyVASTcontrol(gmrf="sep") )
  expect_equal( my1$opt, my2$opt, tolerance=0.001 )

  # Predicted sample-weighted total
  integrate_output( my1,
                    newdata = subset(Data,time==t) )
})
