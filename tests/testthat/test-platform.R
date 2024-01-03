# context("Testing cross platform and R version compatibility")
# platform test
test_that("tinyVAST example is working ", {
  # Simulate a 2D AR1 spatial process with a cyclic confounder w
  n_x = n_y = 25
  n_w = 10
  R_xx = exp(-0.4 * abs(outer(1:n_x, 1:n_x, FUN="-")) )
  R_yy = exp(-0.4 * abs(outer(1:n_y, 1:n_y, FUN="-")) )
  z = mvtnorm::rmvnorm(1, sigma=kronecker(R_xx,R_yy) )

  # Simulate nuisance parameter z from oscillatory (day-night) process
  w = sample(1:n_w, replace=TRUE, size=length(z))
  Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), w=w, z=as.vector(z) + cos(w/n_w*2*pi))
  Data$n = Data$z + rnorm(nrow(Data), sd=1)

  # Add columns for multivariate and temporal dimensions
  Data$var = "n"
  Data$time = 2020

  # make mesh
  mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], n=100 )

  # fit model
  out = fit( data = Data,
             formula = n ~ s(w),
             spatial_graph = mesh,
             control = tinyVASTcontrol(quiet=TRUE, trace=0),
             sem = "" )
  expect_s3_class(out, "tinyVAST")
})

