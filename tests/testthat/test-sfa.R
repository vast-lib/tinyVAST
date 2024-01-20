

test_that("Basic spatial factor analysis works", {
  library(fmesher)
  set.seed(101)

  # Simulate settings
  theta_xy = 0.4
  n_x = n_y = 10
  n_c = 5
  rho = 0.8
  resid_sd = 0.5

  # Simulate GMRFs
  R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
  R_ss = kronecker(X=R_s, Y=R_s)
  delta_fs = mvtnorm::rmvnorm(n_c, sigma=R_ss )

  #
  L_cf = matrix( rnorm(n_c^2), nrow=n_c )
  L_cf[,3:5] = 0
  L_cf = L_cf + resid_sd * diag(n_c)

  #
  d_cs = L_cf %*% delta_fs

  # Shape into longform data-frame and add error
  Data = data.frame( expand.grid(species=1:n_c, x=1:n_x, y=1:n_y),
                     "var"="logn", "z"=exp(as.vector(d_cs)) )
  Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$z, phi=0.5, power=1.5 )
  mean(Data$n==0)

  # make mesh
  mesh = fm_mesh_2d( Data[,c('x','y')] )

  #
  sem = "
    f1 -> 1, l1
    f1 -> 2, l2
    f1 -> 3, l3
    f1 -> 4, l4
    f1 -> 5, l5
    f2 -> 2, l6
    f2 -> 3, l7
    f2 -> 4, l8
    f2 -> 5, l9
    f1 <-> f1, NA, 1
    f2 <-> f2, NA, 1
    1 <-> 1, NA, 0
    2 <-> 2, NA, 0
    3 <-> 3, NA, 0
    4 <-> 4, NA, 0
    5 <-> 5, NA, 0
  "

  # fit model
  out = tinyVAST( sem = sem,
             data = Data,
             formula = n ~ 0 + factor(species),
             spatial_graph = mesh,
             family = tweedie(),
             variables = c( "f1", "f2", 1:n_c ),
             space_columns = c("x","y"),
             variable_column = "species",
             time_column = "time",
             distribution_column = "dist",
             control = tinyVASTcontrol(gmrf="proj") )
  #
  summary( out, what="sem" )
})
