
test_that("Basic sfnetworks works", {
  library(sf)
  library(sfnetworks)

  stream <- st_read( file.path(system.file("stream_network",package="tinyVAST"),
                     "East_Fork_Lewis_basin.shp"), quiet=TRUE )
  stream = as_sfnetwork(stream)

  # Rescale
  graph = sfnetwork_mesh( stream )
  graph$table$dist = graph$table$dist / 1000  # Convert distance scale

  # Parameters
  alpha = 2
  kappa = 0.05

  # simulate
  omega_s = simulate_sfnetwork( n=1, sfnetwork_mesh=graph, theta=kappa)[,1]

  # sample locations along network
  extrap = st_union( st_line_sample( activate(stream,"edges"), density=1/10000))
  extrap = st_cast( extrap, "POINT" )

  # Project to sampled locations
  A_is = sfnetwork_evaluator( stream = graph$stream,
                                  loc = st_coordinates(extrap) )
  omega_i = (A_is %*% omega_s)[,1]

  # Simulate sampling
  #Count = rpois( n=graph$n, lambda=exp(alpha + omega) )
  Count_i = rnorm( n=length(omega_i), mean=alpha + omega_i, sd=0.5 )

  # Format into long-form data frame expected by tinyVAST
  Data = data.frame( Count = Count_i,
                     st_coordinates(extrap),
                     var = "species",  # Univariate model so only one value
                     time = "2020",    # no time-dynamics, so only one value
                     dist = "obs" )    # only one type of sampling in data

#  # fit model
#  out = tinyVAST( data = Data,
#             formula = Count ~ 1,
#             spatial_graph = graph,
#             space_column = c("X","Y"),
#             variable_column = "var",
#             time_column = "time",
#             distribution_column = "dist",
#             sem = "" )
#
#  #
#  integrate_output( out,
#                    newdata = Data,
#                    bias.correct = TRUE )
#  integrate_output( out,
#                    newdata = Data,
#                    bias.correct = FALSE,
#                    apply.epsilon = TRUE,
#                    intern = TRUE )
})
