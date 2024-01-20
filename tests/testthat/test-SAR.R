


test_that("SAR model with spatio-temporal dynamics works", {
  library(igraph)
  options("tinyVAST.verbose" = FALSE)
  data( salmon_returns )

  # Transform data
  salmon_returns$Biomass_nozeros = ifelse( salmon_returns$Biomass==0,
                                           NA, salmon_returns$Biomass )
  Data = na.omit(salmon_returns)

  # Define graph for SAR process
  unconnected_graph = make_empty_graph( nlevels(Data$Region) )
  V(unconnected_graph)$name = levels(Data$Region)

  # Define SEM for AR2 process
  dsem = "
    sockeye -> sockeye, -1, lag1_sockeye
    sockeye -> sockeye, -2, lag2_sockeye

    pink -> pink, -1, lag1_pink
    pink -> pink, -2, lag2_pink

    chum -> chum, -1, lag1_chum
    chum -> chum, -2, lag2_chum
  "

  # Fit tinyVAST model
  mytiny0 = tinyVAST(
       formula = Biomass_nozeros ~ 0 + Species + Region,
       data = Data,
       dsem = dsem,
       variable_column = "Species",
       time_column = "Year",
       space_column = "Region",
       distribution_column = "Species",
       family = list( "chum" = lognormal(),
                            "pink" = lognormal(),
                            "sockeye" = lognormal() ),
       spatial_graph = unconnected_graph,
       control = tinyVASTcontrol( profile="alpha_j" ) )

  # Summarize output
  Summary = summary(mytiny0, what="dsem")
})
