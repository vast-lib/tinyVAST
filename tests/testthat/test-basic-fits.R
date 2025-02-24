sim_dat1 <- function(seed = 192838) {
  set.seed(seed)
  n_x <- n_y <- 25
  n_w <- 10
  R_xx <- exp(-0.4 * abs(outer(1:n_x, 1:n_x, FUN = "-")))
  R_yy <- exp(-0.4 * abs(outer(1:n_y, 1:n_y, FUN = "-")))
  z <- mvtnorm::rmvnorm(1, sigma = kronecker(R_xx, R_yy))
  # Simulate nuisance parameter z from oscillatory (day-night) process
  w <- sample(1:n_w, replace = TRUE, size = length(z))
  Data <- data.frame(expand.grid(x = 1:n_x, y = 1:n_y), w = w, z = as.vector(z) + cos(w / n_w * 2 * pi))
  Data$n <- Data$z + rnorm(nrow(Data), sd = 1)
  # Add columns for multivariate and temporal dimensions
  Data$var <- "n"
  Data$time <- 2020
  Data
}
dat <- sim_dat1()

test_that("Basic tinyVAST works", {
  mesh <- fmesher::fm_mesh_2d(dat[, c("x", "y")], n = 100)
  out <- tinyVAST(
    data = dat,
    formula = n ~ s(w),
    spatial_domain = mesh,
    space_term = "",
    control = tinyVASTcontrol(
      verbose = TRUE,
      newton_loops = 1,
      silent = FALSE
    )
  )
  expect_s3_class(out, "tinyVAST")
  s <- summary(out$sdrep)
  expect_true(sum(is.na(s[,2])) == 0L)
})

# test_that("data_colnames are robust", {
#
#   expect_error({out <- tinyVAST(
#     data = dat,
#     formula = n ~ s(w),
#     spatial_domain = mesh,
#     control = tinyVASTcontrol(quiet = TRUE, trace = 0),
#     data_colnames = list(
#       space = c("x", "y"),
#       variable = "banana",
#       time = "time",
#       distribution = "dist"
#     ),
#     space_term = ""
#   )
#   }, regexp = "variable")
  #
#  mesh <- fmesher::fm_mesh_2d(dat[, c("x", "y")], n = 100)
#  expect_error({out <- tinyVAST(
#    data = dat,
#    formula = n ~ s(w),
#    spatial_domain = mesh,
#    space_columns = c("x", "y"),
#    variable_column = "var",
#    time_column = "time",
#    distribution_column = "dist",
#    space_term = ""
#  )
#  }, regexp = "data_colnames")
# })
