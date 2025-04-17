test_that("Basic splines work", {
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  dat <- mgcv::gamSim(1, n = 400, scale = 2)
  dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
  dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
  dat$y <- as.numeric(dat$y)

  # m_g4 <- gamm4::gamm4(y ~ s(x2), data = dat, random = ~ (1 | fac), REML = FALSE)
  # m_gt <- glmmTMB::glmmTMB(y ~ s(x2) + (1 | fac), data = dat, REML = FALSE)
  m_m <- mgcv::gam(y ~ s(x2) + s(fac, bs = "re"), data = dat, REML = FALSE)
  m_s <- sdmTMB::sdmTMB(y ~ s(x2) + (1 | fac), data = dat, reml = FALSE, spatial = "off")
  m_v <- tinyVAST(formula = y ~ s(x2) + s(fac, bs = "re"), data = dat)
  m_wrong <- tinyVAST(formula = y ~ s(x2), data = dat)

  # logLik(m_g4$mer)
  # logLik(m_gt)
  logLik(m_s)
  logLik(m_m)
  logLik(m_v)

  # p_g4 <- predict(m_g4$gam)
  # p_gt <- predict(m_gt)
  p_m <- predict(m_m)
  p_s <- predict(m_s)$est
  p_v <- predict(m_v)
  p_wrong <- predict(m_wrong)
  expect_gt(cor(p_m, p_v), 0.999)
  expect_gt(cor(p_s, p_v), 0.999)
  # cor( p_v, p_wrong )

  # m_g4 <- gamm4::gamm4(y ~ s(x2), data = dat, REML = FALSE)
  # m_gt <- glmmTMB::glmmTMB(y ~ s(x2), data = dat)
  m_m <- mgcv::gam(y ~ s(x2), data = dat)
  m_s <- sdmTMB::sdmTMB(y ~ s(x2), data = dat, spatial = "off")
  m_v <- tinyVAST(formula = y ~ s(x2), data = dat)

  # logLik(m_g4$mer)
  # logLik(m_gt)
  logLik(m_s)
  logLik(m_m)
  logLik(m_v)

  # p_g4 <- predict(m_g4$gam)
  # p_gt <- predict(m_gt)
  p_m <- predict(m_m)
  p_s <- predict(m_s)$est
  p_v <- predict(m_v)
  expect_gt(cor(p_m, p_v), 0.999)
  expect_gt(cor(p_s, p_v), 0.999)
})

test_that("t2 splines throw an error for now", {
  set.seed(1)
  dat <- mgcv::gamSim(1, n = 400, scale = 2)
  dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
  dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
  dat$y <- as.numeric(dat$y)

  m_m <- mgcv::gam(y ~ t2(x1, x2), data = dat)
  m_s <- sdmTMB::sdmTMB(y ~ t2(x1, x2), data = dat, spatial = "off")
  expect_error(m_v <- tinyVAST(formula = y ~ t2(x1, x2), data = dat), regexp = "t2")
})

test_that("ti and te splines work", {
  set.seed(1)
  dat <- mgcv::gamSim(1, n = 400, scale = 2)
  dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
  dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
  dat$y <- as.numeric(dat$y)

  form = y ~ s(x1) + s(x2) + ti(x1, x2)
  m_m <- mgcv::gam( form, data = dat, method="ML")
  m_v <- tinyVAST( form, data = dat)
  p_m <- predict(m_m)
  p_v <- predict(m_v)
  expect_gt(cor(p_v, p_m), 0.98)
  cbind( "tinyVAST" = m_v$internal$parlist$log_lambda, "mgcv" = log(m_m$sp) )

  #
  form = y ~ te(x1, x2)
  m_m <- mgcv::gam( form, data = dat, method="ML")
  m_v <- tinyVAST( form, data = dat)
  p_m <- predict(m_m)
  p_v <- predict(m_v)
  expect_gt(cor(p_v, p_m), 0.98)
  cbind( "tinyVAST" = m_v$internal$parlist$log_lambda, "mgcv" = log(m_m$sp) )
})

# modified from sdmTMB tests:
test_that("A model with s(x, bs = 'fs') works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  d = as.data.frame(d)
  d$yearf <- as.factor(d$year)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, by = year, bs = "fs"), data = d, method = "REML")
  m_v <- tinyVAST(
    data = d,
    formula = log(density) ~ s(depth_scaled, by = year, bs = "fs")
  )
  p2 <- predict(m_mgcv)
  p3 <- predict(m_v)
  expect_gt(stats::cor(p3, p2), 0.995)
})

test_that("An fx=TRUE smoother errors out", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  expect_error(m_v <- tinyVAST(
    data = d,
    formula = log(density) ~ s(depth_scaled, fx = TRUE)
  )) # TODO issue a meaningful error!
})

test_that("A model with s(x, bs = 'gp') works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  d = as.data.frame(d)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "gp"), data = d, method = "REML")
  m_v <- tinyVAST(log(density) ~ s(depth_scaled, bs = "gp"), data = d, method = "REML")
  m_s <- sdmTMB::sdmTMB(log(density) ~ s(depth_scaled, bs = "gp"), data = d, spatial = "off")
  p <- predict(m_mgcv)
  p_s <- predict(m_s)$est
  p_v <- predict(m_v)
  plot(p_v, p)
  plot(p_v, p_s)
  expect_gt(stats::cor(p_s, p_v), 0.999)
  expect_gt(stats::cor(p_v, p), 0.999)
})

test_that("A model with s(x, bs = 'ds') works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  d = as.data.frame(d)
  m_s <- sdmTMB::sdmTMB(
    log(density) ~ s(depth_scaled, bs = "ds"),
    data = d, spatial = "off"
  )
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "ds"), data = d)
  m_v <- tinyVAST(formula = log(density) ~ s(depth_scaled, bs = "ds"), data = d)
  p_m <- predict(m_mgcv)
  p_v <- predict(m_v)
  p_s <- predict(m_s)$est
  plot(p_s, p_m)
  plot(p_s, p_v)
  expect_gt(stats::cor(p_s, p_v), 0.999)
  expect_gt(stats::cor(p_m, p_v), 0.999)
})

test_that("A model with s(x, bs = 'cr') works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  d = as.data.frame(d)
  m_s <- sdmTMB::sdmTMB(log(density) ~ s(depth_scaled, bs = "cr"), data = d, spatial = "off")
  m_m <- mgcv::gam(data = d, formula = log(density) ~ s(depth_scaled, bs = "cr"))
  m_v <- tinyVAST(data = d, formula = log(density) ~ s(depth_scaled, bs = "cr"))
  p_v <- predict(m_v)
  p_m <- predict(m_m)
  p_s <- predict(m_s)$est
  plot(p_v, p_s)
  plot(p_v, p_m)
  expect_gt(stats::cor(p_s, p_v), 0.999)
  expect_gt(stats::cor(p_v, p_m), 0.999)
})

test_that("A model with s(x, bs = 'cs') works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  d = as.data.frame(d)
  m_s <- sdmTMB::sdmTMB(log(density) ~ s(depth_scaled, bs = "cs"), data = d, spatial = "off")
  m_m <- mgcv::gam(data = d, formula = log(density) ~ s(depth_scaled, bs = "cs"))
  m_v <- tinyVAST(data = d, formula = log(density) ~ s(depth_scaled, bs = "cs"))
  p_v <- predict(m_v)
  p_m <- predict(m_m)
  p_s <- predict(m_s)$est
  plot(p_v, p_s)
  plot(p_v, p_m)
  expect_gt(stats::cor(p_s, p_v), 0.999)
  expect_gt(stats::cor(p_v, p_m), 0.999)
})

test_that("A model with s(x, bs = 'cc') works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(sdmTMB::pcod, density > 0)
  d = as.data.frame(d)
  m_s <- sdmTMB::sdmTMB(log(density) ~ s(depth_scaled, bs = "cc"), data = d, spatial = "off")
  m_m <- mgcv::gam(data = d, formula = log(density) ~ s(depth_scaled, bs = "cc"))
  m_v <- tinyVAST(data = d, formula = log(density) ~ s(depth_scaled, bs = "cc"))
  p_v <- predict(m_v)
  p_m <- predict(m_m)
  p_s <- predict(m_s)$est
  plot(p_v, p_s)
  plot(p_v, p_m);abline(0, 1)
  expect_gt(stats::cor(p_s, p_v), 0.999)
  expect_gt(stats::cor(p_v, p_m), 0.999)
})

test_that("A model with s() by variables works", {
  set.seed(1)
  # For some reason, it doesn't work on CI but does locally
  skip_on_ci()
  skip_on_cran()

  dat <- mgcv::gamSim(4)
  m_mgcv <- mgcv::gam(y ~ fac + s(x2, by = fac) + s(x0), data = dat)
  p_mgcv <- predict(m_mgcv)

  m_s <- sdmTMB::sdmTMB(formula = y ~ fac + s(x2, by = fac) + s(x0), data = dat, spatial = "off")

  m_v <- tinyVAST(formula = y ~ fac + s(x2, by = fac) + s(x0), data = dat)
  expect_s3_class(m_v, "tinyVAST")

  p_m <- predict(m_mgcv)
  p_v <- predict(m_v)
  p_s <- predict(m_s)$est

  plot(p_m, p_v);abline(a = 0, b = 1)
  expect_gt(cor(p_v, p_m), 0.99)
})
