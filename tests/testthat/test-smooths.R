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
  m_v <- fit(formula = y ~ s(x2) + s(fac, bs = "re"), data = dat)

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

  # m_g4 <- gamm4::gamm4(y ~ s(x2), data = dat, REML = FALSE)
  # m_gt <- glmmTMB::glmmTMB(y ~ s(x2), data = dat)
  m_m <- mgcv::gam(y ~ s(x2), data = dat)
  m_s <- sdmTMB::sdmTMB(y ~ s(x2), data = dat, spatial = "off")
  m_v <- fit(formula = y ~ s(x2), data = dat)

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

test_that("2D splines work (or throw an error for now)", {
  set.seed(1)
  dat <- mgcv::gamSim(1, n = 400, scale = 2)
  dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
  dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5
  dat$y <- as.numeric(dat$y)

  m_m <- mgcv::gam(y ~ t2(x1, x2), data = dat)
  m_s <- sdmTMB::sdmTMB(y ~ t2(x1, x2), data = dat, spatial = "off")
  expect_error(m_v <- fit(formula = y ~ t2(x1, x2), data = dat), regexp = "t2")
  expect_error(m_v <- fit(formula = y ~ te(x1, x2), data = dat), regexp = "te")
})

