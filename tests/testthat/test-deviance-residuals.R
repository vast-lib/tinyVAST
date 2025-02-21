
test_that("deviance residuals for Gamma match glm", {
  skip_on_cran()
  #skip_on_ci()

  set.seed(101)
  x = rnorm(100)
  mu = exp(1 + 0.5*x)
  y = rgamma( n=100, shape=1/0.5^2, scale=mu*0.5^2 )

  # simulated model
  mytiny = tinyVAST( y ~ 1 + x,
            data = data.frame(y=y, x=x),
            family = Gamma(link = "log") )
  resid1 = residuals(mytiny, type="deviance")

  # Null model
  mytiny0 = tinyVAST( y ~ 1,
            data = data.frame(y=y),
            family = Gamma(link = "log") )
  resid0 = residuals(mytiny0, type="deviance")

  #
  myglm = glm( y ~ 1 + x, family=Gamma(link="log"))
  resid2 = residuals( myglm, type="deviance" )
  expect_equal( as.numeric(resid1), as.numeric(resid2),
                tolerance=1e-3 )

  # Compare percent-deviance-explained
  PDEglm = with(summary(myglm), 1 - deviance/null.deviance)
  PDEtiny = (sum(resid0^2)-sum(resid1^2)) / sum(resid0^2)
  expect_equal( PDEglm, PDEtiny,
                tolerance=1e-3 )
})

test_that("nbinom1 and nbinom2 family matches glmmTMB", {
  skip_on_cran()
  #skip_on_ci()
  library(glmmTMB)

  set.seed(101)
  X = rnorm(100)
  eps = 0.5 * rnorm(100)
  p = 2 + 1 * X + eps
  Y = rpois( length(X), lambda = exp(p) )
  data = data.frame(X=X, Y=Y)

  #
  glmmTMB1 = glmmTMB( Y ~ 1 + X, family = nbinom1(), data = data )
  tinyVAST1 = tinyVAST( Y ~ 1 + X, family = nbinom1(), data = data )
  expect_equal( as.numeric(glmmTMB1$fit$par), as.numeric(tinyVAST1$opt$par),
                tolerance=1e-3 )

  #
  glmmTMB2 = glmmTMB( Y ~ 1 + X, family = nbinom2(), data = data )
  tinyVAST2 = tinyVAST( Y ~ 1 + X, family = nbinom2(), data = data )
  expect_equal( as.numeric(glmmTMB2$fit$par), as.numeric(tinyVAST2$opt$par),
                tolerance=1e-3 )
  tmp = mgcv::gam( Y ~ 1 + X, family = nb() )
  if( FALSE ){
    y = Y
    mu = predict(tinyVAST2)
    theta = exp(tinyVAST2$opt$par[3])
    c1 = y * log(y / mu)
    c2 = ( y + 1.0/theta ) * log( (1.0 + theta*y) / (1.0 + theta*mu) )
    pow = function(a,b) a^b
    deviance = 2.0 * (c1 - c2)
    devresid = sign( y - mu ) * pow( deviance, 0.5 )
  }

})

test_that("deviance residuals for lognormal match glm", {
  skip_on_cran()
  #skip_on_ci()

  set.seed(101)
  x = rnorm(100)
  logmu = 1 + 0.5*x
  y = rlnorm( n=100, meanlog=logmu, sdlog=0.5 )

  # simulated model
  mytiny = tinyVAST( y ~ 1 + x,
            data = data.frame(y=y, x=x),
            family = lognormal(link = "log") )
  resid1 = residuals(mytiny, type="deviance")

  # Null model
  mytiny0 = tinyVAST( y ~ 1,
            data = data.frame(y=y),
            family = lognormal(link = "log") )
  resid0 = residuals(mytiny0, type="deviance")

  #
  myglm = glm( log(y) ~ 1 + x, family=gaussian(link="identity"))
  resid2 = residuals( myglm, type="deviance" )
  expect_equal( as.numeric(resid1), as.numeric(resid2),
                tolerance=1e-3 )

  # Compare percent-deviance-explained
  PDEglm = with(summary(myglm), 1 - deviance/null.deviance)
  PDEtiny = (sum(resid0^2)-sum(resid1^2)) / sum(resid0^2)
  expect_equal( PDEglm, PDEtiny,
                tolerance=1e-3 )
})

test_that("deviance residuals for tweedie match mgcv", {
  skip_on_cran()
  #skip_on_ci()
  skip_if_not_installed("tweedie")
  skip_if_not_installed("mgcv")

  set.seed(101)
  y = tweedie::rtweedie( n=100, mu=2, phi=1, power=1.5 )

  #
  mytiny = tinyVAST( y ~ 1,
            data = data.frame(y=y),
            family = tweedie(link = "log") )
  resid1 = residuals(mytiny, type="deviance")

  #
  library(mgcv)
  mygam = gam( y ~ 1, family=tw(link="log"))
  resid2 = residuals( mygam, type="deviance" )
  expect_equal( as.numeric(resid1), as.numeric(resid2),
                tolerance=1e-3 )
})

test_that("deviance residuals for poisson works", {
  skip_on_cran()
  #skip_on_ci()

  # simulate data
  set.seed(101)
  x = rnorm(100)
  mu = exp(1 + 0.5*x)
  y = rpois( n=100, lambda=mu )

  # delta-gamma in tinyVAST
  mytiny = tinyVAST( y ~ 1 + x,
                     data = data.frame(x=x, y=y),
                     family = poisson("log") )
  resid1 = residuals( mytiny, type="deviance" )
  # resid1^2
  # mu = mytiny$rep$mu_i; 
  # 2*y*log(y/mu) - (y-mu)
  # sign(y - mu) * pow(2*(y*log((1e-10 + y)/mu) - (y-mu)), 0.5)
  # resid1
  # pow = function(a,b) a^b
  # Type = function(z)z
  # sign(y - mu) * pow(2*(y*log((1e-10 + y)/mu) - (y-mu)), 0.5)
  
  #
  myglm = glm( y ~ 1 + x,
                     data = data.frame(x=x, y=y),
                     family = poisson("log") )
  resid2 = residuals( myglm, type="deviance" )
  expect_equal( as.numeric(resid1), as.numeric(resid2),
                tolerance=1e-3 )
})

test_that("deviance residuals for Bernoulli works", {
  skip_on_cran()
  #skip_on_ci()

  # simulate data
  set.seed(101)
  x = rnorm(100)
  mu = plogis(1 + 0.5*x)
  y = rbinom( n=100, prob=mu, size=1 )

  # delta-gamma in tinyVAST
  mytiny = tinyVAST( y ~ 1 + x,
                     data = data.frame(x=x, y=y),
                     family = binomial("logit") )
  resid1 = residuals( mytiny, type="deviance" )

  #
  myglm = glm( y ~ 1 + x,
                     data = data.frame(x=x, y=y),
                     family = binomial("logit") )
  resid2 = residuals( myglm, type="deviance" )
  expect_equal( as.numeric(resid1), as.numeric(resid2),
                tolerance=1e-3 )
})

test_that("delta-gamma works", {
  skip_on_cran()
  #skip_on_ci()

  # simulate data
  set.seed(101)
  x = rnorm(100)
  mu = exp(1 + 0.5*x)
  y = rgamma( n=100, shape=1/0.5^2, scale=mu*0.5^2 ) *
      rbinom( n=100, size=1, prob=0.5 )

  # delta-gamma in tinyVAST
  mytiny = tinyVAST( y ~ 1 + x,
                     data = data.frame(x=x, y=y),
                     family = delta_gamma(),
                     delta_options = list(
                       delta_formula = ~ 1 + x
                     ) )

  # separate GLMs
  myglm1 = glm( present ~ 1 + x,
                data = data.frame(x=x, present=y>0),
                family = binomial(link="logit") )
  myglm2 = glm( y ~ 1 + x,
                data = subset(data.frame(x=x, y=y),y>0),
                family = Gamma(link="log") )
  
  # relative deviance
  expect_equal( mytiny$rep$deviance, 
                as.numeric(myglm1$deviance + myglm2$deviance),
                tolerance=1e-2 )

  # Compare them
  glmpar = c( coef(myglm1), coef(myglm2) )
  tinypar = mytiny$opt$par[1:4]
  expect_equal( as.numeric(glmpar), as.numeric(tinypar),
                tolerance=1e-3 )
})

test_that("Poisson-link delta-gamma works", {
  skip_on_cran()
  #skip_on_ci()
  skip_if_not_installed("sdmTMB")

  # simulate data
  set.seed(101)
  x = rnorm(100)
  mu = exp(1 + 0.5*x)
  y = rgamma( n=100, shape=1/0.5^2, scale=mu*0.5^2 ) *
      rbinom( n=100, size=1, prob=0.5 )

  # delta-gamma in tinyVAST
  mytiny = tinyVAST( y ~ 1 + x,
                     data = data.frame(x=x, y=y),
                     family = delta_gamma(type="poisson-link"),
                     delta_options = list(
                       delta_formula = ~ 1 + x
                     ),
                     control = tinyVASTcontrol(verbose=FALSE) )

  # delta-gamma in sdmTMB
  mysdmTMB = sdmTMB::sdmTMB( list( y ~ 1 + x, y ~ 1 + x ),
                     data = data.frame(x=x, y=y),
                     family = delta_gamma(type="poisson-link"),
                     spatial = "off" )

  # separate GLMs
  myglm1 = glm( present ~ 1 + x,
                data = data.frame(x=x, present=y>0),
                family = binomial(link="cloglog") )

  # Compare them
  glmpar = c( coef(myglm1) )
  tinypar = mytiny$opt$par[1:4]
  sdmTMBpar = mysdmTMB$model$par[1:4]
  expect_equal( as.numeric(tinypar), as.numeric(sdmTMBpar),
                tolerance=1e-3 )
})

