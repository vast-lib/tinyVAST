
test_that("deviance residuals for Gamma match glm", {
  skip_on_cran()
  skip_on_ci()

  set.seed(101)
  x = rnorm(100)
  mu = exp(1 + 0.5*x)
  y = rgamma( n=100, shape=1/0.5^2, scale=mu*0.5^2 )
  
  # simulated model
  mytiny = tinyVAST( y ~ 1 + x, 
            data = data.frame(y=y),
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

test_that("deviance residuals for tweedie match mgcv", {
  skip_on_cran()
  skip_on_ci()
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

test_that("delta-gamma works", {
  skip_on_cran()
  skip_on_ci()

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

  # Compare them
  glmpar = c( coef(myglm1), coef(myglm2) )
  tinypar = mytiny$opt$par[1:4]
  expect_equal( as.numeric(glmpar), as.numeric(tinypar),
                tolerance=1e-3 )
})

test_that("Poisson-link delta-gamma works", {
  skip_on_cran()
  skip_on_ci()
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
                     family = delta_poisson_link_gamma(),
                     delta_options = list(
                       delta_formula = ~ 1 + x
                     ),
                     control = tinyVASTcontrol(verbose=FALSE) )

  # delta-gamma in sdmTMB
  mysdmTMB = sdmTMB( list( y ~ 1 + x, y ~ 1 + x ),
                     data = data.frame(x=x, y=y),
                     family = delta_poisson_link_gamma(),
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

