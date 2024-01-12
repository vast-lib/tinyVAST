
test_that("deviance residuals for Gamma match glm", {
  skip_on_cran()
  skip_on_ci()

  set.seed(101)
  x = rnorm(100)
  mu = 1 + 0.5*x
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

