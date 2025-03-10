
#' @title print summary of tinyVAST model
#' @param x output from \code{tinyVAST}
#' @param ... not used
#' @method print tinyVAST
#' @import methods
#' @return invisibly returns a named list of key model outputs and summary
#'  statements
#' @export
print.tinyVAST <-
function( x,
          ... ){
  out = NULL

  #
  cat( "Call: \n")
  out$call = print(x$call)
  cat( "\n")

  cat( "Run time: \n")
  out$run_time = print(x$run_time)
  cat( "\n")

  cat( "Family: \n")
  out$family = print(x$internal$family)
  cat( "\n")

  if( !is.null(x$sdrep) ){
    cat( "\n")
    out$sdrep = print(x$sdrep)
    cat( "\n")
  }
  if( !is.null(x$deviance_explained) ){
    cat( "Proportion conditional deviance explained: \n")
    out$deviance_explained = print(x$deviance_explained)
    cat( "\n")
  }

  out$space_term = summary( x, "space_term")
  if( nrow(out$space_term)>0 ){
    cat( "space_term: \n")
    print(out$space_term)
    cat( "\n")
  }

  out$time_term = summary( x, "time_term")
  if( nrow(out$time_term)>0 ){
    cat( "time_term: \n")
    print(out$time_term)
    cat( "\n")
  }

  out$spacetime_term = summary( x, "spacetime_term")
  if( nrow(out$spacetime_term)>0 ){
    cat( "spacetime_term: \n")
    print(out$spacetime_term)
    cat( "\n")
  }

  out$fixed = summary( x, "fixed")
  if( nrow(out$fixed)>0 ){
    cat( "Fixed terms: \n")
    print(out$fixed)
    cat( "\n")
  }
  return(invisible(out))
}


#' @title summarize tinyVAST
#'
#' @description summarize parameters from a fitted tinyVAST
#'
#' @details
#' \code{tinyVAST} includes three components:
#' \describe{
#'   \item{Space-variable interaction}{a separable Gaussian Markov random field (GMRF)
#'   constructed from a structural equation model (SEM) and a spatial variable}
#'   \item{Space-variable-time interaction}{a separable GMRF constructed from a
#'   a dynamic SEM (a nonseparable time-variable interaction) and a spatial variable}
#'   \item{Additive variation}{a generalized additive model (GAM), representing exogenous
#'   covariates }
#' }
#' Each of these are summarized and interpreted differently, and \code{summary.tinyVAST}
#' facilitates this.
#'
#' Regarding the DSEM componennt, tinyVAST includes an "arrow and lag"
#' notation, which specifies the set of
#' path coefficients and exogenous variance parameters to be estimated. Function \code{tinyVAST}
#' then estimates the maximum likelihood value for those coefficients and parameters
#' by maximizing the log-marginal likelihood.
#'
#' However, many users will want to associate individual parameters and standard errors
#' with the path coefficients that were specified using the "arrow and lag" notation.
#' This task is complicated in
#' models where some path coefficients or variance parameters are specified to share a single value a priori,
#' or were assigned a name of NA and hence assumed to have a fixed value a priori (such that
#' these coefficients or parameters have an assigned value but no standard error).
#' The \code{summary} function therefore compiles the MLE for coefficients (including duplicating
#' values for any path coefficients that assigned the same value) and standard error
#' estimates, and outputs those in a table that associates them with the user-supplied path and parameter names.
#' It also outputs the z-score and a p-value arising from a two-sided Wald test (i.e.
#' comparing the estimate divided by standard error against a standard normal distribution).
#'
#' @param object Output from [tinyVAST()]
#' @param what What component to summarize, whether \code{space_term}, \code{spacetime_term}, or
#'        \code{fixed} for the fixed effects included in the GAM formula
#' @param predictor whether to get the 1st or 2nd linear predictor (the latter is
#'        only applicable in delta models)
#' @param ... Not used
#'
#' @return
#' A data-frame containing the estimate (and standard errors, two-sided Wald-test
#' z-value, and associated p-value if the standard errors are available) for
#' model parameters, including the fixed-effects specified via \code{formula},
#' or the path coefficients for the spatial SEM specified via \code{space_term},
#' the dynamic SEM specified via \code{time_term}, or the spatial dynamic SEM
#' specified via \code{spacetime_term}
#'
#' @method summary tinyVAST
#' @export
summary.tinyVAST <-
function( object,
          what = c("space_term","time_term","spacetime_term","fixed"),
          predictor = c( "one", "two" ),
          ... ){

  #
  what = match.arg(what)
  predictor = match.arg(predictor)

  ParHat = object$obj$env$parList()
  if( !is.null(object$sdrep) ){
    SE = as.list( object$sdrep, report=FALSE, what="Std. Error")
  }

  # SEM component
  if( what=="space_term" ){
    if(predictor=="one"){
      term = "space_term_ram"
      parname = "theta_z"
    }else{
      term = "delta_space_term_ram"
      parname = "theta2_z"
    }

    model = object$internal[[term]]$output$ram
    if(nrow(model)>0){
      model$to = as.character(object$internal$variables)[model$to]
      model$from = as.character(object$internal$variables)[model$from]
    }

    #
    coefs = data.frame( model, "Estimate"=c(NA,ParHat[[parname]])[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,'start']), coefs$Estimate )
    if( !is.null(object$sdrep) ){
      coefs = data.frame( coefs, "Std_Error"=c(NA,SE[[parname]])[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
  }

  # DSEM component
  if( what=="time_term" ){
    if(predictor=="one"){
      term = "time_term_ram"
      parname = "nu_z"
    }else{
      term = "delta_time_term_ram"
      parname = "nu2_z"
    }

    model = object$internal[[term]]$output$model
    model = data.frame( heads = model[,'direction'],
                        to = model[,'second'],
                        from = model[,'first'],
                        parameter = model[,'parameter'],
                        start = model[,'start'],
                        lag = model[,'lag'] )

    #
    coefs = data.frame( model, "Estimate"=c(NA,ParHat[[parname]])[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,'start']), coefs$Estimate )
    if( !is.null(object$sdrep) ){
      coefs = data.frame( coefs, "Std_Error"=c(NA,SE[[parname]])[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
  }

  # DSEM component
  if( what=="spacetime_term" ){
    if(predictor=="one"){
      term = "spacetime_term_ram"
      parname = "beta_z"
    }else{
      term = "delta_spacetime_term_ram"
      parname = "beta2_z"
    }

    if( is(object$internal[[term]]$output,"dsem_ram") ){
      model = object$internal[[term]]$output$model
      model = data.frame( heads = model[,'direction'],
                          to = model[,'second'],
                          from = model[,'first'],
                          parameter = model[,'parameter'],
                          start = model[,'start'],
                          lag = model[,'lag'] )
    }else if( is(object$internal[[term]]$output,"eof_ram") ){
      model = object$internal[[term]]$output$model
      vars = object$internal[[term]]$output$variances
      model = data.frame( heads = c( rep(1,nrow(model)), rep(2,nrow(vars)) ),
                          to = c( as.character(model$to), as.character(vars$to) ),
                          from = c( as.character(model$from), as.character(vars$from) ),
                          parameter = c( model$parameter, vars$parameter ),
                          start = NA,
                          lag = NA )
    }else{
      stop("Class for what=`spacetime_term` not implemented for `summary(.)`")
    }

    #
    coefs = data.frame( model, "Estimate"=c(NA,ParHat[[parname]])[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,'start']), coefs$Estimate )
    if( !is.null(object$sdrep) ){
      coefs = data.frame( coefs, "Std_Error"=c(NA,SE[[parname]])[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
  }

  if( what=="fixed" ){
    if(predictor=="one"){
      term = "X_ij"
      parname = "alpha_j"
    }else{
      term = "X2_ij"
      parname = "alpha2_j"
    }

    coefs = data.frame(
      Estimate = ParHat[[parname]]
    )
    if( !is.null(object$sdrep) ){
      coefs = data.frame( coefs, "Std_Error"=SE[[parname]] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
    rownames(coefs) = colnames(object$tmb_input$tmb_data[[term]])
  }

  return(coefs)
}

#' Calculate residuals
#'
#' @title Calculate deviance or response residuals for tinyVAST
#'
#' @param object Output from [tinyVAST()]
#' @param type which type of residuals to compute (only option is `"deviance"` or `"response"` for now)
#' @param ... Note used
#'
#' @method residuals tinyVAST
#' @return a vector residuals, associated with each row of \code{data} supplied during fitting
#'
#' @export
residuals.tinyVAST <-
function( object,
          type = c("deviance","response"),
          ... ){

  # https://stats.stackexchange.com/questions/1432/what-do-the-residuals-in-a-logistic-regression-mean
  # Normal deviance residuals
  #if( FALSE ){
  #  x = rnorm(10)
  #  y = x + rnorm(10)
  #  Glm = glm( y ~ 1 + x, family="gaussian")
  #  mu = predict(Glm,type="response")
  #  r1 = y - mu
  #  r1 - resid(Glm)
  #}
  ## Poisson deviance residuals
  # Poisson: sign(y - mu) * sqrt(2*(ifelse(y==0, 0, y*log(y/mu)) - (y-mu)))
  #if( FALSE ){
  #  x = rnorm(10)
  #  y = rpois(10, exp(x))
  #  Glm = glm( y ~ 1 + x, family="poisson")
  #  mu = predict(Glm,type="response")
  #  # https://stats.stackexchange.com/questions/398098/formula-for-deviance-residuals-for-poisson-model-with-identity-link-function
  #  r1 = sign(y - mu) * sqrt(2*(y*log((y+1e-10)/mu) - (y-mu)))
  #  r1 - resid(Glm)
  #}
  ## Binomial deviance residuals
  # Binomial:  -2 * ((1-y)*log(1-mu) + y*log(mu))
  #if( FALSE ){
  #  p = 0.5
  #  y = rbinom(10, prob=p, size=1)
  #  Glm = glm( y ~ 1, family="binomial")
  #  mu = predict(Glm, type="response")
  #  r1 = sign(y - mu) * sqrt(-2*(((1-y)*log(1-mu) + y*log(mu))))
  #  r1 - resid(Glm)
  #}
  ## Gamma deviance residuals
  # Gamma: 2 * ( (y-mu)/mu - log(y/mu) )
  #if( FALSE ){
  #  mu = 1
  #  cv = 0.8
  #  y = rgamma( n=10, shape=1/cv^2, scale=mu*cv^2 )
  #  Glm = glm( y ~ 1, family=Gamma(link='log'))
  #  mu = predict(Glm, type="response")
  #  r1 = sign(y - mu) * sqrt(2 * ( (y-mu)/mu - log(y/mu) ))
  #  r1 - resid(Glm)
  #}

  # Easy of use
  mu = object$rep$mu
  Y = object$tmb_inputs$tmb_data$Y
  #familycode_j = object$tmb_inputs$data$familycode_j
  report = object$rep

  #
  type = match.arg(type)
  if( type == "deviance" ){
    resid = report$devresid
  }
  if( type == "response" ){
    resid = Y - mu
  }

  return(resid)
}

#' Extract the (marginal) log-likelihood of a tinyVAST model
#'
#' @param object output from \code{tinyVAST}
#' @param ... not used
#'
#' @return object of class \code{logLik} with attributes
#'   \item{val}{log-likelihood}
#'   \item{df}{number of parameters}
#' @importFrom stats logLik
#' @export
logLik.tinyVAST <- function(object, ...) {
  val = -1 * object$opt$objective
  # Get df including profiled parameters
  df = length( object$opt$par ) +
       sum(names(object$obj$env$last.par) %in% object$internal$control$profile)
  # S3 object "logLik"
  out = structure( val,
             df = df,
             class = "logLik")
  return(out)
}

#' Get fitted values from a tinyVAST model
#'
#' @param object The fitted tinyVAST model object
#' @param ... Not used
#' @importFrom stats fitted
#' @export
#' @return a vector of fitted values for each row of \code{data}
#' @noRd
fitted.tinyVAST <- function(object, ...) {
  predict( object, what="mu_g" )
}

#' Extract Variance-Covariance Matrix
#'
#' extract the covariance of fixed effects, or both fixed and random effects.
#'
#' @param object output from [tinyVAST()]
#' @param which whether to extract the covariance among fixed effects, random effects, or both
#' @param ... ignored, for method compatibility
#' @importFrom stats vcov
#' @method vcov tinyVAST
#' @return
#' A square matrix containing the estimated covariances among the parameter estimates in the model.
#' The dimensions dependend upon the argument \code{which}, to determine whether fixed, random effects,
#' or both are outputted.
#'
#' @export
vcov.tinyVAST <-
function( object,
          which = c("fixed", "random", "both"),
          ...) {

  which = match.arg(which)

  if( which=="fixed" ){
    V = object$sdrep$cov.fixed
    if(is.null(V)){
      warning("Please re-run `tinyVAST` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$sdrep$jointPrecision
    if(is.null(H)){
      warning("Please re-run `tinyVAST` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

#' Simulate new data from a fitted model
#'
#' \code{simulate.tinyVAST} is an S3 method for producing a matrix of simulations from
#' a fitted model. It can be used with the \pkg{DHARMa} package
#' among other uses.  Code is modified from the version in sdmTMB
#'
#' @param object output from [tinyVAST()]
#' @param nsim how many simulations to do
#' @param seed random seed
#' @param type How parameters should be treated. `"mle-eb"`: fixed effects
#'   are at their maximum likelihood (MLE) estimates  and random effects are at
#'   their empirical Bayes (EB) estimates. `"mle-mvn"`: fixed effects are at
#'   their MLEs but random effects are taken from a single approximate sample.
#'   This latter option is a suggested approach if these simulations will be
#'   used for goodness of fit testing (e.g., with the DHARMa package).
#' @param ... not used
#'
#' @return
#' A matrix with row for each row of \code{data} in the fitted model and \code{nsim}
#' columns, containing new samples from the fitted model.
#'
#' @method simulate tinyVAST
#' @importFrom stats simulate
#'
#' @examples
#' set.seed(101)
#' x = seq(0, 2*pi, length=100)
#' y = sin(x) + 0.1*rnorm(length(x))
#' fit = tinyVAST( data=data.frame(x=x,y=y), formula = y ~ s(x) )
#' simulate(fit, nsim=100, type="mle-mvn")
#'
#' if(requireNamespace("DHARMa")){
#'   # simulate new data conditional on fixed effects
#'   # and sampling random effects from their predictive distribution
#'   y_iz = simulate(fit, nsim=500, type="mle-mvn")
#'
#'   # Visualize using DHARMa
#'   res = DHARMa::createDHARMa( simulatedResponse = y_iz,
#'                       observedResponse = y,
#'                       fittedPredictedResponse = fitted(fit) )
#'   plot(res)
#' }
#' @export
simulate.tinyVAST <-
function( object,
          nsim = 1L,
          seed = sample.int(1e+06, 1L),
          type = c("mle-eb", "mle-mvn"),
          ... ) {

  type = match.arg(type)
  set.seed(seed)
  par_iz = outer( object$obj$env$last.par.best, rep(1,nsim) )

  if( (type == "mle-mvn") & (length(object$obj$env$random)>0) ){
    tmp = object$obj$env$MC( n = nsim, keep = TRUE, antithetic = FALSE )
    par_iz[object$obj$env$lrandom(),] = attr(tmp, "samples")
  }

  y_iz = apply( par_iz,
                MARGIN = 2,
                FUN = function(par_i) object$obj$simulate(par_i)$y_i )
  return( y_iz )
}

