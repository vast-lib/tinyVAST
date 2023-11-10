#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set, and widely used interface for objects
#'
#' @inheritParams dsem::dsem
#' @inheritParams dsem::make_ram
#'
#' @param data Data-frame of predictor and response variables.
#' @param formula Formula with response on left-hand-side and predictors on right-hand-side,
#'        parsed by \code{mgcv} and hence allowing \code{s(.)} for splines
#' @param family_link Vector of length-two, indicating the distribution and link-function
#' @param spatial_graph Object that represents spatial relationships, either using \code{fmesher}
#'        to apply the SPDE method, \code{igraph} to apply a simultaneous autoregressive (SAR)
#'        process, or \code{NULL} to specify a single site.
#' @param control Output from \code{\link{tinyVASTcontrol}}, used to define user
#'        settings, and see documentation for that function for details.
#'
#' @details
#' [tinyVAST] includes four basic inputs that specify the model structure:
#' * [formula] specifies covariates and splines in a Generalized Additive Model;
#' * [sem] specifies interactions among variables and over time
#' * [spatial_graph] specifies spatial correlations
#'
#' the default [sem=NULL] turns off all multivariate and temporal indexing, such
#' that \code{spatial_graph} is then ignored, and the model collapses
#' to a standard model using \code{mgcv::gam}.  To specify a univeriate spatial model,
#' the user must specify both \code{spatial_graph} and \code{sem=""}, where the latter
#' is then parsed to include a single exogenous variance for the single variable
#'
#' | Model type | How to specify |
#' | --- | --- |
#' | Spatial model | specify both \code{spatial_graph} and \code{sem=""}, where the latter is then parsed to include a single exogenous variance for the single variable |
#'
#' @importFrom dsem make_ram classify_variables parse_path
#' @importFrom igraph as_adjacency_matrix
#'
#' @examples
#' methods(class="tinyVAST")
#'
#' @useDynLib tinyVAST, .registration = TRUE
#' @export
fit <-
function( data,
          formula,
          sem = NULL,
          estimate_delta0 = FALSE,
          family_link = c(0,0),
          #spatial_sem,
          #spatiotemporal_sem,
          data_colnames = list("spatial"=c("x","y"), "variable"="var", "time"="time", "distribution"="dist"),
          times = sort(unique(data[,data_colnames$time])),
          variables = unique(data[,data_colnames$variable]),
          spatial_graph = NULL,
          control = tinyVASTcontrol(),
          ... ){
  # https://roxygen2.r-lib.org/articles/rd-formatting.html#tables for roxygen formatting
  start_time = Sys.time()

  # General error checks
  if( class(control) != "tinyVASTcontrol" ) stop("`control` must be made by `tinyVASTcontrol()`")

  ##############
  # SEM constructor
  ##############

  # (I-Rho)^-1 * Gamma * (I-Rho)^-1
  if( is.null(sem) ){
    ram_output = list(
      ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
      model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
    )
    times = numeric(0)
    variables = numeric(0)
    # Allow user to avoid specifying these if is.null(sem)
    if( !(data_colnames$time %in% colnames(data)) ){
      data = cbind( data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$time)) )
    }
    if( !(data_colnames$variable %in% colnames(data)) ){
      data = cbind(data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$variable)))
    }
  }else{
    ram_output = make_ram( sem, times=times, variables=variables, quiet=control$quiet, covs=variables )
  }
  ram = ram_output$ram

  # Identify arrow-type for each beta_j estimated in RAM
  which_nonzero = which(ram[,4]>0)
  beta_type = tapply( ram[which_nonzero,1], INDEX=ram[which_nonzero,4], max)

  # Error checks
  if( any((ram_output$model[,'direction']==2) & (ram_output$model[,2]!=0)) ){
    stop("All two-headed arrows should have lag=0")
  }
  if( !all(c(ram_output$model[,'first'],ram_output$model[,'second']) %in% variables) ){
    stop("Some variable in `sem` is not in `tsdata`")
  }

  ##############
  # Spatial domain constructor
  ##############

  if( "fm_mesh_2d" %in% class(spatial_graph) ){
    # SPDE
    n_s = spatial_graph$n
    spatial_method_code = 1
    spatial_list = fm_fem( spatial_graph )
    spatial_list = list("M0"=spatial_list$c0, "M1"=spatial_list$g1, "M2"=spatial_list$g2)
    A_is = fm_evaluator( spatial_graph, loc=as.matrix(data[,data_colnames$spatial]) )$proj$A
  }else if( "igraph" %in% class(spatial_graph) ) {
    # SAR
    spatial_method_code = 2
    Adj = as_adjacency_matrix( spatial_graph, sparse=TRUE )
    n_s = nrow(Adj)
    Match = match( data[,data_colnames$spatial], rownames(Adj) )
    if(any(is.na(Match))) stop("Check `spatial_graph` for SAR")
    A_is = sparseMatrix( i=1:nrow(data), j=Match, x=rep(1,nrow(data)) )
  }else {      # if( !is.null(sem) )
    # Single-site
    spatial_method_code = 3
    n_s = 1
    A_is = matrix(1, nrow=nrow(data), ncol=1)    # dgCMatrix
    A_is = as(Matrix(A_is),"dgCMatrix")
    spatial_list = list( "M0" = as(Matrix(1,nrow=1,ncol=1),"dgCMatrix"),
                         "M1" = as(Matrix(0,nrow=1,ncol=1),"dgCMatrix"),
                         "M2" = as(Matrix(0,nrow=1,ncol=1),"dgCMatrix") )
  }
  Atriplet = Matrix::mat2triplet(A_is)

  ##############
  # Formula constructor
  ##############

  # Initial constructor of splines
  gam_setup = gam( formula, data = data, fit=FALSE ) # select doesn't do anything in this setup
  y_i = model.response(gam_setup$mf)  # OR USE: model.extract(gam_setup$mf, "response")

  # Extrtact and combine penelization matrices
  S_list = lapply( seq_along(gam_setup$smooth), \(x) gam_setup$smooth[[x]]$S[[1]] )
  S_kk = .bdiag(S_list)         # join S's in sparse matrix
  Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
  if(is.null(Sdims)) Sdims = vector(length=0)

  # Get covariates
  which_se = grep( pattern="s(", x=gam_setup$term.names, fixed=TRUE )
  X_ij = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
  Z_ik = gam_setup$X[,which_se,drop=FALSE]

  ##############
  # distribution/link
  ##############

  # Check for errors, and telescope
  if( !(data_colnames$distribution %in% colnames(data)) ){
    data = cbind( data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$distribution)) )
  }
  if( (any(data[,data_colnames$distribution]>nrow(family_link))) | (any(data[,data_colnames$distribution]<1)) ){
    stop("`data[,data_colnames$distribution]` has some bad value")
  }

  # Telescope format
  if( !is.matrix(family_link) ){
    family_link = matrix( family_link, nrow=1, ncol=2 )
  }

  # Construct log_sigma based on family_link
  remove_last = \(x) x[-length(x)]
  Nsigma_e = sapply( as.character(family_link[,1]), "0"=1, "1"=2, FUN=switch )
  log_sigma = rep( 0, sum(Nsigma_e) )
  Edims_ez = cbind( "start"=remove_last(cumsum(c(0,Nsigma_e))), "length"=Nsigma_e )

  ##############
  # Build inputs
  ##############

  # Turn of t_i and c_i when times and variables are missing, so that delta_k isn't built
  if( length(times) > 0 ){
    t_i = match( data[,data_colnames$time], times )
  }else{ t_i = integer(0) }
  if( length(variables) > 0 ){
    c_i = match( data[,data_colnames$var], variables )
  }else{ c_i = integer(0) }

  # Drop rows from Aistc_zz and Axi_z (e.g., when times and/or variables are empty because sem=NULL)
  Aistc_zz = cbind(Atriplet$i, Atriplet$j, t_i[Atriplet$i], c_i[Atriplet$i]) - 1
  which_Arows = which(apply( Aistc_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) ))

  # make dat
  tmb_data = list(
    spatial_method_code = spatial_method_code,
    y_i = y_i,
    X_ij = X_ij,
    Z_ik = Z_ik,
    t_i = t_i - 1, # -1 to convert to CPP index
    c_i = c_i - 1, # -1 to convert to CPP index
    f_ez = family_link,
    e_i = data[,data_colnames$distribution] - 1, # -1 to convert to CPP index
    Edims_ez = Edims_ez,
    S_kk = S_kk,
    Sdims = Sdims,
    Aistc_zz = Aistc_zz[which_Arows,,drop=FALSE],     # Index form, i, s, t
    Axi_z = Atriplet$x[which_Arows],
    RAM = as.matrix(na.omit(ram[,1:4])),
    RAMstart = as.numeric(ram[,5]),
    X_gj = matrix(0,ncol=ncol(X_ij),nrow=0),
    Z_gk = matrix(0,ncol=ncol(Z_ik),nrow=0),
    Agstc_zz = matrix(0,nrow=0,ncol=4),
    Axg_z = numeric(0),
    t_g = integer(0),
    c_g = integer(0),
    e_g = integer(0),
    W_gz = matrix(0,nrow=0,ncol=2),
    V_gz = matrix(0,nrow=0,ncol=2)
  )
  if( spatial_method_code %in% c(1,3,4) ){
    tmb_data$spatial_list = spatial_list
  }else if( spatial_method_code %in% 2 ){
    tmb_data$Adj = Adj
  }

  # make params
  tmb_par = list(
    log_kappa = log(1),
    alpha_j = rep(0,ncol(X_ij)),  # Spline coefficients
    gamma_k = rep(0,ncol(Z_ik)),  # Spline coefficients
    #omega = rep(0, spatial_graph$n),
    beta_z = ifelse(beta_type==1, 0.01, 1),
    log_lambda = rep(0,length(Sdims)), #Log spline penalization coefficients
    log_sigma = log_sigma,
    delta0_c = rep(0, length(variables)),
    #x_tc = matrix(0, nrow=tmb_data$n_t, ncol=tmb_data$n_c)
    epsilon_stc = array(0, dim=c(n_s, length(times), length(variables))),
    eps = numeric(0)
  )

  # Turn off initial conditions
  if( estimate_delta0==FALSE ){
    tmb_par$delta0_c = numeric(0)
  }

  # Turn of log_kappa when not needed
  tmb_map = list()
  if( tmb_data$spatial_method_code %in% c(3,4) ){
    tmb_map$log_kappa = factor(NA)
  }

  # User-supplied parameters
  if( !is.null(control$tmb_par) ){
    # Check shape but not numeric values, and give informative error
    attr(tmb_par,"check.passed") = attr(control$tmb_par,"check.passed")
    if( isTRUE(all.equal(control$tmb_par, tmb_par, tolerance=Inf)) ){
      tmb_par = control$tmb_par
    }else{
      stop("Not using `control$tmb_par` because it has some difference from `tmb_par` built internally")
    }
  }

  ##############
  # Fit model
  ##############

  if( FALSE ){
    setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
    dyn.unload(dynlib("tinyVAST"))
    compile("tinyVAST.cpp" , framework = "TMBad" )
    dyn.load(dynlib("tinyVAST"))
  }
  obj = MakeADFun( data = tmb_data,
                   parameters = tmb_par,
                   map = tmb_map,
                   random = c("gamma_k","epsilon_stc"),
                   DLL = "tinyVAST",
                   profile = control$profile )  #
  #openmp( ... , DLL="tinyVAST" )
  obj$env$beSilent()

  # Optimize
  opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isFALSE(control$quiet) ) message("Running nlminb_loop #", i)
    opt = nlminb( start = opt$par,
                  obj = obj$fn,
                  gr = obj$gr,
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isFALSE(control$quiet) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(opt$par) )
    h = optimHess(opt$par, fn=obj$fn, gr=obj$gr)
    opt$par = opt$par - solve(h, g)
    opt$objective = obj$fn(opt$par)
  }

  # Run sdreport
  if( isTRUE(control$getsd) ){
    if( isFALSE(control$quiet) ) message("Running sdreport")
    Hess_fixed = optimHess( par=opt$par, fn=obj$fn, gr=obj$gr )
    sdrep = sdreport( obj, hessian.fixed=Hess_fixed )
  }else{
    Hess_fixed = sdrep = NULL
  }

  # bundle and return output
  internal = list(
    ram_output = ram_output,
    sem = sem,
    data_colnames = data_colnames,
    times = times,
    variables = variables,
    parlist = obj$env$parList(par=obj$env$last.par.best),
    Hess_fixed = Hess_fixed,
    control = control
  )
  out = structure( list(
    formula = formula,
    data = data,
    gam_setup = gam_setup,
    obj = obj,
    opt = opt,
    rep = obj$report(obj$env$last.par.best),
    sdrep = sdrep,
    tmb_inputs = list(tmb_data=tmb_data, tmb_par=tmb_par, tmb_map=tmb_map),
    call = match.call(),
    spatial_graph = spatial_graph,
    data_colnames = data_colnames,
    run_time = Sys.time() - start_time,
    internal = internal
  ), class="tinyVAST" )
  return(out)
}

#' @title Control parameters for tinyVAST
#'
#' @inheritParams stats::nlminb
#' @inheritParams TMB::MakeADFun
#'
#' @param getsd Boolean indicating whether to call \code{\link[TMB]{sdreport}}
#' @param newton_loops Integer number of newton steps to do after running \code{nlminb}
#' @param tmb_par list of parameters for starting values, with shape identical to
#'        \code{fit(...)$internal$parlist}
#'
#' @export
tinyVASTcontrol <-
function( nlminb_loops = 1,
          newton_loops = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          quiet = FALSE,
          trace = 1,
          profile = c(),
          tmb_par = NULL ){

  # Return
  structure( list(
    nlminb_loops = nlminb_loops,
    newton_loops = newton_loops,
    eval.max = eval.max,
    iter.max = iter.max,
    getsd = getsd,
    quiet = quiet,
    trace = trace,
    profile = profile,
    tmb_par = tmb_par
  ), class = "tinyVASTcontrol" )
}

#' @title Print fitted tinyVAST object
#'
#' @description Prints output from fitted tinyVAST model
#'
#' @method print tinyVAST
#' @export
print.tinyVAST <-
function( x,
          ... ){
  print(x[c('call','opt','sdrep','run_time')])
}

#' Calculate residuals
#'
#' @title Calculate deviance or response residuals for tinyVAST
#'
#' @param object Output from \code{\link{fit}}
#' @param type which type of residuals to compute (only option is \code{"deviance"} or \code{"response"} for now)
#' @param ... Note used
#'
#' @method residuals tinyVAST
#' @export
residuals.tinyVAST <-
function( object,
          type = c("deviance","response"),
          ... ){

  # https://stats.stackexchange.com/questions/1432/what-do-the-residuals-in-a-logistic-regression-mean
  # Normal deviance residuals
  if( FALSE ){
    x = rnorm(10)
    y = x + rnorm(10)
    Glm = glm( y ~ 1 + x, family="gaussian")
    mu = predict(Glm,type="response")
    r1 = y - mu
    r1 - resid(Glm)
  }
  # Poisson deviance residuals
  if( FALSE ){
    x = rnorm(10)
    y = rpois(10, exp(x))
    Glm = glm( y ~ 1 + x, family="poisson")
    mu = predict(Glm,type="response")
    # https://stats.stackexchange.com/questions/398098/formula-for-deviance-residuals-for-poisson-model-with-identity-link-function
    r1 = sign(y - mu) * sqrt(2*(y*log((y+1e-10)/mu) - (y-mu)))
    r1 - resid(Glm)
  }
  # Binomial deviance residuals
  if( FALSE ){
    p = 0.5
    y = rbinom(10, prob=p, size=1)
    Glm = glm( y ~ 1, family="binomial")
    mu = predict(Glm, type="response")
    r1 = sign(y - mu) * sqrt(-2*(((1-y)*log(1-mu) + y*log(mu))))
    r1 - resid(Glm)
  }
  # Gamma deviance residuals
  if( FALSE ){
    mu = 1
    cv = 0.8
    y = rgamma( n=10, shape=1/cv^2, scale=mu*cv^2 )
    Glm = glm( y ~ 1, family=Gamma(link='log'))
    mu = predict(Glm, type="response")
    r1 = sign(y - mu) * sqrt(2 * ( (y-mu)/mu - log(y/mu) ))
    r1 - resid(Glm)
  }

  # Poisson: sign(y - mu) * sqrt(2*(ifelse(y==0, 0, y*log(y/mu)) - (y-mu)))
  # Binomial:  -2 * ((1-y)*log(1-mu) + y*log(mu))
  # Gamma: 2 * ( (y-mu)/mu - log(y/mu) )

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

# Extract the (marginal) log-likelihood of a tinyVAST model
#
# @return object of class \code{logLik} with attributes
#   \item{val}{log-likelihood}
#   \item{df}{number of parameters}
#' @importFrom stats logLik
#' @export
logLik.tinyVAST <- function(object, ...) {
  val = -1 * object$opt$objective
  df = length( object$opt$par )
  out = structure( val,
             df = df,
             class = "logLik")
  return(out)
}

#' Extract Variance-Covariance Matrix
#'
#' extract the covariance of fixed effects, or both fixed and random effects.
#'
#' @param object output from \code{fit}
#' @param which whether to extract the covariance among fixed effects, random effects, or both
#' @param ... ignored, for method compatibility
#' @importFrom stats vcov
#' @method vcov tinyVAST
#' @export
vcov.tinyVAST <-
function( object,
          which = c("fixed", "random", "both"),
          ...) {

  which = match.arg(which)

  if( which=="fixed" ){
    V = object$sdrep$cov.fixed
    if(is.null(V)){
      warning("Please re-run `dsem` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$sdrep$jointPrecision
    if(is.null(H)){
      warning("Please re-run `dsem` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

