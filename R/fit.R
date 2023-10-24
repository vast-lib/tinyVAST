#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set, and widely used interface for objects
#'
#' @inheritParams dsem::dsem
#'
#' @param data Data-frame of predictor and response variables.
#' @param formula Formula with response on left-hand-side and predictors on right-hand-side,
#'          parsed by \code{mgcv} and hence allowing \code{s(.)} for splines
#' @param family_link Vector of length-two, indicating the distribution and link-function
#' @param spatial_graph Object that represents spatial relationships, either using \code{fmesher}
#'        to apply the SPDE method, or \code{igraph} to apply a simultaneous autoregressive (SAR)
#'        process.
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
          data_colnames = list("spatial"=c("x","y"), "variable"="var", "time"="time"),
          times = sort(unique(data[,data_colnames$time])),
          variables = unique(data[,data_colnames$variable]),
          spatial_graph = NULL,
          quiet = FALSE,
          ... ){

  start_time = Sys.time()
  # SEM stuff
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
      data = cbind(data, "time"=1)
    }
    if( !(data_colnames$variable %in% colnames(data)) ){
      data = cbind(data, "var"=1)
    }
  }else{
    ram_output = make_ram( sem, times=times, variables=variables, quiet=quiet, covs=variables )
  }
  ram = ram_output$ram
  # Error checks
  if( any((ram_output$model[,'direction']==2) & (ram_output$model[,2]!=0)) ){
    stop("All two-headed arrows should have lag=0")
  }
  if( !all(c(ram_output$model[,'first'],ram_output$model[,'second']) %in% variables) ){
    stop("Some variable in `sem` is not in `tsdata`")
  }

  # Spatial domain constructor
  if( "fm_mesh_2d" %in% class(spatial_graph) ){
    n_s = spatial_graph$n
    spatial_method_code = 1
    spatial_list = fm_fem( spatial_graph )
    spatial_list = list("M0"=spatial_list$c0, "M1"=spatial_list$g1, "M2"=spatial_list$g2)
    A_is = fm_evaluator( spatial_graph, loc=as.matrix(data[,data_colnames$spatial]) )$proj$A
  } else if( "igraph" %in% class(spatial_graph) ) {
    spatial_method_code = 2
    Adj = as_adjacency_matrix( spatial_graph, sparse=TRUE )
    n_s = nrow(Adj)
    Match = match( data[,data_colnames$spatial], rownames(Adj) )
    if(any(is.na(Match))) stop("Check `spatial_graph` for SAR")
    A_is = sparseMatrix( i=1:nrow(data), j=Match, x=rep(1,nrow(data)) )
  }else if( !is.null(sem) ){
    spatial_method_code = 3
    n_s = 1
    A_is = matrix(1, nrow=nrow(data), ncol=1)    # dgCMatrix
    A_is = as(Matrix(A_is),"dgCMatrix")
    spatial_list = list( "M0" = as(Matrix(1,nrow=1,ncol=1),"dgCMatrix"),
                         "M1" = as(Matrix(0,nrow=1,ncol=1),"dgCMatrix"),
                         "M2" = as(Matrix(0,nrow=1,ncol=1),"dgCMatrix") )
  } else {
    spatial_method_code = 4
    n_s = 0
    A_is = Matrix(nrow=nrow(data), ncol=0)    # dgCMatrix
    A_is = as(A_is,"dgCMatrix")
    spatial_list = list( "M0" = as(Matrix(nrow=0,ncol=0),"dgCMatrix"),
                         "M1" = as(Matrix(nrow=0,ncol=0),"dgCMatrix"),
                         "M2" = as(Matrix(nrow=0,ncol=0),"dgCMatrix") )
  }
  Atriplet = Matrix::mat2triplet(A_is)

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

  # Turn of t_i and c_i when times and variables are missing, so that delta_k isn't built
  if( length(times>0) ){
    t_i = match( data[,data_colnames$time], times )
  }else{ t_i = integer(0) }
  if( length(variables>0) ){
    c_i = match( data[,data_colnames$var], variables )
  }else{ c_i = integer(0) }

  # make dat
  tmb_data = list(
    spatial_method_code = spatial_method_code,
    y_i = y_i,
    X_ij = X_ij,
    Z_ik = Z_ik,
    t_i = t_i - 1, # -1 to convert to CPP index
    c_i = c_i - 1, # -1 to convert to CPP index
    f_z = family_link,
    S_kk = S_kk,
    Sdims = Sdims,
    Aistc_zz = cbind(Atriplet$i, Atriplet$j, t_i[Atriplet$i], c_i[Atriplet$i]) - 1,     # Index form, i, s, t
    Axi_z = Atriplet$x,
    RAM = as.matrix(na.omit(ram[,1:4])),
    RAMstart = as.numeric(ram[,5]),
    X_gj = matrix(0,ncol=ncol(X_ij),nrow=0),
    Z_gk = matrix(0,ncol=ncol(Z_ik),nrow=0),
    Agstc_zz = matrix(0,nrow=0,ncol=4),
    Axg_z = numeric(0),
    t_g = integer(0),
    c_g = integer(0),
    W_gz = matrix(0,nrow=0,ncol=2),
    V_gz = matrix(0,nrow=0,ncol=2)
  )
  if( spatial_method_code %in% c(1,3,4) ){
    tmb_data$spatial_list = spatial_list
  }else if( spatial_method_code %in% 2 ){
    tmb_data$Adj = Adj
  }

  # Scale starting values with higher value for two-headed than one-headed arrows
  which_nonzero = which(ram[,4]>0)
  beta_type = tapply( ram[which_nonzero,1], INDEX=ram[which_nonzero,4], max)

  #
  if( family_link[1]==0 ){
    log_sigma = c(0)
  }else if( family_link[1]==1 ){
    log_sigma = c(0, 0)
  }else{ stop("Check `family_link`")}

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
    epsilon_stc = array(0, dim=c(n_s, length(times), length(variables)) )
  )
  # Turn off initial conditions
  if( estimate_delta0==FALSE ){
    tmb_par$delta0_c = numeric(0)
  }

  #
  tmb_map = list()
  if( tmb_data$spatial_method_code %in% c(3,4) ){
    tmb_map$log_kappa = factor(NA)
    #tmb_map$log_tau = factor(NA)
  }

  #Fit model
  if( FALSE ){
    setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
    dyn.unload(dynlib("tinyVAST"))
    compile("tinyVAST.cpp")
    dyn.load(dynlib("tinyVAST"))
  }
  obj = MakeADFun( data=tmb_data, parameters=tmb_par, map=tmb_map, random=c("gamma_k","epsilon_stc"), DLL="tinyVAST" )  #
  openmp( ... , DLL="tinyVAST" )
  # Experiment with Q_jonit
  if( FALSE ){
    rep = obj$report()
    error = list( Q_spatial=rep$Q_spatial, Q_kk=rep$Q_kk, epsilon_sk=rep$epsilon_sk )
    saveRDS( error, "error.RDS" )
  }
  if( FALSE ){
    rep = obj$report()
    Matrix::image(rep$Q_joint)
  }
  if( FALSE ){
    obj = MakeADFun( data=tmb_data, parameters=tmb_par, map=tmb_map, DLL="tinyVAST" )
  }
  if( FALSE ){
    H = obj$env$spHess(random=TRUE)
    Matrix::image(H)
  }
  obj$env$beSilent()
  opt = nlminb( start=obj$par, obj=obj$fn, gr=obj$gr,
                control=list(eval.max=1e4, iter.max=1e4, trace=ifelse(isTRUE(quiet),0,1)) )
  Hess_fixed = optimHess( par=opt$par, fn=obj$fn, gr=obj$gr )
  sdrep = sdreport( obj, hessian.fixed=Hess_fixed )

  # bundle and return output
  internal = list(
    ram_output = ram_output,
    sem = sem,
    data_colnames = data_colnames,
    times = times,
    variables = variables,
    parlist = obj$env$parList(par=obj$env$last.par.best),
    Hess_fixed = Hess_fixed
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
  ), class="tinyVAST"
  )
  return(out)
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

