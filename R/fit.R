#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set, and widely used interface for objects
#'
#' @inheritParams dsem::make_dsem_ram
#'
#' @param sem Specification for structural equation model structure for
#'        constructing a space-variable interaction.
#'        See code{\link{make_sem_ram}} for more description.
#' @param dsem Specification for time-series structural equation model structure
#'        including lagged or simultaneous effects for
#'        constructing a space-variable interaction.  See
#'        \code{\link{make_dsem_ram}}  or \code{\link{make_eof_ram}}
#'        for more description
#' @param data Data-frame of predictor and response variables.
#' @param formula Formula with response on left-hand-side and predictors on right-hand-side,
#'        parsed by `mgcv` and hence allowing `s(.)` for splines or `offset(.)` for
#'        an offset.
#' @param family_link Matrix of length-two, indicating the distribution and link-function
#'        for each level of \code{data$dist}
#' @param spatial_graph Object that represents spatial relationships, either using `fmesher`
#'        to apply the SPDE method, `igraph` to apply a simultaneous autoregressive (SAR)
#'        process, or `NULL` to specify a single site.
#' @param control Output from [tinyVASTcontrol()], used to define user
#'        settings, and see documentation for that function for details.
#'
#' @details
#' `tinyVAST` includes four basic inputs that specify the model structure:
#' * `formula` specifies covariates and splines in a Generalized Additive Model;
#' * `dsem` specifies interactions among variables and over time, constructing
#'   the space-time-variable interaction.
#' * `sem` specifies interactions among variables and over time, constructing
#'   the space-variable interaction.
#' * `spatial_graph` specifies spatial correlations
#'
#' the default `dsem=NULL` turns off all multivariate and temporal indexing, such
#' that `spatial_graph` is then ignored, and the model collapses
#' to a standard model using `mgcv::gam`.  To specify a univeriate spatial model,
#' the user must specify both `spatial_graph` and `dsem=""`, where the latter
#' is then parsed to include a single exogenous variance for the single variable
#'
#' | \strong{Model type} | \strong{How to specify} |
#' | --- | --- |
#' | Generalized additive model | specify `spatial_graph=NULL` and `dsem=""`, and then use `formula` to specify splines and covariates |
#' | Dynamic structural equation model (including vector autoregressive, dynamic factor analysis, ARIMA, and structural equation models) | specify `spatial_graph=NULL` and use `dsem` to specify interactions among variables and over time |
#' | Univeriate spatial model | specify `spatial_graph` and `dsem=""`, where the latter is then parsed to include a single exogenous variance for the single variable |
#' | Multivariate spatial model | specify `spatial_graph` and use `dsem` (without any lagged effects) to specify spatial interactions |
#' | Vector autoregressive spatio-temporal model | specify `spatial_graph` and use `dsem=""` to specify interactions among variables and over time, where spatio-temporal variables are constructed via the separable interaction of `dsem` and `spatial_graph` |
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom sem specifyModel specifyEquations
#' @importFrom corpcor pseudoinverse
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
          dsem = NULL,
          family_link = rbind( "obs"=c(0,0) ),
          data_colnames = list("spatial"=c("x","y"), "variable"="var", "time"="time", "distribution"="dist"),
          times = seq(min(data[,data_colnames$time]),max(data[,data_colnames$time])),
          variables = unique(data[,data_colnames$variable]),
          spatial_graph = NULL,
          control = tinyVASTcontrol(),
          ... ){

  # https://roxygen2.r-lib.org/articles/rd-formatting.html#tables for roxygen formatting
  start_time = Sys.time()

  # General error checks
  if( class(control) != "tinyVASTcontrol" ) stop("`control` must be made by `tinyVASTcontrol()`")
  if( !is.data.frame(data) ) stop("`data` must be a data frame")

  ##############
  # SEM / DSEM telescoping
  ##############

  # Allow user to avoid specifying these if is.null(dsem)
  if( is.null(dsem) ){
    times = numeric(0)
    if( !(data_colnames$time %in% colnames(data)) ){
      data = data.frame( data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$time)) )
    }
  }

  # Allow user to avoid specifying these if is.null(dsem) AND is.null(sem)
  if( is.null(dsem) & is.null(sem) ){
    variables = numeric(0)
    if( !(data_colnames$variable %in% colnames(data)) ){
      data = data.frame(data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$variable)))
    }
  }

  ##############
  # DSEM RAM constructor
  ##############

  # (I-Rho)^-1 * Gamma * (I-Rho)^-1
  if( is.null(dsem) ){
    dsem_ram_output = list(
      ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
      model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
    )
  }else if( isTRUE(is.character(dsem)) ){
    dsem_ram_output = make_dsem_ram( dsem, times=times, variables=variables, quiet=control$quiet, covs=variables )
  }else if( class(dsem) %in% c("dsem_ram","eof_ram") ){
    dsem_ram_output = dsem
  }else{
    stop("`dsem` must be either `NULL` or a character-string")
  }
  ram_dsem = dsem_ram_output$ram

  # Identify arrow-type for each beta_j estimated in RAM
  which_nonzero = which(ram_dsem[,4]>0)
  beta_type = tapply( ram_dsem[which_nonzero,1], INDEX=ram_dsem[which_nonzero,4], max)

  # Error checks
  if( class(dsem_ram_output) %in% c("dsem_ram") ){
    if( any((dsem_ram_output$model[,'direction']==2) & (dsem_ram_output$model[,2]!=0)) ){
      stop("All two-headed arrows should have lag=0")
    }
    if( !all(c(dsem_ram_output$model[,'first'],dsem_ram_output$model[,'second']) %in% variables) ){
      stop("Some variable in `dsem` is not in `tsdata`")
    }
  }

  # Check for rank-deficient precision from RAM
  ram2 = subset( data.frame(dsem_ram_output$ram), heads==2 )
  total_variance_h = tapply( as.numeric(ifelse(is.na(ram2$start), 1, ram2$start)),
          INDEX = ram2$from, FUN=\(x)sum(abs(x)) )
  ram1 = subset( data.frame(dsem_ram_output$ram), heads==1 )
  total_effect_h = tapply( as.numeric(ifelse(is.na(ram1$start), 1, ram1$start)),
          INDEX = ram1$from, FUN=\(x)sum(abs(x)) )
  if( any(total_variance_h==0) & control$gmrf_parameterization=="separable" ){
    stop("Must use gmrf_parameterization=`projection` for the dsem RAM supplied")
  }

  ##############
  # SEM RAM constructor
  ##############

  if( is.null(sem) ){
    sem_ram_output = list(
      ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
      model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
    )
  }else if( isTRUE(is.character(sem)) ){
    sem_ram_output = make_sem_ram( sem, variables=as.character(variables), quiet=control$quiet, covs=as.character(variables) )
  } else {
    stop("`sem` must be either `NULL` or a character-string")
  }
  ram_sem = sem_ram_output$ram

  # Identify arrow-type for each beta_j estimated in RAM
  which_nonzero = which(ram_sem[,4]>0)
  theta_type = tapply( ram_sem[which_nonzero,1], INDEX=ram_sem[which_nonzero,4], max)

  # Check for rank-deficient precision from RAM
  ram2 = subset( data.frame(sem_ram_output$ram), heads==2 )
  total_variance_h = tapply( as.numeric(ifelse( is.na(ram2$start), 1, ram2$start)),
          INDEX = ram2$from, FUN=\(x)sum(abs(x)) )
  ram1 = subset( data.frame(sem_ram_output$ram), heads==1 )
  total_effect_h = tapply( as.numeric(ifelse(is.na(ram1$start), 1, ram1$start)),
          INDEX = ram1$from, FUN=\(x)sum(abs(x)) )
  if( any(total_variance_h==0) & control$gmrf_parameterization=="separable" ){
    stop("Must use options$gmrf_parameterization=`projection` for the sem RAM supplied")
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
  }else{      # if( !is.null(sem) )
    # Single-site
    spatial_method_code = 3
    n_s = 1
    A_is = matrix(1, nrow=nrow(data), ncol=1)    # dgCMatrix
    A_is = as(Matrix(A_is),"CsparseMatrix")
    spatial_list = list( "M0" = as(Matrix(1,nrow=1,ncol=1),"CsparseMatrix"),
                         "M1" = as(Matrix(0,nrow=1,ncol=1),"CsparseMatrix"),
                         "M2" = as(Matrix(0,nrow=1,ncol=1),"CsparseMatrix") )
  }
  Atriplet = Matrix::mat2triplet(A_is)

  ##############
  # Formula constructor
  ##############

  # Initial constructor of splines
  gam_setup = gam( formula, data = data, fit=FALSE ) # select doesn't do anything in this setup
  y_i = model.response(gam_setup$mf)  # OR USE: model.extract(gam_setup$mf, "response")
  offset_i = gam_setup$offset

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

  # Telescope
  if( !(data_colnames$distribution %in% colnames(data)) ){
    if(nrow(family_link)>1) stop("Must supply `dist` if using multiple `family_link` options")
    data = data.frame( data, matrix(rownames(family_link)[1], nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$distribution)) )
  }

  # Build e_i ... always has length nrow(data)
  e_i = match( data[,data_colnames$distribution], rownames(family_link) )

  # Check for errors
  if( (any(is.na(e_i))) ){
    stop("`data[,data_colnames$distribution]` has values that don't match `rownames(family_link)`")
  }
  if( any( (family_link[,2]==0) & (family_link[,1] %in% c(1,2)) ) ){
    warning("Using an identify link with a Tweedie or lognormal distribution is not advised.")
  }

  # Telescope format
  if( !is.matrix(family_link) ){
    family_link = matrix( family_link, nrow=1, ncol=2 )
  }

  # Construct log_sigma based on family_link
  remove_last = \(x) x[-length(x)]
  Nsigma_e = sapply( as.character(family_link[,1]), "0"=1, "1"=2, "2"=1, FUN=switch )
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

  #
  Aepsilon_zz = cbind(Atriplet$i, Atriplet$j, t_i[Atriplet$i], c_i[Atriplet$i])
  which_Arows = which(apply( Aepsilon_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) & any(x>0) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( nrow(ram_dsem)==0 ){
    which_Arows = numeric(0)
  }
  Aepsilon_zz = Aepsilon_zz[which_Arows,,drop=FALSE]
  Aepsilon_z = Atriplet$x[which_Arows]

  #
  Aomega_zz = cbind(Atriplet$i, Atriplet$j, c_i[Atriplet$i])
  which_Arows = which(apply( Aomega_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( nrow(ram_sem)==0 ){
    which_Arows = numeric(0)
  }
  Aomega_zz = Aomega_zz[which_Arows,,drop=FALSE]
  Aomega_z = Atriplet$x[which_Arows]

  # make dat
  tmb_data = list(
    spatial_options = c(spatial_method_code, ifelse(control$gmrf_parameterization=="separable",0,1) ),
    y_i = y_i,
    X_ij = X_ij,
    Z_ik = Z_ik,
    t_i = t_i - 1, # -1 to convert to CPP index
    c_i = c_i - 1, # -1 to convert to CPP index
    offset_i = offset_i,
    f_ez = family_link,
    e_i = e_i - 1, # -1 to convert to CPP index
    Edims_ez = Edims_ez,
    S_kk = S_kk,
    Sdims = Sdims,
    Aepsilon_zz = Aepsilon_zz - 1,     # Index form, i, s, t
    Aepsilon_z = Aepsilon_z,
    Aomega_zz = Aomega_zz - 1,     # Index form, i, s, t
    Aomega_z = Aomega_z,
    ram_sem = as.matrix(na.omit(ram_sem[,1:4])),
    ram_sem_start = as.numeric(ram_sem[,5]),
    ram_dsem = as.matrix(na.omit(ram_dsem[,1:4])),
    ram_dsem_start = as.numeric(ram_dsem[,5]),
    X_gj = matrix(0,ncol=ncol(X_ij),nrow=0),
    Z_gk = matrix(0,ncol=ncol(Z_ik),nrow=0),
    AepsilonG_zz = matrix(0,nrow=0,ncol=4),
    AepsilonG_z = numeric(0),
    AomegaG_zz = matrix(0,nrow=0,ncol=4),
    AomegaG_z = numeric(0),
    t_g = integer(0),
    c_g = integer(0),
    offset_g = integer(0),
    e_g = integer(0),
    W_gz = matrix(0,nrow=0,ncol=2),
    V_gz = matrix(0,nrow=0,ncol=2)
  )
  if( spatial_method_code %in% c(1,3) ){
    tmb_data$spatial_list = spatial_list
  }else if( spatial_method_code %in% 2 ){
    tmb_data$Adj = Adj
  }

  # make params
  tmb_par = list(
    alpha_j = rep(0,ncol(X_ij)),  # Spline coefficients
    gamma_k = rep(0,ncol(Z_ik)),  # Spline coefficients
    beta_z = as.numeric(ifelse(beta_type==1, 0.01, 1)),  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
    theta_z = as.numeric(ifelse(theta_type==1, 0.01, 1)),
    log_lambda = rep(0,length(Sdims)), #Log spline penalization coefficients
    log_sigma = log_sigma,
    delta0_c = rep(0, length(variables)),
    epsilon_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega_sc = array(0, dim=c(n_s, length(variables))),
    eps = numeric(0),
    log_kappa = log(1)
  )

  # Telescoping
  if( nrow(ram_dsem)==0 ){
    tmb_par$epsilon_stc = tmb_par$epsilon_stc[,numeric(0),,drop=FALSE]   # Keep c original length so n_c is detected correctly
  }
  if( nrow(ram_sem)==0 ){
    tmb_par$omega_sc = tmb_par$omega_sc[,numeric(0),drop=FALSE]
  }

  # Turn off initial conditions
  if( control$estimate_delta0==FALSE ){
    tmb_par$delta0_c = numeric(0)
  }

  # Turn of log_kappa when not needed
  if( spatial_method_code %in% c(3) ){
    tmb_par = tmb_par[-match("log_kappa",names(tmb_par))]
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

  # Check for obvious issues ... no NAs except in RAMstart
  if( any(is.na(tmb_data[-match(c("ram_sem_start","ram_dsem_start"),names(tmb_data))])) ){
    stop("Check `tmb_data` for NAs")
  }

  # Empty map, but leaving for future needs
  tmb_map = list()

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
                   random = c("gamma_k","epsilon_stc","omega_sc"),
                   DLL = "tinyVAST",
                   profile = control$profile )  #
  #openmp( ... , DLL="tinyVAST" )
  obj$env$beSilent()
  # L = rep$IminusRho_hh %*% rep$Gamma_hh

  # Optimize
  #start_time = Sys.time()
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
  #Sys.time() - start_time
  #opt

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
    dsem_ram_output = dsem_ram_output,
    sem_ram_output = sem_ram_output,
    dsem = dsem,
    sem = sem,
    data_colnames = data_colnames,
    times = times,
    variables = variables,
    parlist = obj$env$parList(par=obj$env$last.par.best),
    Hess_fixed = Hess_fixed,
    control = control,
    family_link = family_link
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
#' @param getsd Boolean indicating whether to call [TMB::sdreport()]
#' @param newton_loops Integer number of newton steps to do after running `nlminb`
#' @param tmb_par list of parameters for starting values, with shape identical to
#'        `fit(...)$internal$parlist`
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
          tmb_par = NULL,
          gmrf_parameterization = c("separable","projection"),
          estimate_delta0 = FALSE ){

  gmrf_parameterization = match.arg(gmrf_parameterization)

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
    tmb_par = tmb_par,
    gmrf_parameterization = gmrf_parameterization,
    estimate_delta0 = estimate_delta0
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


#' @title summarize tinyVAST
#'
#' @description summarize parameters from a fitted tinyVAST
#'
#' @details
#' tinyVAST includes an "arrow and lag" notation, which specifies the set of
#' path coefficients and exogenous variance parameters to be estimated. Function \code{fit}
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
#' @param object Output from \code{\link{fit}}
#' @param what What component to summarize
#' @param ... Not used
#'
#' @method summary tinyVAST
#' @export
summary.tinyVAST <-
function( object,
          what = c("sem"),
          ... ){

  # Easy of use
  what = match.arg(what)
  if( what=="sem" ){
    model = object$internal$sem_ram_output$ram
    model$to = as.character(object$internal$variables)[model$to]
    model$from = as.character(object$internal$variables)[model$from]
    ParHat = object$obj$env$parList()

    #
    coefs = data.frame( model, "Estimate"=c(NA,ParHat$theta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,5]), coefs$Estimate )
    if( "sdrep" %in% names(object) ){
      SE = as.list( object$sdrep, report=FALSE, what="Std. Error")
      coefs = data.frame( coefs, "Std_Error"=c(NA,SE$theta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
    #coefs
  }

  return(coefs)
}

#' Calculate residuals
#'
#' @title Calculate deviance or response residuals for tinyVAST
#'
#' @param object Output from [fit()]
#' @param type which type of residuals to compute (only option is `"deviance"` or `"response"` for now)
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
#' @param object output from `fit`
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
      warning("Please re-run `tinyVAS` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$sdrep$jointPrecision
    if(is.null(H)){
      warning("Please re-run `tinyVAS` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

