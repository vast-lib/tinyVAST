#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set, and widely used interface for objects
#'
#' @importFrom dsem make_ram classify_variables parse_path
#'
#' @useDynLib tinyVAST, .registration = TRUE
#' @export
fit <-
function( data,
          formula,
          sem = NULL,
          estimate_delta0 = FALSE,
          #spatial_sem,
          #spatiotemporal_sem,
          times,
          variables,
          data_colnames = list("spatial"=c("x","y"), "variable"="var", "time"="time"),
          spatial_graph = NULL,
          quiet = FALSE ){

  # SEM stuff
  # (I-Rho)^-1 * Gamma * (I-Rho)^-1
  if( is.null(sem) ){
    out = list(
      ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
      model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
    )
    times = numeric(0)
    variables = numeric(0)
    # Allow user to avoid specifying these if is.null(sem)
    if( !(data_colnames$time %in% colnames(data)) ){
      data = cbind(data, 1)
      colnames(data)[ncol(data)] = data_colnames$time
    }
    if( !(data_colnames$variable %in% colnames(data)) ){
      data = cbind(data, 1)
      colnames(data)[ncol(data)] = data_colnames$variable
    }
  }else{
    out = make_ram( sem, times=times, variables=variables, quiet=quiet, covs=variables )
  }
  ram = out$ram
  # Error checks
  if( any((out$model[,'direction']==2) & (out$model[,2]!=0)) ){
    stop("All two-headed arrows should have lag=0")
  }
  if( !all(c(out$model[,'first'],out$model[,'second']) %in% variables) ){
    stop("Some variable in `sem` is not in `tsdata`")
  }

  # Initial constructor of splines
  gam_setup = gam( formula, data = data, fit=FALSE ) # select doesn't do anything in this setup
  y_i = model.response(gam_setup$mf)  # OR USE: model.extract(gam_setup$mf, "response")

  # Extrtact and combine penelization matrices
  S_list = lapply( seq_along(gam_setup$smooth), \(x) gam_setup$smooth[[x]]$S[[1]] )
  S_combined = .bdiag(S_list)         # join S's in sparse matrix
  Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
  if(is.null(Sdims)) Sdims = vector(length=0)

  # Get covariates
  which_se = grep( pattern="s(", x=gam_setup$term.names, fixed=TRUE )
  X = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
  Z = gam_setup$X[,which_se,drop=FALSE]

  # Spatial domain constructor
  if( "fm_mesh_2d" %in% class(spatial_graph) ){
    spatial_list = fm_fem( spatial_graph )
    spatial_list = list("M0"=spatial_list$c0, "M1"=spatial_list$g1, "M2"=spatial_list$g2)
    A_is = fm_evaluator( spatial_graph, loc=as.matrix(data[,data_colnames$spatial]) )$proj$A
  } else {
    spatial_graph = list( "n"=0 )
    A_is = Matrix(nrow=nrow(data), ncol=0)    # dgCMatrix
    A_is = as(A_is,"dgCMatrix")
    spatial_list = list( "M0" = as(Matrix(nrow=0,ncol=0),"dgCMatrix"),
                         "M1" = as(Matrix(nrow=0,ncol=0),"dgCMatrix"),
                         "M2" = as(Matrix(nrow=0,ncol=0),"dgCMatrix") )
  }

  # make dat
  tmb_data = list(
    n_t = length(times),
    n_c = length(variables),
    Y = y_i,
    X = X,
    Z = Z,
    t_i = match( data[,data_colnames$time], times ) - 1,  # -1 to convert to CPP index
    c_i = match( data[,data_colnames$var], variables ) - 1, # -1 to convert to CPP index
    S = S_combined,
    Sdims = Sdims,
    spatial_list = spatial_list,
    A_is = A_is,
    RAM = as.matrix(na.omit(ram[,1:4])),
    RAMstart = as.numeric(ram[,5]),
    Xpred = matrix(ncol=ncol(X),nrow=0),
    Zpred = matrix(ncol=ncol(Z),nrow=0),
    Apred_is = matrix(nrow=0,ncol=ncol(A_is))
  )

  # Scale starting values with higher value for two-headed than one-headed arrows
  which_nonzero = which(ram[,4]>0)
  beta_type = tapply( ram[which_nonzero,1], INDEX=ram[which_nonzero,4], max)

  # make params
  tmb_par = list(
    log_kappa = log(1),
    log_tau = log(1),
    alpha = rep(0,ncol(X)),  # Spline coefficients
    gamma = rep(0,ncol(Z)),  # Spline coefficients
    omega = rep(0, spatial_graph$n),
    beta_z = ifelse(beta_type==1, 0.01, 1),
    log_lambda = rep(0,length(Sdims)), #Log spline penalization coefficients
    log_sigma = 0,
    delta0_c = rep(0, tmb_data$n_c),
    x_tc = matrix(0, nrow=tmb_data$n_t, ncol=tmb_data$n_c)
  )
  # Turn off initial conditions
  if( estimate_delta0==FALSE ){
    tmb_par$delta0_c = numeric(0)
  }

  #
  tmb_map = list()
  if( spatial_graph$n==0 ){
    tmb_map$log_kappa = factor(NA)
    tmb_map$log_tau = factor(NA)
  }

  #Fit model
  if( FALSE ){
    setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
    dyn.unload(dynlib("tinyVAST"))
    compile("tinyVAST.cpp")
    dyn.load(dynlib("tinyVAST"))
  }
  obj = MakeADFun( data=tmb_data, parameters=tmb_par, map=tmb_map, random=c("gamma","omega","x_tc"), DLL="tinyVAST" )
  obj$env$beSilent()
  opt = nlminb( start=obj$par, obj=obj$fn, gr=obj$gr,
                control=list(eval.max=1e4, iter.max=1e4, trace=ifelse(isTRUE(quiet),0,1)) )
  sdrep = sdreport(obj)

  # bundle and return output
  out = structure( list(
    formula = formula,
    data = data,
    gam_setup = gam_setup,
    obj = obj,
    opt = opt,
    rep = obj$report(obj$env$last.par.best),
    sdrep = sdrep,
    tmb_inputs = list(tmb_data=tmb_data, tmb_par=tmb_par),
    call = match.call(),
    spatial_graph = spatial_graph,
    data_colnames = data_colnames
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
  print(x$opt)
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
  print(x$opt)
}
