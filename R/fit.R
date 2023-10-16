#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set, and widely used interface for objects
#'
#' @useDynLib tinyVAST, .registration = TRUE
#' @export
fit <-
function( data,
          formula,
          #spatial_sem,
          #spatiotemporal_sem,
          #times,
          #variables,
          data_colnames = list("spatial"=c("x","y"), "variables"="var", "time"="time"),
          spatial_graph = NULL,
          quiet = FALSE ){

  # Initial constructor of splines
  gam_setup = gam( formula, data = data, fit=FALSE ) # select doesn't do anything in this setup
  y_i = model.response(gam_setup$mf)  # model.extract(gam_setup$mf, "response")    # formula.tools::lhs

  # Extrtact and combine penelization matrices
  S_list = lapply( seq_along(gam_setup$smooth), \(x) gam_setup$smooth[[x]]$S[[1]] )
  S_combined = .bdiag(S_list)         # join S's in sparse matrix
  Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
  if(is.null(Sdims)) Sdims = vector(length=0)

  # Get covariates
  which_se = grep( pattern="s(", x=gam_setup$term.names, fixed=TRUE )
  X = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
  Z = gam_setup$X[,which_se,drop=FALSE]

  #
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
    Y = y_i,
    X = X,
    Z = Z,
    S = S_combined,
    Sdims = Sdims,
    spatial_list = spatial_list,
    A_is = A_is,
    Xpred = matrix(ncol=ncol(X),nrow=0),
    Zpred = matrix(ncol=ncol(Z),nrow=0),
    Apred_is = matrix(nrow=0,ncol=ncol(A_is))
  )

  # make params
  tmb_par = list(
    log_kappa = log(1),
    log_tau = log(1),
    beta = rep(0,ncol(X)),  # Spline coefficients
    gamma = rep(0,ncol(Z)),  # Spline coefficients
    omega = rep(0, spatial_graph$n),
    log_lambda = rep(0,length(Sdims)), #Log spline penalization coefficients
    log_sigma = 0
  )

  #
  tmb_map = list()
  if( spatial_graph$n==0 ){
    tmb_map$log_kappa = factor(NA)
    tmb_map$log_tau = factor(NA)
  }

  #Fit model
  #if( FALSE ){
  # setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
  # dyn.unload(dynlib("tinyVAST2"))
  # compile("tinyVAST2.cpp")
  # dyn.load(dynlib("tinyVAST2"))
  #}
  obj = MakeADFun( data=tmb_data, parameters=tmb_par, map=tmb_map, random=c("gamma","omega") )
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
