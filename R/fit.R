#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set, and widely used interface for objects
#'
#' @useDynLib tinyVAST, .registration = TRUE
#' @export
fit <-
function( data,
          formula ){

  # Initial constructor of splines
  gam_setup = gam( formula, data = data, fit=FALSE ) # select doesn't do anything in this setup
  y_i = model.response(gam_setup$mf)  # model.extract(gam_setup$mf, "response")    # formula.tools::lhs

  #Extrtact penelization matrices
  S_list = lapply( seq_along(gam_setup$smooth), \(x) gam_setup$smooth[[x]]$S[[1]] )
  S_combined = .bdiag(S_list)         # join S's in sparse matrix
  Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
  if(is.null(Sdims)) Sdims = vector(length=0)

  #
  which_se = grep( pattern="s(", x=gam_setup$term.names, fixed=TRUE )
  X = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
  Z = gam_setup$X[,which_se,drop=FALSE]

  # make dat
  tmb_data = list(
    Y = y_i,
    X = X,
    Z = Z,
    S = S_combined,
    Sdims = Sdims,
    Xpred = matrix(ncol=ncol(X),nrow=0),
    Zpred = matrix(ncol=ncol(Z),nrow=0)
  )

  # make params
  tmb_par = list(
    beta = rep(0,ncol(X)),  # Spline coefficients
    gamma = rep(0,ncol(Z)),  # Spline coefficients
    log_lambda = rep(0,length(Sdims)), #Log spline penalization coefficients
    log_sigma = 0
  )

  #Fit model
  obj = MakeADFun( data=tmb_data, parameters=tmb_par, random="gamma" )
  obj$env$beSilent()
  opt = nlminb(obj$par,obj$fn,obj$gr, control=list(trace=1))
  sdrep = sdreport(obj)

  out = structure( list(
    gam_setup = gam_setup,
    obj = obj,
    opt = opt,
    rep = obj$report(obj$env$last.par.best),
    sdrep = sdrep,
    tmb_inputs = list(tmb_data=tmb_data, tmb_par=tmb_par),
    call = match.call()
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
