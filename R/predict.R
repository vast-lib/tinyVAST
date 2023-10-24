#' @title Predict using vector autoregressive spatio-temporal model
#'
#' @description Predicts values given new covariates using a tinyVAST model
#'
#' @inheritParams add_predictions
#'
#' @param ... Not used.
#'
#' @method predict tinyVAST
#' @export
predict.tinyVAST <-
function( object,
          newdata,
          ... ){

  # extract original X and Z
  if(missing(newdata)) newdata = object$data

  # Build new
  tmb_data2 = add_predictions( object=object, newdata=newdata )

  # Area-expanded sum
  #tmb_data2$W_gz = matrix(1, nrow=nrow(newdata), ncol=2)
  #tmb_data2$V_gz = matrix(0, nrow=nrow(newdata), ncol=2)

  # Abundance-weighted z
  #tmb_data2$W_gz = cbind( 1, newdata$x )
  #tmb_data2$V_gz = matrix(1, nrow=nrow(newdata), ncol=2)

  # Rebuild object
  newobj = MakeADFun( data = tmb_data2,
                      parameters = object$internal$parlist,
                      map = object$tmb_inputs$tmb_map,
                      random = c("gamma_k","epsilon_stc"),
                      DLL = "tinyVAST" )
  out = newobj$report()$mu_g

  return(out)
}

#' @title Monte Carlo integration for abundance
#'
#' @description Applies Monte Carlo integration to approximate area-expanded abundance
#'
#' @inheritParams TMB::sdreport
#' @inheritParams add_predictions
#'
#' @param area value used for area-weighted expansion of estimated density surface
#'     for each row of \code{newdata}.
#'
#' @export
integrate_output <-
function( object,
          newdata,
          area,
          bias.correct = TRUE ){

  # extract original X and Z
  if(missing(newdata)) newdata = object$data
  if(missing(area)) area = rep(1, nrow(newdata))
  if(length(area)==1) area = rep(area, nrow(newdata))

  # Build new
  tmb_data2 = add_predictions( object=object, newdata=newdata )

  # Area-expanded sum
  tmb_data2$W_gz = matrix(area, nrow=nrow(newdata), ncol=2)
  tmb_data2$V_gz = matrix(0, nrow=nrow(newdata), ncol=2)

  # Abundance-weighted z
  #tmb_data2$W_gz = cbind( 1, newdata$x )
  #tmb_data2$V_gz = matrix(1, nrow=nrow(newdata), ncol=2)

  # Rebuild object
  newobj = MakeADFun( data = tmb_data2,
                      parameters = object$internal$parlist,
                      map = object$tmb_inputs$tmb_map,
                      random = c("gamma_k","epsilon_stc"),
                      DLL = "tinyVAST" )
  newobj$env$beSilent()

  # Run sdreport
  newsd = sdreport( obj = newobj,
                    par.fixed = object$opt$par,
                    hessian.fixed = object$internal$Hess_fixed,
                    bias.correct = bias.correct )
  out = summary(newsd, "report")['Metric',]

  return(out)
}

#' @title Add predictions to data-list
#'
#' @description Given user-provided `newdata`, expand the object `tmb_data`
#'    to include predictions corresponding to those new observations
#'
#' @param object Output from [fit()].
#' @param newdata New data-frame of independent variables used to predict the response.
#'
#' @export
add_predictions <-
function( object,
          newdata ){

  tmb_data = object$tmb_inputs$tmb_data
  gam = object$gam_setup
  origdata = object$data

  # Check newdata for missing variables and/or factors
  pred_set = all.vars( object$gam_setup$pred.formula )
  for(pred in pred_set ){
    if( !(pred %in% colnames(newdata)) ){
      stop("Missing ", pred, " in newdata")
    }
    if( is.factor(origdata[,pred]) ){
      newdata[,pred] = factor(newdata[,pred], levels=levels(origdata[,pred]))
    }
  }

  # Assemble predZ
  Z_gk = lapply( seq_along(gam$smooth),
                  FUN = \(x) PredictMat(gam$smooth[[x]], data=newdata) )
  Z_gk = do.call( cbind, Z_gk )
  if(is.null(Z_gk)) Z_gk = matrix(0, nrow=nrow(newdata),ncol=ncol(tmb_data$Z_ik))

  # Assemble predX
  formula_no_sm = remove_s_and_t2(gam$formula)
  mf1 = model.frame(formula_no_sm, origdata)
  terms1 = attr(mf1, "terms")
  terms2 = stats::terms(formula_no_sm)
  attr(terms2, "predvars") = attr(terms1, "predvars")
  terms2 = stats::delete.response(terms2)
  mf2 = model.frame(terms2, newdata, xlev=.getXlevels(terms1,mf1))
  X_gj = model.matrix(terms2, mf2)

  # Assemble Apred_is
  if( "fm_mesh_2d" %in% class(object$spatial_graph) ){
    predA_is = fm_evaluator( object$spatial_graph, loc=as.matrix(newdata[,object$internal$data_colnames$spatial]) )$proj$A
  }else if( "igraph" %in% class(object$spatial_graph) ) {
    Match = match( newdata[,object$internal$data_colnames$spatial], rownames(object$tmb_inputs$tmb_data$Adj) )
    if(any(is.na(Match))) stop("Check `spatial_graph` for SAR")
    predA_is = sparseMatrix( i=1:nrow(newdata), j=Match, x=rep(1,nrow(newdata)) )
  }else if( !is.null(object$internal$sem) ){
    predA_is = matrix(1, nrow=nrow(newdata), ncol=1)    # dgCMatrix
    predA_is = as(Matrix(predA_is),"dgCMatrix")
  } else {
    predA_is = Matrix(nrow=nrow(newdata), ncol=0)    # dgCMatrix
    predA_is = as(predA_is,"dgCMatrix")
  }
  predAtriplet = Matrix::mat2triplet(predA_is)

  # Turn of t_i and c_i when times and variables are missing, so that delta_k isn't built
  if( length(object$internal$times>0) ){
    t_g = match( newdata[,object$internal$data_colnames$time], object$internal$times )
  }else{ t_g = integer(0) }
  if( length(object$internal$variables>0) ){
    c_g = match( newdata[,object$internal$data_colnames$var], object$internal$variables )
  }else{ c_g = integer(0) }

  # Error checks
  tmb_data2 = object$tmb_inputs$tmb_data
  if( ncol(tmb_data2$X_ij) != ncol(X_gj) ) stop("Check X_gj")
  if( ncol(tmb_data2$Z_ik) != ncol(Z_gk) ) stop("Check Z_gk")

  # Swap in new predictive stuff
  tmb_data2$X_gj = X_gj
  tmb_data2$Z_gk = Z_gk
  tmb_data2$Agstc_zz = cbind(predAtriplet$i, predAtriplet$j, t_g[predAtriplet$i], c_g[predAtriplet$i]) - 1    # Triplet form, i, s, t
  tmb_data2$Axg_z = predAtriplet$x
  tmb_data2$t_g = t_g - 1 # Convert to CPP indexing
  tmb_data2$c_g = c_g - 1 # Convert to CPP indexing

  return( tmb_data2 )
}
