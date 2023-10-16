#' @title Predict using vector autoregressive spatio-temporal model
#'
#' @description Predicts values given new covariates using a tinyVAST model
#'
#' @method predict tinyVAST
#' @export
predict.tinyVAST <-
function( object,
          newdata,
          ... ){

  # extract original X and Z
  origdata = object$data
  if(missing(newdata)) newdata = origdata
  tmb_data = object$tmb_inputs$tmb_data
  gam = object$gam_setup

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

  # Assemble Zpred
  Zpred = lapply( seq_along(gam$smooth),
                  FUN = \(x) PredictMat(gam$smooth[[x]], data=newdata) )
  Zpred = do.call( cbind, Zpred )
  if(is.null(Zpred)) Zpred = matrix(0, nrow=nrow(newdata),ncol=ncol(tmb_data$Z))

  # Assemble Xpred
  formula_no_sm = remove_s_and_t2(gam$formula)
  mf1 = model.frame(formula_no_sm, origdata)
  terms1 = attr(mf1, "terms")
  terms2 = stats::terms(formula_no_sm)
  attr(terms2, "predvars") = attr(terms1, "predvars")
  terms2 = stats::delete.response(terms2)
  mf2 = model.frame(terms2, newdata, xlev=.getXlevels(terms1,mf1))
  Xpred = model.matrix(terms2, mf2)

  # Assemble Apred_is
  if( "fm_mesh_2d" %in% class(object$spatial_graph) ){
    Apred_is = fm_evaluator( object$spatial_graph, loc=as.matrix(newdata[,object$data_colnames$spatial]) )$proj$A
    Apred_is = as.matrix(Apred_is)  # Apred_is must be dense!
  } else {
    Apred_is = matrix(0, nrow=nrow(newdata), ncol=0)
  }

  #
  object$obj$env$data$Xpred = Xpred          # object$tmb_inputs$tmb_data$X
  object$obj$env$data$Zpred = Zpred          # object$tmb_inputs$tmb_data$Z
  object$obj$env$data$Apred_is = Apred_is    # object$tmb_inputs$tmb_data$Z
  out = object$obj$report()$mu_pred
  return(out)
}
