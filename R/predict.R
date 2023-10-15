#' @title Predict using vector autoregressive spatio-temporal model
#'
#' @description Predicts values given new covariates using a tinyVAST model
#'
#' @method predict tinyVAST
#' @export
predict.tinyVAST <-
function( object,
          newdata ){

  # extract original X and Z
  origdata = eval(object$call$data)
  if(missing(newdata)) newdata = origdata
  tmb_data = object$tmb_inputs$tmb_data
  gam = object$gam_setup

  # Check newdata for missing variables and/or factors
  #response_name = as.character(gam$formula)[2]
  #pred_set = setdiff( colnames(gam$mf), response_name )
  pred_set = strsplit( as.character(object$gam_setup$pred.formula)[2], split=" + ", fixed=TRUE )[[1]]
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
  formula_no_sm <- remove_s_and_t2(gam$formula)
  mf <- model.frame(formula_no_sm, origdata)
  terms = attr(mf, "terms")
  f2 <- remove_s_and_t2(gam$formula)
  tt <- stats::terms(f2)
  attr(tt, "predvars") <- attr(terms, "predvars")
  Terms <- stats::delete.response(tt)
  mf2 <- model.frame(Terms, newdata, xlev = .getXlevels(Terms,mf))
  Xpred <- model.matrix(Terms, mf2)

  #
  object$obj$env$data$Xpred = Xpred # object$tmb_inputs$tmb_data$X
  object$obj$env$data$Zpred = Zpred # object$tmb_inputs$tmb_data$Z
  out = object$obj$report()$mu_pred
  return(out)
}
