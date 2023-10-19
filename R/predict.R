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

  # Assemble predZ
  predZ = lapply( seq_along(gam$smooth),
                  FUN = \(x) PredictMat(gam$smooth[[x]], data=newdata) )
  predZ = do.call( cbind, predZ )
  if(is.null(predZ)) predZ = matrix(0, nrow=nrow(newdata),ncol=ncol(tmb_data$Z))

  # Assemble predX
  formula_no_sm = remove_s_and_t2(gam$formula)
  mf1 = model.frame(formula_no_sm, origdata)
  terms1 = attr(mf1, "terms")
  terms2 = stats::terms(formula_no_sm)
  attr(terms2, "predvars") = attr(terms1, "predvars")
  terms2 = stats::delete.response(terms2)
  mf2 = model.frame(terms2, newdata, xlev=.getXlevels(terms1,mf1))
  predX = model.matrix(terms2, mf2)

  # Assemble Apred_is
  if( "fm_mesh_2d" %in% class(object$spatial_graph) ){
    predA_is = fm_evaluator( object$spatial_graph, loc=as.matrix(newdata[,object$internal$data_colnames$spatial]) )$proj$A
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
    predt_i = match( newdata[,object$internal$data_colnames$time], object$internal$times )
  }else{ predt_i = integer(0) }
  if( length(object$internal$variables>0) ){
    predc_i = match( newdata[,object$internal$data_colnames$var], object$internal$variables )
  }else{ predc_i = integer(0) }

  #
  #object$obj$env$data$predX = predX          # object$tmb_inputs$tmb_data$X
  #object$obj$env$data$predZ = predZ          # object$tmb_inputs$tmb_data$Z
  #object$obj$env$data$Apred_is = Apred_is    # object$tmb_inputs$tmb_data$Z
  #out = object$obj$report()$mu_pred

  # Error checks
  tmb_data2 = object$tmb_inputs$tmb_data
  if( ncol(tmb_data2$X) != ncol(predX) ) stop("Check predX")
  if( ncol(tmb_data2$Z) != ncol(predZ) ) stop("Check predZ")

  # Swap in new predictive stuff
  tmb_data2$predX = predX
  tmb_data2$predZ = predZ
  tmb_data2$predAistc = cbind(predAtriplet$i, predAtriplet$j, predt_i[predAtriplet$i], predc_i[predAtriplet$i]) - 1    # Triplet form, i, s, t
  tmb_data2$predAx = predAtriplet$x
  tmb_data2$predt_i = predt_i
  tmb_data2$predc_i = predc_i

  # Rebuild object
  newobj = MakeADFun( data = tmb_data2,
                      parameters = object$internal$parlist,
                      map = object$tmb_inputs$tmb_map,
                      random = c("gamma","epsilon_stc"),
                      DLL = "tinyVAST" )
  out = newobj$report()$mu_pred

  return(out)
}
