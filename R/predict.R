#' @title Predict using vector autoregressive spatio-temporal model
#'
#' @description Predicts values given new covariates using a tinyVAST model
#'
#' @inheritParams add_predictions
#'
#' @param remove_origdata Whether to speed `predict(.)` by elminating original data
#'        from TMB object, thereby speeding TMB object construction.  However, this
#'        also eliminates information about random-effect variance, and is not
#'        appropriate when requesting predictive standard errors or epsilon
#'        bias-correction.
#' @param ... Not used.
#'
#' @method predict tinyVAST
#' @export
predict.tinyVAST <-
function( object,
          newdata,
          remove_origdata = FALSE,
          ... ){

  # extract original X and Z
  if(missing(newdata)) newdata = object$data

  # Build new
  tmb_data2 = add_predictions( object = object,
                               newdata = newdata,
                               remove_origdata = remove_origdata )

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
#'     for each row of `newdata`.
#'
#' @export
integrate_output <-
function( object,
          newdata,
          bias.correct = TRUE,
          apply.epsilon = FALSE,
          intern = FALSE,
          W_gz,
          V_gz,
          area ){

  # extract original X and Z
  if(missing(newdata)) newdata = object$data
  # Build new .. object$data must be same as used for fitting to get SEs / skewness of random effects
  tmb_data2 = add_predictions( object = object,
                               newdata = newdata,
                               remove_origdata = isFALSE(apply.epsilon) & isFALSE(bias.correct) )

  # Area-expanded sum
  if(missing(W_gz)){
    # Default for area
    if(missing(area)){
      area = rep(1, nrow(newdata))
    }else if(length(area)==1){
      area = rep(area, nrow(newdata))
    }else if( length(area)!=nrow(newdata) ){
      stop("Check length of `area`")
    }
    tmb_data2$W_gz = cbind(area, 0)
  }else{
    tmb_data2$W_gz = W_gz
  }
  if(missing(V_gz)){
    tmb_data2$V_gz = cbind( rep(1,nrow(newdata)), 0 )
  }else{
    tmb_data2$V_gz = V_gz
  }

  # Abundance-weighted z
  #tmb_data2$W_gz = cbind( 1, newdata$x )
  #tmb_data2$V_gz = matrix(1, nrow=nrow(newdata), ncol=2)

  #
  tmb_par2 = object$internal$parlist
  if( isTRUE(apply.epsilon) ){
    tmb_par2$eps = 0
    inner.control = list(sparse=TRUE, lowrank=TRUE, trace=FALSE)
  }else{
    inner.control = list(sparse=TRUE, lowrank=FALSE, trace=FALSE)
  }

  # Rebuild object
  newobj = MakeADFun( data = tmb_data2,
                      parameters = tmb_par2,
                      map = object$tmb_inputs$tmb_map,
                      random = c("gamma_k","epsilon_stc","omega_sc"),
                      DLL = "tinyVAST",
                      intern = intern,
                      inner.control = inner.control,
                      profile = object$internal$control$profile )
  newobj$env$beSilent()

  # Run sdreport
  if( isTRUE(apply.epsilon) ){
    #Metric = newobj$report()$Metric
    grad = newobj$gr( newobj$par )[which(names(newobj$par)=="eps")]
    out = c( "Estimate"=NA, "Std. Error"=NA, "Est. (bias.correct)"=grad, "Std. (bias.correct)"=NA )
  }else if( isTRUE(bias.correct) ){
    newsd = sdreport( obj = newobj,
                      par.fixed = object$opt$par,
                      hessian.fixed = object$internal$Hess_fixed,
                      bias.correct = bias.correct )
    out = summary(newsd, "report")['Metric',]
  }else{
    rep = newobj$report()
    out = c( "Estimate"=rep$Metric, "Std. Error"=NA, "Est. (bias.correct)"=NA, "Std. (bias.correct)"=NA )
  }
  return(out)
}

#' @title Add predictions to data-list
#'
#' @description Given user-provided `newdata`, expand the object `tmb_data`
#'    to include predictions corresponding to those new observations
#'
#' @param object Output from [fit()].
#' @param newdata New data-frame of independent variables used to predict the response.
#' @param remove_origdata Whether to remove original-data to allow faster evaluation.
#'        \code{remove_origdata=TRUE} eliminates information about the distribution
#'        for random effects, and cannot be combined with epsilon bias-correction.
#'        WARNING:  feature is experimental and subject to change.
#'
#' @export
add_predictions <-
function( object,
          newdata,
          remove_origdata = FALSE ){

  tmb_data = object$tmb_inputs$tmb_data
  gam = object$gam_setup
  origdata = object$data
  data_colnames = object$internal$data_colnames

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
  if(is.null(Z_gk)) Z_gk = matrix( 0, nrow=nrow(newdata), ncol=ncol(tmb_data$Z_ik) )

  # Assemble predX
  formula_no_sm = remove_s_and_t2(gam$formula)
  mf1 = model.frame(formula_no_sm, origdata)
  terms1 = attr(mf1, "terms")
  terms2 = stats::terms(formula_no_sm)
  attr(terms2, "predvars") = attr(terms1, "predvars")
  terms2 = stats::delete.response(terms2)
  mf2 = model.frame(terms2, newdata, xlev=.getXlevels(terms1,mf1))
  X_gj = model.matrix(terms2, mf2)
  offset_g = model.offset(mf2)
  if(is.null(offset_g)) offset_g = rep(0,nrow(newdata))

  # Assemble Apred_is
  if( "fm_mesh_2d" %in% class(object$spatial_graph) ){
    predA_is = fm_evaluator( object$spatial_graph, loc=as.matrix(newdata[,data_colnames$spatial]) )$proj$A
  }else if( "igraph" %in% class(object$spatial_graph) ) {
    Match = match( newdata[,data_colnames$spatial], rownames(object$tmb_inputs$tmb_data$Adj) )
    if(any(is.na(Match))) stop("Check `spatial_graph` for SAR")
    predA_is = sparseMatrix( i=1:nrow(newdata), j=Match, x=rep(1,nrow(newdata)) )
  }else if( !is.null(object$internal$sem) ){
    predA_is = matrix(1, nrow=nrow(newdata), ncol=1)    # dgCMatrix
    #predA_is = as(Matrix(predA_is),"dgCMatrix")
    predA_is = as(Matrix(predA_is),"CsparseMatrix")
  } else {
    predA_is = Matrix(nrow=nrow(newdata), ncol=0)    # dgCMatrix
    #predA_is = as(predA_is,"dgCMatrix")
    predA_is = as(predA_is,"CsparseMatrix")
  }
  predAtriplet = Matrix::mat2triplet(predA_is)

  # Turn of t_i and c_i when times and variables are missing, so that delta_k isn't built
  if( length(object$internal$times) > 0 ){
    t_g = match( newdata[,data_colnames$time], object$internal$times )
  }else{ t_g = integer(0) }
  if( length(object$internal$variables) > 0 ){
    c_g = match( newdata[,data_colnames$var], object$internal$variables )
  }else{ c_g = integer(0) }

  #
  AepsilonG_zz = cbind(predAtriplet$i, predAtriplet$j, t_g[predAtriplet$i], c_g[predAtriplet$i])
  which_Arows = which(apply( AepsilonG_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) & any(x>0) ))
  which_Arows = which_Arows[ which(predAtriplet$x[which_Arows] > 0) ]
  if( nrow(object$internal$dsem_ram_output$ram)==0 ){
    which_Arows = numeric(0)
  }
  AepsilonG_zz = AepsilonG_zz[which_Arows,,drop=FALSE]
  AepsilonG_z = predAtriplet$x[which_Arows]

  #
  AomegaG_zz = cbind(predAtriplet$i, predAtriplet$j, c_g[predAtriplet$i])
  which_Arows = which(apply( AomegaG_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) ))
  which_Arows = which_Arows[ which(predAtriplet$x[which_Arows] > 0) ]
  if( nrow(object$internal$sem_ram_output$ram)==0 ){
    which_Arows = numeric(0)
  }
  AomegaG_zz = AomegaG_zz[which_Arows,,drop=FALSE]
  AomegaG_z = predAtriplet$x[which_Arows]

  # Telescope
  if( !(data_colnames$distribution %in% colnames(newdata)) ){
    if(nrow(object$internal$family_link)>1) stop("Must supply `dist` if using multiple `family_link` options")
    newdata = data.frame( newdata, matrix(rownames(object$internal$family_link)[1], nrow=nrow(newdata), ncol=1, dimnames=list(NULL,data_colnames$distribution)) )
  }

  # Build e_i ... always has length nrow(data)
  e_g = match( newdata[,data_colnames$distribution], rownames(object$internal$family_link) )

  #
  if( !(data_colnames$distribution %in% colnames(newdata)) ){
    newdata = cbind( newdata, matrix(1, nrow=nrow(newdata), ncol=1, dimnames=list(NULL,data_colnames$distribution)) )
  }

  # Error checks
  tmb_data2 = object$tmb_inputs$tmb_data
  if( ncol(tmb_data2$X_ij) != ncol(X_gj) ) stop("Check X_gj")
  if( ncol(tmb_data2$Z_ik) != ncol(Z_gk) ) stop("Check Z_gk")

  # Swap in new predictive stuff
  tmb_data2$X_gj = X_gj
  tmb_data2$Z_gk = Z_gk
  tmb_data2$AepsilonG_zz = AepsilonG_zz - 1    # Triplet form, i, s, t
  tmb_data2$AepsilonG_z = AepsilonG_z
  tmb_data2$AomegaG_zz = AomegaG_zz - 1    # Triplet form, i, s, t
  tmb_data2$AomegaG_z = AomegaG_z
  tmb_data2$t_g = t_g - 1 # Convert to CPP indexing
  tmb_data2$c_g = c_g - 1 # Convert to CPP indexing
  tmb_data2$offset_g = offset_g
  tmb_data2$e_g = e_g - 1 # -1 to convert to CPP index

  # Simplify by eliminating observations ... experimental
  if( isTRUE(remove_origdata) ){
    warning("`remove_origdata` is experimental")
    tmb_data2$y_i = tmb_data2$y_i[numeric(0)]
    tmb_data2$X_ij = tmb_data2$X_ij[numeric(0),,drop=FALSE]
    tmb_data2$Z_ik = tmb_data2$Z_ik[numeric(0),,drop=FALSE]
    tmb_data2$t_i = tmb_data2$t_i[numeric(0)]
    tmb_data2$c_i = tmb_data2$c_i[numeric(0)]
    tmb_data2$e_i = tmb_data2$e_i[numeric(0)]
    tmb_data2$Aepsilon_zz = tmb_data2$Aepsilon_zz[numeric(0),,drop=FALSE]
    tmb_data2$Aepsilon_z = tmb_data2$Aepsilon_z[numeric(0)]
    tmb_data2$Aomega_zz = tmb_data2$Aomega_zz[numeric(0),,drop=FALSE]
    tmb_data2$Aomega_z = tmb_data2$Aomega_z[numeric(0)]
  }

  # Check for obvious issues ... no NAs except in RAMstart
  if( any(is.na(tmb_data2[-match("RAMstart",names(tmb_data2))])) ){
    stop("Check output of `add_predictions` for NAs")
  }
  # Check for obvious issues ... length of inputs
  if( !all( sapply(tmb_data2[c('t_g','c_g','e_g','offset_g')],FUN=length) %in% c(0,nrow(newdata)) ) ){
    stop("Check output of `add_predictions` for variables with unexpected length")
  }
  if( any( sapply(tmb_data2[c('X_gj','Z_gk')],FUN=nrow) != nrow(newdata) ) ){
    stop("Check output of `add_predictions` for NAs")
  }

  return( tmb_data2 )
}
