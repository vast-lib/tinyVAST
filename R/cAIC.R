#' Calculate conditional AIC
#'
#' Calculates the conditional Akaike Information criterion (cAIC).
#'
#' @param object Output from [tinyVAST()].
#'
#' @details cAIC is designed to optimize the expected out-of-sample predictive
#' performance for new data that share the same random effects as the in-sample
#' (fitted) data, e.g., spatial interpolation.  In this sense, it should be a
#' fast approximation to optimizing the model structure based on k-fold
#' cross-validation.
#'
#' By contrast, [AIC()] calculates the marginal Akaike Information Criterion,
#' which is designed to optimize expected predictive performance for new data
#' that have new random effects, e.g., extrapolation, or inference about
#' generative parameters.
#'
#' Both cAIC and EDF are calculated using Eq. 6 of Zheng, Cadigan, and Thorson
#' (2024).
#'
#' For models that include profiled fixed effects, these profiles are turned
#' off.
#'
#' @return
#' cAIC value
#'
#' @references
#' Zheng, N., Cadigan, N., & Thorson, J. T. (2024).
#' A note on numerical evaluation of conditional Akaike information for
#' nonlinear mixed-effects models (arXiv:2411.14185). arXiv.
#' \doi{10.48550/arXiv.2411.14185}
#'
#' @examples
#' data( red_snapper )
#' red_snapper = droplevels(subset(red_snapper, Data_type=="Biomass_KG"))
#'
#' # Define mesh
#' mesh = fmesher::fm_mesh_2d( red_snapper[,c('Lon','Lat')],
#'                            cutoff = 1 )
#'
#' # define formula with a catchability covariate for gear
#' formula = Response_variable ~ factor(Year) + offset(log(AreaSwept_km2))
#'
#' # make variable column
#' red_snapper$var = "logdens"
#
#' # fit using tinyVAST
#' fit = tinyVAST( data = red_snapper,
#'                 formula = formula,
#'                 sem = "logdens <-> logdens, sd_space",
#'                 space_columns = c("Lon",'Lat'),
#'                 spatial_graph = mesh,
#'                 family = tweedie(link="log"),
#'                 variable_column = "var",
#'                 control = tinyVASTcontrol( getsd = FALSE,
#'                                            profile = "alpha_j" ) )
#'
#' cAIC(fit)
#'
cAIC <-
function( object ){

  #what = match.arg(what)
  require(Matrix)
  data = object$tmb_inputs$tmb_data

  # Make sure profile = NULL
  if( is.null(object$internal$control$profile) ){
    obj = object$obj
  }else{
    obj = TMB::MakeADFun( data = data,
                        parameters = object$internal$parlist,
                        map = object$tmb_inputs$tmb_map,
                        random = object$tmb_inputs$tmb_random,
                        DLL = "tinyVAST",
                        profile = NULL )
  }

  # Make obj_new
  data$weights_i[] = 0
  obj_new = TMB::MakeADFun( data = data,
                      parameters = object$internal$parlist,
                      map = object$tmb_inputs$tmb_map,
                      random = object$tmb_inputs$tmb_random,
                      DLL = "tinyVAST",
                      profile = NULL )

  #
  par = obj$env$parList()
  parDataMode <- obj$env$last.par
  indx = obj$env$lrandom()
  q = sum(indx)
  p = length(object$opt$par)

  ## use - for Hess because model returns negative loglikelihood;
  #cov_Psi_inv = -Hess_new[indx,indx]; ## this is the marginal prec mat of REs;
  Hess_new = -Matrix(obj_new$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  Hess_new = Hess_new[indx,indx]

  ## Joint hessian etc
  Hess = -Matrix(obj$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  Hess = Hess[indx,indx]
  negEDF = diag(solve(Hess, Hess_new))

  #if(what == "cAIC"){
    jnll = obj$env$f(parDataMode)
    cnll = jnll - obj_new$env$f(parDataMode)
    cAIC = 2*cnll + 2*(p+q) - 2*sum(negEDF)
    return(cAIC)
  #}
  #if(what == "EDF"){
  #  group = colnames(object$tmb_inputs$tmb_data$Z_ik)
  #  group = sapply(group, FUN=\(char)strsplit(char,split=".",fixed=TRUE)[[1]][1])
  #  group = factor(group)
  #  EDF = tapply(negEDF,INDEX=group,FUN=length) - tapply(negEDF,INDEX=group,FUN=sum)
  #  return(EDF)
  #}
}
