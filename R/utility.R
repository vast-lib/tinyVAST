#' @title Multivariate Normal Random Deviates using Sparse Precision
#'
#' @description This function provides a random number generator for
#'              the multivariate normal distribution with mean equal
#'              to `mean` and sparse precision matrix `Q`.
#'
#' @param Q sparse precision (inverse-covariance) matrix.
#' @param n number of observations.
#' @param mean mean vector.
#'
#' @return a matrix with dimension \code{length(mean)} by
#'         \code{n}, containing realized draws from the specified
#'         mean and precision
#'
#' @export
rmvnorm_prec <-
function( Q,
          n = 1,
          mean = rep(0,nrow(Q)) ) {

  # Simulate values
  z0 = matrix( rnorm(length(mean) * n), ncol=n)

  # Q = t(P) * L * t(L) * P
  L = Matrix::Cholesky(Q, super=TRUE)

  # Calculate t(P) * solve(t(L)) * z0 in two steps
  z = Matrix::solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = Matrix::solve(L, z, system = "Pt") # z = Pt    * z
  return(mean + as.matrix(z))
}

#' @title Rotate factors to match Principal-Components Analysis
#'
#' @description Rotate lower-triangle loadings matrix
#'              to order factors from largest to smallest variance.
#'
#' @param L_tf Loadings matrix with dimension \eqn{T \times F}.
#' @param x_sf Spatial response with dimensions \eqn{S \times F}.
#' @param order Options for resolving label-switching via reflecting
#'        each factor to achieve a given order across dimension \eqn{T}.
#'
#' @return List containing the rotated loadings \code{L_tf},
#'         the inverse-rotated response matrix \code{x_sf},
#'         and the rotation \code{H}
#'
#' @export
rotate_pca <-
function( L_tf,
          x_sf = matrix(0, nrow=0, ncol=ncol(L_tf)),
          order = c("none","increasing","decreasing") ){

  # Eigen-decomposition
  Cov_tmp = L_tf %*% t(L_tf)
  Cov_tmp = 0.5*Cov_tmp + 0.5*t(Cov_tmp) # Ensure symmetric
  Eigen = eigen(Cov_tmp)

  # Desired loadings matrix
  L_tf_rot = (Eigen$vectors%*%diag(sqrt(Eigen$values)))[,seq_len(ncol(L_tf)),drop=FALSE]

  # My new factors
  H = pseudoinverse(L_tf_rot) %*% L_tf
  x_sf = t(H %*% t(x_sf))

  # Get all loadings matrices to be increasing or decreasing
  order = match.arg(order)
  if( !is.null(order) ){
    for( f in seq_len(ncol(L_tf)) ){
      Lm = lm( L_tf_rot[,f] ~ 1 + I(seq_len(nrow(L_tf))) )
      Sign = sign(Lm$coef[2]) * ifelse(order=="decreasing", -1, 1)
      L_tf_rot[,f] = L_tf_rot[,f] * Sign
      x_sf[,f] = x_sf[,f] * Sign
    }
  }

  # return
  out = list( "L_tf"=L_tf_rot, "x_sf"=x_sf, "H"=H)
  return(out)
}


#' @title Parse path
#'
#' @description \code{parse_path} is copied from \code{sem::parse.path}
#'
#' @param path character string indicating a one-headed or two-headed path
#'        in a structural equation model
#'
#' @details
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @return Tagged-list defining variables and direction for a specified path coefficient
parse_path <-
function( path ){
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  out = list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
  return(out)
}


#' @title Classify variables path
#'
#' @description \code{classify_variables} is copied from \code{sem:::classifyVariables}
#'
#' @param model syntax for structural equation model
#'
#' @details
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @return Tagged-list defining exogenous and endogenous variables
classify_variables <-
function( model ){

    variables <- logical(0)
    for (paths in model[, 1]) {
        vars <- gsub(pattern=" ", replacement="", x=paths)
        vars <- sub("-*>", "->", sub("<-*", "<-", vars))
        if (grepl("<->", vars)) {
            vars <- strsplit(vars, "<->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
        }
        else if (grepl("->", vars)) {
            vars <- strsplit(vars, "->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            variables[vars[2]] <- TRUE
        }
        else if (grepl("<-", vars)) {
            vars <- strsplit(vars, "<-")[[1]]
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
            variables[vars[1]] <- TRUE
        }
        else stop("incorrectly specified model", call. = FALSE)
    }
    list(endogenous = names(variables[variables]), exogenous = names(variables[!variables]))
}

#' Reload a previously fitted model
#'
#' \code{reload_model} allows a user to save a fitted model, reload it in a new
#'      R terminal, and then relink the DLLs so that it functions as expected.
#'
#' @param x Output from \code{\link{tinyVAST}}, potentially with DLLs not linked
#' @param check_gradient Whether to check the gradients of the reloaded model
#'
#' @return Output from \code{\link{tinyVAST}} with DLLs relinked
#'
#' @export
reload_model <-
function( x,
          check_gradient = TRUE ){

  # Retape
  #obj = x$obj
  #obj$retape()

  # rebuild
  obj = MakeADFun( data = x$tmb_inputs$tmb_data,
                      parameters = x$internal$parlist,
                      map = x$tmb_inputs$tmb_map,
                      random = x$tmb_inputs$tmb_random,
                      DLL = "tinyVAST",
                      profile = x$internal$control$profile )
  obj$env$beSilent()

  # Ensure that last.par and last.par.best are right
  nll_new = obj$fn(x$opt$par)
  if( abs(nll_new-x$opt$obj) > 0.01 ){
    stop("Model fit is not identical to recorded value: Something is not working as expected")
  }

  # Check gradient
  if( isTRUE(check_gradient) ){
    Gr = obj$gr(x$opt$par)
    if( max(abs(Gr))>1 ){
      warning("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": does not seem converged")
    }else if( max(abs(Gr))>0.01 ){
      warning("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": might not be converged")
    }else{
      message("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": No evidence of non-convergence")
    }
  }

  x$obj = obj
  return(x)
}



#' @title
#' Sample from predictive distribution of a variable
#'
#' @description
#' \code{sample_variable} samples from the joint distribution of random and fixed effects to approximate the predictive distribution for a variable
#'
#' Using \code{sample_fixed=TRUE} (the default) in \code{\link{sample_variable}} propagates variance in both fixed and random effects, while
#'       using \code{sample_fixed=FALSE} does not.
#'       Sampling fixed effects will sometimes cause numerical under- or overflow (i.e., output values of \code{NA}) in cases when
#'       variance parameters are estimated imprecisely.  In these cases, the multivariate normal approximation being used is a poor
#'       representation of the tail probabilities, and results in some samples with implausibly high (or negative) variances,
#'       such that the associated random effects then have implausibly high magnitude.
#'
#' @param object output from `\code{tinyVAST()}`
#' @param newdata data frame of new data, used to sample model components for predictions e.g., \code{mu_g}
#' @param variable_name name of variable available in report using \code{Obj$report()} or parameters using \code{Obj$env$parList()}
#' @param n_samples number of samples from the joint predictive distribution for fixed and random effects.  Default is 100, which is slow.
#' @param seed integer used to set random-number seed when sampling variables, as passed to \code{set.seed(.)}
#' @param sample_fixed whether to sample fixed and random effects, \code{sample_fixed=TRUE} as by default, or just sample random effects, \code{sample_fixed=FALSE}
#'
#' @return
#' A matrix with a row for each \code{data} supplied during fitting, and
#' \code{n_samples} columns, where each column in a vector of samples
#' for a requested quantity given sampled uncertainty in fixed and/or random effects
#'
#' @examples
#'  set.seed(101)
#'  x = runif(n = 100, min = 0, max = 2*pi)
#'  y = 1 + sin(x) + 0.1 * rnorm(100)
#'
#'  # Do fit with getJointPrecision=TRUE
#'  fit = tinyVAST( formula = y ~ s(x),
#'                  data = data.frame(x=x,y=y) )
#'
#'  # samples from distribution for the mean
#'  # excluding fixed effects due to CRAN checks
#'  samples = sample_variable(fit, sample_fixed = FALSE)
#'
#' @export
sample_variable <-
function( object,
          newdata = NULL,
          variable_name = "mu_i",
          n_samples = 100,
          sample_fixed = TRUE,
          seed = 123456 ){

  # Rebuild object with newdata
  if( !is.null(newdata) ){
    tmb_data2 = add_predictions( object = object, newdata = newdata )
    newobj = MakeADFun( data = tmb_data2,
                        parameters = object$internal$parlist,
                        map = object$tmb_inputs$tmb_map,
                        random = object$tmb_inputs$tmb_random,
                        profile = object$internal$control$profile,
                        DLL = "tinyVAST" )
    object = list(
      obj = newobj,
      sdrep = object$sdrep,
      rep = newobj$rep()
    )
  }

  # Combine Report and ParHat, and check for issues
  ParHat = object$obj$env$parList()
  Intersect = intersect(names(object$rep), names(ParHat))
  if( isFALSE(all.equal(object$rep[Intersect],ParHat[Intersect])) ){
    stop("Duplicate entries in `Obj$report()` and `Obj$env$parList()` are not identical when calling `sample_variable`")
  }
  Output = c( object$rep, ParHat )
  # Check that variable_name is available
  if( isFALSE(variable_name %in% names(Output)) ){
    stop( variable_name, " not found in `Obj$report()` or `Obj$env$parList()`; please choose check your requested variable name from available list: ", paste(names(Output),collapse=", ") )
  }

  #### Local function
  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims, seed) {
    set.seed(seed)
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- Cholesky(prec, super=TRUE)
    z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    return(mu + z)
  }

  # Sample from joint distribution
  if( sample_fixed==TRUE ){
    # Informative error messages
    if( !("jointPrecision" %in% names(object$sdrep)) ){
      stop("jointPrecision not present in object$sdrep; please re-run with `getJointPrecision=TRUE`")
    }
    u_zr = rmvnorm_prec( mu=object$obj$env$last.par.best, prec=object$sdrep$jointPrecision, n.sims=n_samples, seed=seed)
    # apply( u_zr, MARGIN=2, FUN=function(vec){sum(abs(vec)==Inf)})
    # u_zr[-object$obj$env$random,1]
  }else{
    u_zr = object$obj$env$last.par.best %o% rep(1, n_samples)
    Q = object$obj$env$spHess( random = TRUE, par = object$obj$env$last.par.best )
    MC = rmvnorm_prec( mu=object$obj$env$last.par.best[object$obj$env$lrandom()], prec=Q, n.sims=n_samples, seed=seed)
    #MC = object$obj$env$MC( keep=TRUE, n=n_samples, antithetic=FALSE )
    #u_zr[object$obj$env$random,] = attr(MC, "samples")
    u_zr[object$obj$env$lrandom(),] = MC
  }

  # Extract variable for each sample
  message( "# Obtaining samples from predictive distribution for variable ", variable_name )
  for( rI in 1:n_samples ){
    if( rI%%max(1,floor(n_samples/10)) == 0 ){
      message( "  Finished sample ", rI, " of ",n_samples )
    }
    Report = object$obj$report( par=u_zr[,rI] )
    ParHat = object$obj$env$parList( x=u_zr[,rI][object$obj$env$lfixed()], par=u_zr[,rI] )
    if( isFALSE(all.equal(object$rep[Intersect],ParHat[Intersect])) ){
      stop("Duplicate entries in `object$obj$report()` and `object$obj$env$parList()` are not identical when calling `sample_variable`")
    }
    Output = c( Report, ParHat )
    Var = Output[[variable_name]]
    if(is.vector(Var)) Var = as.array(Var)
    if(rI==1) Var_zr = Var
    if(rI>=2){
      Var_zr = abind( Var_zr, Var, along=length(dim(Var))+1 )
    }
  }

  # Return
  return( Var_zr )
}

#' @title
#' Calculate deviance explained
#'
#' @description
#' \code{deviance_explained} fits a null model, calculates the deviance relative to 
#'       a saturated model for both the original and the null model, and uses these
#'       to calculate the proportion of deviance explained.                   
#'
#' This implementation conditions upon the maximum likelihood estimate of fixed effects
#'      and the empirical Bayes ("plug-in") prediction of random effects.  It can 
#'      be described as "conditional deviance explained". A state-space model that 
#'      estimates measurement error variance approaching zero (i.e., collapses to 
#'      a process-error-only model) will have a conditional deviance explained
#'      that approaches 1.0
#'
#' @param x output from `\code{tinyVAST()}`
#' @param null_formula formula for the null model.  If missing, it uses
#'        \code{null_formula = response ~ 1}. For multivariate models, it 
#'        might make sense to use \code{null_formula = response ~ category}
#' @param null_delta_formula formula for the null model for the delta component.  
#'        If missing, it uses
#'        \code{null_formula = response ~ 1}. For multivariate models, it 
#'        might make sense to use \code{null_delta_formula = response ~ category}
#'
#' @return the proportion of conditional deviance explained.
#'
#' @export
deviance_explained <-
function( x, 
          null_formula,
          null_delta_formula = ~ 1 ){

  # null formula
  if(missing(null_formula)){
    null_formula = update.formula( old = x$formula, 
                                   new = . ~ 1 )
  }
  
  # Edit control to save time
  # Make sure to set calculate_deviance_explained = FALSE to avoid recursion in tinyVAST(.) call
  control = x$internal$control
    control$getsd = FALSE
    control$silent = TRUE
    control$trace = 0
    control$verbose = FALSE
    control$calculate_deviance_explained = FALSE
    control$tmb_par = NULL
    control$tmb_map = NULL

  # 
  control_initial = control
    control_initial$nlminb_loops = 0
    control_initial$newton_loops = 0
  
  # Run null model to check that some parameters remain
  null_fit = tinyVAST( data = x$data,
                       formula = null_formula, 
                       control = control_initial,
                       family = x$internal$family,
                       space_columns = x$internal$space_columns,
                       time_column = x$internal$time_column,
                       variable_column = x$internal$variable_column,
                       times = x$internal$times,
                       weights = x$internal$weights,
                       variables = x$internal$variables,
                       delta_options = list( formula = null_delta_formula ),
                       distribution_column = x$internal$distribution_column) 
  null_obj = null_fit$obj

  # Run if some parameters remain
  if( length(null_obj$par)>0 ){
    # null model
    #null_fit = tinyVAST( data = x$data,
    #                     formula = null_formula,
    #                     control = control,
    #                     family = x$internal$family,
    #                     space_columns = x$internal$space_columns,
    #                     time_column = x$internal$time_column,
    #                     variable_column = x$internal$variable_column,
    #                     times = x$internal$times,
    #                     variables = x$internal$variables,
    #                     delta_options = list( delta_formula = null_delta_formula ),
    #                     distribution_column = x$internal$distribution_column)
    null_opt = nlminb( start = null_obj$par,
                       objective = null_obj$fn,
                       gradient = null_obj$gr )
  }else{
    #devexpl = "Not calculating deviance explained. Please try again without using `control$profile`"
    null_obj$fn(null_obj$par)
  }
  null_rep = null_obj$report()

  # Calculate deviance explained
  devexpl = 1 - x$rep$deviance / null_rep$deviance
  if( isTRUE(devexpl<0) | isTRUE(devexpl>1) ){
    warning("Problem detected: deviance explained should be between 0 and 1")
  }
  if( is.na(devexpl) ){
    warning("Deviance explained is NA, probably because it's not implemented for the family used")
  }

  # Return
  return( devexpl )
}

# from glmmTMB:
named_list <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(nm <- names(L))) {
    nm <- snm
  }
  if (any(nonames <- nm == "")) {
    nm[nonames] <- snm[nonames]
  }
  stats::setNames(L, nm)
}



#' @title Extract covariance
#'
#' @description Extract the covariance resulting from a specified path structure
#' and estimated parameters for a SEM or DSEM term in tinyVAST
#'
#' @param object Output from \code{\link{tinyVAST}}
#' @param what Which SEM or DSEM term to extract
#' @param pred Extract the term \code{what} for which linear predictor
#' @param n_times The number of times to include when calculating covariance for a DSEM
#'        component, i.e., \code{time_term} or \code{spacetime_term}.  If missing, the default
#'        is to use the one more than the maximum specified lag (e.g., \code{n_times=2}
#'        by default when the maximum \code{lag=1})
#'
#' @importFrom dsem make_matrices
#'
#' @details
#' tinyVAST constructs the covariance from specified path structure and estimated parameters
#'
#' @return
#' The covariance matrix among variables
#'
#' @examples
#' # Extract covariance for spatial factor analysis (too slow for CRAN)
#' \donttest{
#' # Simulate settings
#' set.seed(101)
#' theta_xy = 0.4
#' n_x = n_y = 10
#' n_c = 3           # Number of species
#' n_f = 1           # Number of factors
#' rho = 0.8
#' resid_sd = 0.5
#'
#' # Simulate GMRFs
#' R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
#' R_ss = kronecker(X=R_s, Y=R_s)
#' delta_fs = mvtnorm::rmvnorm(n_c, sigma=R_ss )
#'
#' # Simulate loadings for two factors
#' L_cf = matrix( rnorm(n_c^2), nrow=n_c )
#' L_cf[,seq(from=n_f+1, to=n_c)] = 0
#' L_cf = L_cf + resid_sd * diag(n_c)
#'
#' # Simulate correlated densities
#' d_cs = L_cf %*% delta_fs
#'
#' # Shape into longform data-frame and add error
#' Data = data.frame( expand.grid(species=1:n_c, x=1:n_x, y=1:n_y),
#'                    "var"="logn", "z"=exp(as.vector(d_cs)) )
#' Data$n = rnorm( n=nrow(Data), mean=Data$z, sd=1 )
#'
#' # make mesh
#' mesh = fmesher::fm_mesh_2d( Data[,c('x','y')] )
#'
#' # Specify factor model with two factors and additional independent variance with shared SD
#' sem = "
#'   # Loadings matrix
#'   f1 -> 1, l1
#'   f1 -> 2, l2
#'   f1 -> 3, l3
#'
#'   # Factor variance = 1
#'   f1 <-> f1, NA, 1
#'
#'   # Shared residual variance
#'   1 <-> 1, sd, 1
#'   2 <-> 2, sd, 1
#'   3 <-> 3, sd, 1
#' "
#'
#' # fit model
#' out = tinyVAST( space_term = sem,
#'            data = Data,
#'            formula = n ~ 0 + factor(species),
#'            spatial_domain = mesh,
#'            variables = c( "f1", "f2", 1:n_c ),
#'            space_columns = c("x","y"),
#'            variable_column = "species",
#'            time_column = "time",
#'            distribution_column = "dist" )
#'
#' # Extract covariance among species and factors, where
#' # estimated covariance is obtained by ignoring factors
#' V = term_covariance( out, what = "space_term", pred = "one" )
#' }
#'
#' @export
term_covariance <-
function( object,
          what = c("space_term", "time_term", "spacetime_term"),
          pred = c("one", "two"),
          n_times = NULL ){

  pred = match.arg(pred)
  what = match.arg(what)

  # Extract stuff
  variables = object$internal$variables
  if(what=="space_term"){
    if( pred == "one" ){
      parname = "theta_z"
      ram = object$internal$space_term_ram$output$ram
    }else{
      parname = "theta2_z"
      ram = object$internal$delta_space_term_ram$output$ram
    }
    model = data.frame( "direction" = ram$heads,
                        "lag" = 0,
                        "first" = variables[ram$from],
                        "second" = variables[ram$to],
                        "start" = ram$start,
                        "parameter" = ram$parameter )

    if(is.null(n_times)){
      times = 1
    }else{
      times = seq_len(n_times)
    }
    var_names = variables
  }
  if(what=="time_term"){
    if( pred == "one" ){
      parname = "nu_z"
      model = as.data.frame(object$internal$time_term_ram$output$model)
    }else{
      parname = "nu2_z"
      model = as.data.frame(object$internal$delta_time_term_ram$output$model)
    }
    if(is.null(n_times)){
      times = seq_len( max(as.numeric(model$lag)+1) )
    }else{
      times = seq_len(n_times)
    }
    var_names = apply( expand.grid(paste0("_lag",times-1),variables),
                       FUN=function(vec){paste0(vec[2],vec[1])}, MARGIN=1 )
  }
  if(what=="spacetime_term"){
    if( pred == "one" ){
      parname = "beta_z"
      model = as.data.frame(object$internal$spacetime_term_ram$output$model)
    }else{
      parname = "beta2_z"
      model = as.data.frame(object$internal$delta_spacetime_term_ram$output$model)
    }
    if(is.null(n_times)){
      times = seq_len( max(as.numeric(model$lag)+1) )
    }else{
      times = seq_len(n_times)
    }
    var_names = apply( expand.grid(paste0("_lag",times-1),variables),
                       FUN=function(vec){paste0(vec[2],vec[1])}, MARGIN=1 )
  }

  if( nrow(model) > 0 ){
    # Unpack stuff
    if(is.null(object$internal$parlist)){
      object$internal$parlist = object$obj$env$parList()
    }
    # Extract path matrix
    matrices = dsem::make_matrices(
      beta_p = object$internal$parlist[[parname]],
      model = model,
      times = times,
      variables = variables
    )
    P_kk = matrices$P_kk
    G_kk = matrices$G_kk

    # Calcualte covariance (not precision, in case G is not invertible)
    IminusP_kk = Diagonal(n = nrow(P_kk)) - P_kk
    #Q_kk = t(IminusP_kk) %*% solve(t(G_kk) %*% G_kk) %*% IminusP_kk
    V_kk = solve(IminusP_kk) %*% t(G_kk) %*% G_kk %*% solve(t(IminusP_kk))
    dimnames(V_kk) = list( var_names, var_names )
  }else{
    V_kk = "selected term not present in model"
  }

  # Return
  return(V_kk)
}

