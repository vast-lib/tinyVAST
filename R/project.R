
#' @title Conditional simulation from a GMRF
#'
#' @description
#' Generates samples from a Gaussian Markov random field (GMRF) conditional upon
#' fixed values for some elements.
#'
#' @param Q precision for a zero-centered GMRF.
#' @param observed_idx integer vector listing rows of \code{Q} corresponding to
#'        fixed measurements
#' @param x_obs numeric vector with fixed values for indices \code{observed_idx}
#' @param n_sims integer listing number of simulated values
#' @param what Whether to simulate from the conditional GMRF, or predict the mean
#'        and precision
#'
#' @return
#' A matrix with \code{n_sims} columns and a row for every row of \code{Q} not in
#' \code{observed_idx}, with simulations for those rows
#'
#' @export
conditional_gmrf <-
function( Q,
          observed_idx,
          x_obs,
          n_sims = 1,
          what = c("simulate","predict") ){

  # Required libraries
  #library(Matrix)
  what = match.arg(what)

  # Error checks
  if( !all(observed_idx %in% seq_len(nrow(Q))) ){
    stop("Check `observed_idx` in `conditional_gmrf")
  }
  if( length(observed_idx) != length(x_obs) ){
    stop("Check length of `observed_idx` and `x_obs`")
  }
  if( any(is.na(x_obs)) ){
    stop("`x_obs` cannot include NA values")
  }

  # Calculate conditional mean and variance
  predict_conditional_gmrf <- function(Q, observed_idx, x_obs) {
    all_idx <- seq_len(nrow(Q))
    unobserved_idx <- setdiff(all_idx, observed_idx)

    # Partition Q
    #Q_oo <- Q[observed_idx, observed_idx, drop = FALSE]
    Q_ou <- Q[observed_idx, unobserved_idx, drop = FALSE]
    Q_uo <- Matrix::t(Q_ou)
    Q_uu <- Q[unobserved_idx, unobserved_idx, drop = FALSE]

    # Compute conditional mean and covariance
    #mu_cond <- -Q_uu_inv %*% Q_uo %*% x_obs
    mu_cond <- -1 * Matrix::solve(Q_uu, Q_uo %*% x_obs)

    out = list( mean = as.vector(mu_cond),
                Q_uu = Q_uu,
                unobserved_idx = unobserved_idx )
    return(out)
  }

  res <- predict_conditional_gmrf(Q, observed_idx, x_obs )
  if( what == "predict" ){
    return(res)
  }else{
    y = rmvnorm_prec( n = n_sims, mu = res$mean, prec = res$Q_uu )
    return(y)
  }
}


#' @title Project tinyVAST to future times (EXPERIMENTAL)
#'
#' @description
#' Projects a fitted model forward in time.
#'
#' @inheritParams predict.tinyVAST
#' @param object fitted model from \code{tinyVAST(.)}
#' @param extra_times a vector of extra times, matching values in \code{newdata}
#' @param newdata data frame including new values for \code{time_variable}
#' @param future_var logical indicating whether to simulate future process errors
#'        from GMRFs, or just compute the predictive mean
#' @param past_var logical indicating whether to re-simulate past process errors
#'        from predictive distribution of random effects, thus changing the boundary
#'        condition of the forecast
#' @param parm_var logical indicating whether to re-sample fixed effects from their
#'        predictive distribution, thus changing the GMRF for future process errors
#' @param nsim number of samples
#'
#' @return
#' A vector of values corresponding to rows in \code{newdata}, or a matrix
#' with a column for each sample when \code{nsim > 1}
#'
#' @examples
#' # Convert to long-form
#' set.seed(123)
#' n_obs = 100
#' rho = 0.9
#' sigma_x = 0.2
#' sigma_y = 0.1
#' x = rnorm(n_obs, mean=0, sd = sigma_x)
#' for(i in 2:length(x)) x[i] = rho * x[i-1] + x[i]
#' y = x + rnorm( length(x), mean = 0, sd = sigma_y )
#' data = data.frame( "val" = y, "var" = "y", "time" = seq_along(y) )
#'
#' # Define AR2 time_term
#' time_term = "
#'   y -> y, 1, rho1
#'   y -> y, 2, rho2
#'   y <-> y, 0, sd
#' "
#'
#' # fit model
#' mytiny = tinyVAST(
#'   time_term = time_term,
#'   data = data,
#'   times = unique(data$t),
#'   variables = "y",
#'   formula = val ~ 1,
#'   control = tinyVASTcontrol( getJointPrecision = TRUE )
#' )
#'
#' # Deterministic projection
#' extra_times = length(x) + 1:100
#' n_sims = 10
#' newdata = data.frame( "time" = c(seq_along(x),extra_times), "var" = "y" )
#' Y = project(
#'   mytiny,
#'   newdata = newdata,
#'   extra_times = extra_times,
#'   future_var = FALSE
#' )
#' plot( x = seq_along(Y),
#'       y = Y,
#'       type = "l", lty = "solid", col = "black" )
#'
#' # Stochastic projection with future process errors
#' \dontrun{
#' extra_times = length(x) + 1:100
#' n_sims = 10
#' newdata = data.frame( "time" = c(seq_along(x),extra_times), "var" = "y" )
#' Y = project(
#'   mytiny,
#'   newdata = newdata,
#'   extra_times = extra_times,
#'   future_var = TRUE,
#'   past_var = TRUE,
#'   parm_var = TRUE,
#'   nsim = n_sims
#' )
#' matplot( x = row(Y),
#'          y = Y,
#'          type = "l", lty = "solid", col = "black" )
#' }
#'
#' @export
project <-
function( object,
          extra_times,
          newdata,
          what = "mu_g",
          future_var = TRUE,
          past_var = FALSE,
          parm_var = FALSE,
          nsim = 1 ){


  ##############
  # Step 1: Generate uncertainty from parm_var and past_var,
  #         and load into columns of parmat
  ##############

  if( isFALSE(parm_var) & isFALSE(past_var) ){
    parmat = object$obj$env$last.par.best %o% rep(1, nsim)
  }
  if( isTRUE(parm_var) & isFALSE(past_var) ){
    stop("option not available")
  }
  if( isFALSE(parm_var) & isTRUE(past_var) ){
    parmat = object$obj$env$last.par.best %o% rep(1, nsim)
    MC = object$obj$env$MC( keep=TRUE, n=nsim, antithetic=FALSE )
    parmat[object$obj$env$lrandom(),] = attr(MC, "samples")
  }
  if( isTRUE(parm_var) & isTRUE(past_var) ){
    if(is.null(object$sdrep$jointPrecision)) stop("Rerun with `getJointPrecision=TRUE`")
    parmat = rmvnorm_prec( mu = object$obj$env$last.par.best,
                           prec = object$sdrep$jointPrecision,
                           n = nsim )
  }

  ##############
  # Step 2: Augment objects
  ##############

  all_times = union( object$internal$times, extra_times )

  ##############
  # Step 3: Build object with padded bounds
  ##############

  new_control = object$internal$control
  new_control$run_model = FALSE
  new_control$opt_loops = 0
  new_control$newton_loops = 0
  new_control$getsd = FALSE
  new_control$calculate_deviance_explained = FALSE
  new_control$suppress_user_warnings = TRUE
  new_control$extra_reporting = TRUE

  new_inputs = tinyVAST(
    formula = object$formula,
    data = object$data,
    time_term = object$internal$time_term,
    space_term = object$internal$space_term,
    spacetime_term = object$internal$spacetime_term,
    family = object$internal$family,
    space_columns = object$internal$space_columns,
    spatial_domain = object$spatial_domain,
    time_column = object$internal$time_column,
    times = all_times,
    variable_column = object$internal$variable_column,
    variables = object$internal$variables,
    distribution_column = object$internal$distribution_column,
    delta_options = list( formula = object$internal$delta_formula,
                          space_term = object$internal$delta_space_term,
                          time_term = object$internal$delta_time_term,
                          spacetime_term = object$internal$delta_spacetime_term,
                          spatial_varying = object$internal$delta_spatial_varying ),
    spatial_varying = object$internal$spatial_varying,
    weights = object$internal$weights,
    control = new_control,
    development = object$internal$development
  )

  # REPORT uses expanded parameter lists, so omit both maps and random effects.
  # parList() on the fitted object restores fixed and shared mapped values.
  newobj = object
  newobj$tmb_inputs = new_inputs
  newobj$internal$times = all_times
  newobj$obj = MakeADFun( data = new_inputs$tmb_data,
                          parameters = new_inputs$tmb_par,
                          silent = TRUE,
                          DLL = "tinyVAST" )

  ##############
  # Step 4: Merge ParList and ParList1
  ##############

  augment_epsilon <-
  function( neweps_stc,
            eps_stc,
            #beta_z,
            #model,
            linpred ){

    if( prod(dim(eps_stc)) > 0 ){
    #if( length(beta_z) > 0 ){
      #
      #mats = dsem::make_matrices(
      #  beta_p = beta_z,
      #  model = model,
      #  variables = object$internal$variables,
      #  times = all_times
      #)
      #Q_kk = Matrix::t(mats$IminusP_kk) %*% solve(Matrix::t(mats$G_kk) %*% mats$G_kk) %*% mats$IminusP_kk
      if(linpred==1){
        IminusRho_hh = Matrix::Diagonal(n=nrow(newrep$Rho_hh)) - newrep$Rho_hh
        Q_kk = Matrix::t(IminusRho_hh) %*% Matrix::solve(Matrix::t(newrep$Gamma_hh) %*% newrep$Gamma_hh) %*% IminusRho_hh
      }else{
        IminusRho_hh = Matrix::Diagonal(n=nrow(newrep$Rho_hh)) - newrep$Rho2_hh
        Q_kk = Matrix::t(IminusRho_hh) %*% Matrix::solve(Matrix::t(newrep$Gamma2_hh) %*% newrep$Gamma2_hh) %*% IminusRho_hh
      }

      #
      grid = expand.grid( t = all_times,
                          c = object$internal$variables )
      grid$num = seq_len(nrow(grid))
      observed_idx = subset( grid, t %in% object$internal$times )$num
      unobserved_idx = setdiff( grid$num, observed_idx )

      # Precision is kronecker( Q_kk, Q_ss ) and past times include all sites, so the
      # conditional mean only needs Q_kk and the conditional precision is kronecker( Q_uu, Q_ss )
      Q_uu = Q_kk[unobserved_idx, unobserved_idx, drop = FALSE]
      Q_uo = Q_kk[unobserved_idx, observed_idx, drop = FALSE]
      eps_so = matrix( eps_stc, nrow = dim(eps_stc)[1] )
      simeps_h = -1 * as.matrix( eps_so %*% Matrix::t(Matrix::solve(Q_uu, Q_uo)) )
      if( isTRUE(future_var) ){
        # Columns have precision Q_ss, and are then correlated to have precision Q_uu among columns
        z_su = rmvnorm_prec( prec = Q_ss, n = length(unobserved_idx) )
        simeps_h = simeps_h + t( backsolve(chol(as.matrix(Q_uu)), t(z_su)) )
      }

      # Compile
      #missing_indices = as.matrix(subset( grid, t %in% extra_times )[,1:3])
      #neweps_stc[missing_indices] = simeps_stc[,1]
      tset = match( extra_times, all_times )
      neweps_stc[,tset,] = simeps_h
      #observed_indices = as.matrix(subset( grid, t %in% object$internal$times )[,1:3])
      #neweps_stc[observed_indices] = eps_stc[observed_indices]
      tset = match( object$internal$times, all_times )
      neweps_stc[,tset,] = eps_stc
    }
    return(neweps_stc)
  }
  augment_delta <-
  function( newdelta_tc,
            delta_tc,
            #nu_z,
            #model,
            linpred ){

    if( prod(dim(delta_tc)) > 0 ){
    #if( length(nu_z) > 0 ){
      #
      #mats = dsem::make_matrices(
      #  beta_p = nu_z,
      #  model = model,
      #  variables = object$internal$variables,
      #  times = all_times
      #)
      #Q_kk = Matrix::t(mats$IminusP_kk) %*% solve(Matrix::t(mats$G_kk) %*% mats$G_kk) %*% mats$IminusP_kk
      if(linpred==1){
        IminusRho_hh = Matrix::Diagonal(n=nrow(newrep$Rho_time_hh)) - newrep$Rho_time_hh
        Q_kk = Matrix::t(IminusRho_hh) %*% Matrix::solve(Matrix::t(newrep$Gamma_time_hh) %*% newrep$Gamma_time_hh) %*% IminusRho_hh
      }else{
        IminusRho_hh = Matrix::Diagonal(n=nrow(newrep$Rho2_time_hh)) - newrep$Rho2_time_hh
        Q_kk = Matrix::t(IminusRho_hh) %*% Matrix::solve(Matrix::t(newrep$Gamma2_time_hh) %*% newrep$Gamma2_time_hh) %*% IminusRho_hh
      }

      #
      grid = expand.grid( t = all_times,
                          c = object$internal$variables )
      grid$num = seq_len(prod(dim(newdelta_tc)))
      observed_idx = subset( grid, t %in% object$internal$times )$num

      #
      tmp = conditional_gmrf(
        Q = Q_kk,
        observed_idx = observed_idx,
        x_obs = as.vector( delta_tc ),
        n_sims = 1,
        what = ifelse(future_var, "simulate", "predict")
      )
      if( isTRUE(future_var) ){
        simdelta_k = tmp[,1]
      }else{
        simdelta_k = tmp$mean
      }

      # Compile
      #missing_indices = as.matrix(subset( grid, t %in% extra_times )[,1:2])
      #newdelta_tc[missing_indices] = simdelta_tc[,1]
      tset = match( extra_times, all_times )
      newdelta_tc[tset,] = simdelta_k
      #observed_indices = as.matrix(subset( grid, t %in% object$internal$times )[,1:2])
      #newdelta_tc[observed_indices] = delta_tc[observed_indices]
      tset = match( object$internal$times, all_times )
      newdelta_tc[tset,] = delta_tc
    }
    return(newdelta_tc)
  }

  ##############
  # Step 5: Re-build model
  ##############

  # Build prediction object once, and then REPORT for each sample
  predobj = MakeADFun( data = add_predictions( object = newobj, newdata = newdata ),
                       parameters = newobj$tmb_inputs$tmb_par,
                       type = "Fun",
                       silent = TRUE,
                       DLL = "tinyVAST" )

  ##############
  # Step 6: simulate samples
  ##############

  # Use TMB's parameter order when flattening lists for REPORT.
  par_template = newobj$obj$env$parList()
  same_vars = setdiff( names(par_template), c("epsilon_stc","epsilon2_stc","delta_tc","delta2_tc") )
  pred = matrix( NA, nrow = nrow(newdata), ncol = nsim )
  for( i in seq_len(nsim) ){
    parvec = parmat[,i]
    parlist = object$obj$env$parList( par = parvec )

    # Copy parameters with unchanged dimensions before building future precisions.
    # Process arrays retain their padded dimensions until augmented below.
    new_parlist = par_template
    new_parlist[same_vars] = parlist[same_vars]
    newrep = newobj$obj$report( unlist(new_parlist, use.names = FALSE) )

    #
    # random effects are scaled by tau, so their precision is Q_ss * tau^2
    Q_ss = newrep$Q_ss * exp(2 * newrep$log_tau)

    # Replace epsilon
    new_parlist$epsilon_stc = augment_epsilon(
      #beta_z = parlist$beta_z,
      eps_stc = parlist$epsilon_stc,
      neweps_stc = new_parlist$epsilon_stc,
      #model = object$internal$spacetime_term_ram$output$model,
      linpred = 1
    )
    new_parlist$epsilon2_stc = augment_epsilon(
      #beta_z = parlist$beta2_z,
      eps_stc = parlist$epsilon2_stc,
      neweps_stc = new_parlist$epsilon2_stc,
      #model = object$internal$delta_spacetime_term_ram$output$model,
      linpred = 2
    )

    # Replace delta
    new_parlist$delta_tc = augment_delta(
      #nu_z = parlist$nu_z,
      delta_tc = parlist$delta_tc,
      newdelta_tc = new_parlist$delta_tc,
      #model = object$internal$time_term_ram$output$model,
      linpred = 1
    )
    new_parlist$delta2_tc = augment_delta(
      #nu_z = parlist$nu2_z,
      delta_tc = parlist$delta2_tc,
      newdelta_tc = new_parlist$delta2_tc,
      #model = object$internal$delta_time_term_ram$output$model,
      linpred = 2
    )

    # Include every augmented effect, including effects fitted as fixed parameters.
    pred[,i] = predobj$report( unlist(new_parlist, use.names = FALSE) )[[what]]
  }
  if( nsim == 1 ) pred = pred[,1]
  return(pred)
}
