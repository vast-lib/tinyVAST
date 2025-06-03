
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
    Q_oo <- Q[observed_idx, observed_idx, drop = FALSE]
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

  # Unconditional simulation from precision matrix
  rmvnorm_prec <-
  function( mu, # estimated fixed and random effects
            prec, # estimated joint precision
            n.sims) {

    # Simulate values
    z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    # Q = t(P) * L * t(L) * P
    L = Matrix::Cholesky(prec, super=TRUE)
    # Calcualte t(P) * solve(t(L)) * z0 in two steps
    z = Matrix::solve(L, z0, system = "Lt") # z = Lt^-1 * z
    z = Matrix::solve(L, z, system = "Pt") # z = Pt    * z
    return(mu + as.matrix(z))
  }

  res <- predict_conditional_gmrf(Q, observed_idx, x_obs )
  if( what == "predict" ){
    return(res)
  }else{
    y = rmvnorm_prec( n=n_sims, mu = res$mean, prec = res$Q_uu )
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
#'
#' @return
#' A vector of values corresponding to rows in \code{newdata}
#'
#' @export
project <-
function( object,
          extra_times,
          newdata,
          what = "mu_g",
          future_var = TRUE ){


  ##############
  # Step 1: Generate uncertainty in historical period
  ##############

  # SKIP STEP

  ##############
  # Step 2: Augment objects
  ##############

  all_times = c( object$internal$times, extra_times )

  ##############
  # Step 3: Build object with padded bounds
  ##############

  new_control = object$internal$control
  new_control$run_model = FALSE
  new_control$suppress_user_warnings = TRUE

  newobj = tinyVAST(
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
    spatial_varying = object$internal$spatially_varying,
    weights = object$internal$weights,
    control = new_control
  )

  ##############
  # Step 4: Merge ParList and ParList1
  ##############

  augment_epsilon <-
  function( neweps_stc,
            eps_stc,
            beta_z,
            model ){

    if( length(beta_z) > 0 ){
      #
      mats = dsem::make_matrices(
        beta_p = beta_z,
        model = model,
        variables = object$internal$variables,
        times = all_times
      )
      Q_kk = Matrix::t(mats$IminusP_kk) %*% solve(Matrix::t(mats$G_kk) %*% mats$G_kk) %*% mats$IminusP_kk
      Q_hh = Matrix::kronecker( Q_kk, Q_ss )

      #
      grid = expand.grid( s = seq_len(dim(neweps_stc)[1]),
                          t = all_times,
                          c = object$internal$variables )
      grid$num = seq_len(prod(dim(neweps_stc)))
      observed_idx = subset( grid, t %in% object$internal$times )$num

      #
      tmp = conditional_gmrf(
        Q = Q_hh,
        observed_idx = observed_idx,
        x_obs = as.vector( eps_stc ),
        n_sims = 1,
        what = ifelse(future_var, "simulate", "predict")
      )
      if( future_var == "simulate" ){
        simeps_h = tmp[,1]
      }else{
        simeps_h = tmp$mean
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
            nu_z,
            model ){

    if( length(nu_z) > 0 ){
      #
      mats = dsem::make_matrices(
        beta_p = nu_z,
        model = model,
        variables = object$internal$variables,
        times = all_times
      )
      Q_kk = Matrix::t(mats$IminusP_kk) %*% solve(Matrix::t(mats$G_kk) %*% mats$G_kk) %*% mats$IminusP_kk

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
      if( future_var == "simulate" ){
        simdelta_k = tmp[,1]
      }else{
        simeps_h = tmp$mean
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

  #
  parlist = object$internal$parlist
  new_parlist = newobj$tmb_par
  Q_ss = object$rep$Q_ss

  # Replace epsilon
  new_parlist$epsilon_stc = augment_epsilon(
    beta_z = parlist$beta_z,
    eps_stc = parlist$epsilon_stc,
    neweps_stc = new_parlist$epsilon_stc,
    model = object$internal$spacetime_term_ram$output$model
  )
  new_parlist$epsilon2_stc = augment_epsilon(
    beta_z = parlist$beta2_z,
    eps_stc = parlist$epsilon2_stc,
    neweps_stc = new_parlist$epsilon2_stc,
    model = object$internal$delta_spacetime_term_ram$output$model
  )

  # Replace delta
  new_parlist$delta_tc = augment_delta(
    nu_z = parlist$nu_z,
    delta_tc = parlist$delta_tc,
    newdelta_tc = new_parlist$delta_tc,
    model = object$internal$time_term_ram$output$model
  )
  new_parlist$delta2_tc = augment_delta(
    nu_z = parlist$nu2_z,
    delta_tc = parlist$delta2_tc,
    newdelta_tc = new_parlist$delta2_tc,
    model = object$internal$delta_time_term_ram$output$model
  )

  # Replace other variables that are not changed
  same_vars = setdiff( names(new_parlist), c("epsilon_stc","epsilon2_stc","delta_tc","delta2_tc") )
  new_parlist[same_vars] = parlist[same_vars]

  ##############
  # Step 5: Re-build model
  ##############

  new_control$run_model = TRUE
  new_control$tmb_par = new_parlist
  new_control$nlminb_loops = 0
  new_control$newton_loops = 0
  new_control$getsd = FALSE
  new_control$calculate_deviance_explained = FALSE

  newobj = tinyVAST(
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
    spatial_varying = object$internal$spatially_varying,
    weights = object$internal$weights,
    control = new_control
  )

  ##############
  # Step 6: simulate samples
  ##############

  pred = predict(
    object = newobj,
    newdata = newdata,
    what = what
  )
  return(pred)
}
