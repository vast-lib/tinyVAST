
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
#'
#' @return
#' A matrix with \code{n_sims} columns and a row for every row of \code{Q} not in
#' \code{observed_idx}, with simulations for those rows
#'
#' @export
simulate_conditional_gmrf <-
function( Q,
          observed_idx,
          x_obs,
          n_sims = 1 ){

  # Required libraries
  #library(Matrix)

  #
  if( !all(observed_idx %in% seq_len(nrow(Q))) ){
    stop("Check `observed_idx` in `simulate_conditional_gmrf")
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
  y = rmvnorm_prec( n=n_sims, mu = res$mean, prec = res$Q_uu )
  return(y)
}
