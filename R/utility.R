#' @title Multivariate Normal Random Deviates using Sparse Precision
#'
#' @description This function provide a random number generator for
#'              the multivariate normal distribution with mean equal
#'              to mean and sparse precision matrix Q.
#'
#' @param n number of observations.
#' @param mean mean vector.
#' @param Q sparse precision matrix.
#'
#' @export
rmvnorm_prec <-
function( n,
          mean,
          Q ) {

  # Simulate values
  z0 = matrix( rnorm(length(mean) * n), ncol=n)

  # Q = t(P) * L * t(L) * P
  L = Cholesky(Q, super=TRUE)

  # Calcualte t(P) * solve(t(L)) * z0 in two steps
  z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = solve(L, z, system = "Pt") # z = Pt    * z
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
  L_tf_rot = (Eigen$vectors%*%diag(sqrt(Eigen$values)))[,1:ncol(L_tf),drop=FALSE]

  # My new factors
  H = pseudoinverse(L_tf_rot) %*% L_tf
  x_sf = t(H %*% t(x_sf))

  # Get all loadings matrices to be increasing or decreasing
  order = match.arg(order)
  if( !is.null(order) ){
    for( f in 1:ncol(L_tf) ){
      Lm = lm( L_tf_rot[,f] ~ 1 + I(1:nrow(L_tf)) )
      Sign = sign(Lm$coef[2]) * ifelse(order=="decreasing", -1, 1)
      L_tf_rot[,f] = L_tf_rot[,f] * Sign
      x_sf[,f] = x_sf[,f] * Sign
    }
  }

  # return
  out = list( "L_tf"=L_tf_rot, "x_sf"=x_sf, "H"=H)
  return(out)
}
