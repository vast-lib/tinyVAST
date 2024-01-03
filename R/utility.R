#' @title Multivariate Normal Random Deviates using Sparse Precision
#'
#' @description This function provides a random number generator for
#'              the multivariate normal distribution with mean equal
#'              to `mean` and sparse precision matrix `Q`.
#'
#' @param n number of observations.
#' @param mean mean vector.
#' @param Q sparse precision (inverse-covariance) matrix.
#'
#' @return a matrix with dimension \code{length(mean)} by
#'         \code{n}, containing realized draws from the specified
#'         mean and precision
#'
#' @export
rmvnorm_prec <-
function( n,
          Q,
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


