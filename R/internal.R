.onAttach <- function(libname, pkgname) {
}

ivector_minus_one <- function( ivector ){
  if( any(is.na(ivector)) ) stop("Check ivector for NAs")
  if( length(ivector) > 0 ){
    ivector = as.integer( ivector - 1 )
  }
  return(ivector)
}

# Modified from sdmTMB::make_anisotropy_spde, originally from geostatistical_delta-GLMM
make_anisotropy_spde <-
function( inla_mesh,
          covariates ){

  spde = fm_fem( inla_mesh )
  if( missing(covariates) ){
    loc = inla_mesh$loc[,1:2]
  }else{
    if(nrow(inla_mesh$loc) != nrow(covariates)) stop("Check `covariates` in `make_anisotropy_spde`")
    loc = cbind( inla_mesh$loc[,1:2], covariates )
  }
  TV <- inla_mesh$graph$tv
  V0 <- loc[ TV[,1], ]
  V1 <- loc[ TV[,2], ]
  V2 <- loc[ TV[,3], ]
  E0 <- V2 - V1
  E1 <- V0 - V2
  E2 <- V1 - V0
  TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
  Tri_Area <- rep(NA, nrow(E0))
  for (i in seq_len(length(Tri_Area))){
    Tri_Area[i] <- TmpFn( E0[i,1:2], E1[i,1:2] ) / 2
  }
  ret <- list( n_s = inla_mesh$n,
               n_tri = nrow(TV),
               Tri_Area = Tri_Area,
               E0 = E0,
               E1 = E1,
               E2 = E2,
               TV = TV - 1,
               G0 = spde$c0,
               G0_inv = as(Matrix::diag(1/Matrix::diag(spde$c0)),"TsparseMatrix") )
  return(ret)
}

# Recoded in R from https://github.com/pfmc-assessments/geostatistical_delta-GLMM/blob/master/inst/executables/geo_index_v4b.cpp#L18
# Intended to show logic of geometric anisotropy
# if H = diag(2), G is expected to match `fmesher::fm_fem(mesh)$g1`
make_stiffness <-
function( mesh,
          loc, # Can swap in different values
          H ){

  # local objects to simplify code
  if(missing(loc)) loc = mesh$loc[,1:2]
  if(missing(H)) H = diag( rep(1,ncol(loc)) )
  tv = mesh$graph$tv
  n = mesh$n
  adjH = solve(H) * det(H)

  # Extract edge vectors
  v0 = loc[ tv[,1], ]
  v1 = loc[ tv[,2], ]
  v2 = loc[ tv[,3], ]
  e0 = v2 - v1
  e1 = v0 - v2
  e2 = v1 - v0

  # Loop through triangles
  G = Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(n, n))
  for (i in seq_len(nrow(e0)) ){
    # Get edges
    edgemat = rbind( e0[i,], e1[i,], e2[i,] )
    # Area from just the first two dimensions (x, y)
    triangle_area = abs(det(edgemat[1:2,1:2])) / 2

    # Local stiffness
    G_tri = edgemat %*% adjH %*% t(edgemat)

    # Assemble by summation
    G[ tv[i,], tv[i,] ] = G[ tv[i,], tv[i,] ] + G_tri / (4*triangle_area)
  }
  return(G)
}

# Modified from sdmTMB
check_tinyVAST_version <- function(version) {
  if( utils::packageVersion("tinyVAST") != version ){
    stop("Installed version of `tinyVAST` does not match the version used to
          run model.  Please re-install the same version, or re-run the model.")
  }
}

#rm_wsp <- function (x) {
#  # from brms:::rm_wsp()
#  # VIA sdmTMB smoothers.R
#  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
#  dim(out) <- dim(x)
#  out
#}

#all_terms <- function (x) {
#  # from brms:::all_terms()
#  # VIA sdmTMB smoothers.R
#  if (!length(x)) {
#    return(character(0))
#  }
#  if (!inherits(x, "terms")) {
#    x <- terms(stats::as.formula(x))
#  }
#  rm_wsp(attr(x, "term.labels"))
#}

#get_smooth_terms <- function(terms) {
#  # from brms:::all_terms()
#  # VIA sdmTMB smoothers.R
#  x1 <- grep("s\\(", terms)
#  x2 <- grep("t2\\(", terms)
#  c(x1, x2)
#}

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

#' @title Approximate spatial correlation
#'
#' @description Extract the approximated spatial correlation between one coordinate
#' and other coordinates using a sparse precision and SPDE mesh
#'
#' @param Q sparse precision matrix
#' @param mesh SPDE mesh
#' @param coord vector of length-2 with spatial coordinates for focal point
#' @param pred matrix with two columns and multiple rows, with location
#'        for points to predict correlation
#'
#' @return
#' A vector with length \code{nrow(pred)} giving the spatial correlation
get_cor <-
function( Q,
          mesh,
          coord,
          pred ){

  # Get projection matrices
  if( is(mesh,"fm_mesh_2d" ) ){
    # Projection from mesh to coord
    A_is = fm_evaluator( mesh, loc=matrix(coord,nrow=1) )$proj$A
    A_gs = fm_evaluator( mesh, loc=pred )$proj$A
  }else{
    stop("`get_corr` not implemented for supplied `spatial_domain`")
  }

  # Takahashi-Davis variance for sparsity pattern of Q
  Vsparse = sparseinv::Takahashi_Davis(Q)

  #
  var_s = diag(Vsparse)
  var_g = (A_gs %*% var_s)[,1]
  var_i = (A_is %*% var_s)[,1]

  # Get covariance from coord to pred as weighted average of vertices
  # FIXME:  INCLUDE COVARIANCE AMONG VERTICES
  cov_s = solve( Q, t(A_is) )[,1]
  cov_g = (A_gs %*% cov_s)[,1]

  # Calculate correlation
  cor_g = cov_g / sqrt(var_g) / sqrt(var_i)
  return(cor_g)
}

