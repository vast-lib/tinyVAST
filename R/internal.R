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
function( inla_mesh ){

  spde = fm_fem( inla_mesh )
  Dset <- 1:2
  TV <- inla_mesh$graph$tv
  V0 <- inla_mesh$loc[TV[, 1], Dset]
  V1 <- inla_mesh$loc[TV[, 2], Dset]
  V2 <- inla_mesh$loc[TV[, 3], Dset]
  E0 <- V2 - V1
  E1 <- V0 - V2
  E2 <- V1 - V0
  TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
  Tri_Area <- rep(NA, nrow(E0))
  for (i in seq_len(length(Tri_Area))){
    Tri_Area[i] <- TmpFn(E0[i,], E1[i,])/2
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
          H = diag(2) ){

  # local objects to simplify code
  dims = 1:2
  tv = mesh$graph$tv
  n = mesh$n
  adjH = solve(H) * det(H)

  # Extract edge vectors
  v0 = mesh$loc[ tv[,1], dims]
  v1 = mesh$loc[ tv[,2], dims]
  v2 = mesh$loc[ tv[,3], dims]
  e0 = v2 - v1
  e1 = v0 - v2
  e2 = v1 - v0

  # Loop through triangles
  G = Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(n, n))
  for (i in seq_along(nrow(e0)) ){
    triangle_area = abs(det(edgemat[1:2,])) / 2
    edgemat = rbind( e0[i,], e1[i,], e2[i,] )

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


