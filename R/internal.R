.onAttach <- function(libname, pkgname) {
}

ivector_minus_one <- function( ivector, varname = "" ){
  if( any(is.na(ivector)) ){
    stop("Check ivector for NAs: ", varname )
  }
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

  # Extract vertices
  TV <- inla_mesh$graph$tv
  V0 <- loc[ TV[,1], ]
  V1 <- loc[ TV[,2], ]
  V2 <- loc[ TV[,3], ]

  # Fill NAs as needed
  replace_NA = function(m) ifelse(is.na(m),0,m)
  E0 <- replace_NA(V2 - V1)
  E1 <- replace_NA(V0 - V2)
  E2 <- replace_NA(V1 - V0)

  # Get triangle areas
  #TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
  #Tri_Area <- rep(NA, nrow(E0))
  #for (i in seq_len(length(Tri_Area))){
  #  Tri_Area[i] <- TmpFn( E0[i,1:2], E1[i,1:2] ) / 2
  #}
  Tri_Area <- abs(E0[,1] * E1[,2] - E0[,2] * E1[,1]) / 2
  
  # Return it all
  G0_inv = as(Matrix::Diagonal(n = inla_mesh$n, x = 1/Matrix::diag(spde$c0)),"TsparseMatrix")
  ret <- list( n_s = inla_mesh$n,
               n_tri = nrow(TV),
               Tri_Area = Tri_Area,
               E0 = E0,
               E1 = E1,
               E2 = E2,
               TV = TV - 1,
               G0 = spde$c0,
               G0_inv = G0_inv )
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

# Compute minimum edge length for SPDE mesh
edge_lengths <- 
function( mesh ){
  # mesh is an fm_mesh_2d object
  loc <- mesh$loc           # vertex coordinates (n x 2)
  tv  <- mesh$graph$tv      # triangle vertices (m x 3)
  
  # Build unique edges from triangles
  edges = NULL
  for(i in 1:nrow(tv)){
    x = tv[i,]
    edges = rbind(
      edges,
      sort(c(x[1], x[2])),
      sort(c(x[2], x[3])),
      sort(c(x[3], x[1]))
    )
  }
  edges = unique(edges)
  
  # Compute edge lengths
  edge_lengths <- sqrt(
    (loc[edges[,1], 1] - loc[edges[,2], 1])^2 +
    (loc[edges[,1], 2] - loc[edges[,2], 2])^2
  )
  
  return(edge_lengths)
}

# Check for areal model matching assumptions
is_areal_sf <- function(x) {
  if (inherits(x, "sf")) {
    geom <- st_geometry(x)
  } else if (inherits(x, "sfc")) {
    geom <- x
  } else {
    return(FALSE)
  }
  all(st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON"))
}

# Make list of data containing SCALAR or VECTOR elements representing a NNGP
make_nngp_data =
function( coords, nn, what = c("full","empty") ){
  what = match.arg(what)
  if( what == "empty" ){
    return(list(
      #n = integer(0),
      nn_index_flat = integer(0),
      #nn_start = integer(0),
      nn_len = integer(0),
      #nn_mat_start = integer(0),
      dist_to_nn_flat = numeric(0),
      dist_within_nn_flat = numeric(0),
      gp_order = integer(0)
    ))
  }

  precompute_nngp_structure = function(coords, nn) {
    n = nrow(coords)
    nn_index = vector("list", n)
    dist_to_nn = vector("list", n)
    dist_within_nn = vector("list", n)

    for (i in seq_len(n)) {
      nn_ids = nn[i, -1]
      nn_ids = nn_ids[!is.na(nn_ids)]

      if (length(nn_ids) == 0L) {
        nn_index[[i]] = integer(0)
        dist_to_nn[[i]] = numeric(0)
        dist_within_nn[[i]] = matrix(0, 0, 0)
        next
      }

      nn_index[[i]] = nn_ids
      dist_to_nn[[i]] = sqrt(rowSums((coords[nn_ids, , drop = FALSE] - coords[rep(i, length(nn_ids)), , drop = FALSE])^2))
      dist_within_nn[[i]] = as.matrix(dist(coords[nn_ids, , drop = FALSE]))
    }

    list(
      n = n,
      nn_index = nn_index,
      dist_to_nn = dist_to_nn,
      dist_within_nn = dist_within_nn
    )
  }

  # Re-order
  gp_order = order_maxmin(coords)
  ordered_nn = find_ordered_nn( coords[gp_order,], m = 4L)
  ordered_structure = precompute_nngp_structure( coords[gp_order,], nn = ordered_nn)

  ##############
  # Flatten
  ##############

  #n = length(ordered_structure$nn_index)

  #nn_index_flat = integer(0)
  #dist_to_nn_flat = numeric(0)
  #dist_within_nn_flat = numeric(0)

  #nn_start = integer(n)
  #nn_len   = integer(n)
  #nn_mat_start = integer(n)

  #pos_vec = 0   # position for vector data (nn_index, dist_to_nn)
  #pos_mat = 0   # position for matrix data (dist_within_nn)

  #for (i in seq_len(n)) {
  #  nn_ids = ordered_structure$nn_index[[i]]
  #  k = length(nn_ids)
  #
  #  # Record starts and lengths
  #  #nn_start[i] = pos_vec
  #  nn_len[i] = k
  #  #nn_mat_start[i] = pos_mat
  #
  #  if (k > 0) {
  #    # Flatten neighbor indices (convert to 0-based for TMB)
  #    nn_index_flat = c(nn_index_flat, nn_ids - 1L)
  #
  #    # Distances to neighbors
  #    dist_to_nn_flat = c(
  #      dist_to_nn_flat,
  #      ordered_structure$dist_to_nn[[i]]
  #    )
  #
  #    # Flatten k x k matrix (column-major, consistent with R)
  #    dist_within_nn_flat = c(
  #      dist_within_nn_flat,
  #      as.vector(ordered_structure$dist_within_nn[[i]])
  #    )
  #
  #    # Advance positions
  #    #pos_vec = pos_vec + k
  #    #pos_mat = pos_mat + k * k
  #  }
  #}

  # Equivalent
  nn_len = sapply( ordered_structure$nn_index, FUN = length )
  nn_index_flat = unlist(ordered_structure$nn_index) - 1
  dist_to_nn_flat = unlist(ordered_structure$dist_to_nn)
  dist_within_nn_flat = unlist(lapply(ordered_structure$dist_within_nn, FUN=as.vector))

  # Return
  data = list(
    #n = n,
    nn_index_flat = nn_index_flat,
    #nn_start = nn_start,
    nn_len = nn_len,
    #nn_mat_start = nn_mat_start,
    dist_to_nn_flat = dist_to_nn_flat,
    dist_within_nn_flat = dist_within_nn_flat,
    gp_order = gp_order - 1L
  )
  return(data)
}

