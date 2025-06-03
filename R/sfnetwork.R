
#' @title Construct projection matrix for stream network
#'
#' @description Make sparse matrix to project from stream-network nodes
#'   to user-supplied points
#'
#' @param stream \pkg{sfnetworks} object representing stream network
#' @param loc \pkg{sf} object representing points to which are being projected
#' @param tolerance error-check tolerance
#'
#' @importFrom units drop_units
#' @importFrom sf st_as_sf st_nearest_feature st_distance st_geometry st_length
#'   st_crs
#' @importFrom sfnetworks activate
#' @importFrom Matrix Diagonal
#'
#' @return the sparse interpolation matrix, with rows for each row of \code{data}
#' supplied during fitting and columns for each spatial random effect.
#'
#' @export
sfnetwork_evaluator <-
function( stream,
          loc,
          tolerance = 0.01 ){

  # error checks
  #require(sfnetworks)
  #require(Matrix)
  #require(units)
  #require(sf)
  if(!inherits(loc,"sfc")){
    #stop("Must provide `loc` with class `sf`")
    loc = st_as_sf( as.data.frame(loc), coords=colnames(loc), crs=st_crs(stream) )
  }
  #if(!inherits(stream,"sfnetwork")) stop("Must provide `stream` with class `sfnetwork`")

  # Edge and vertex copies
  edges = st_as_sf(activate(stream,"edges"))
  nodes = st_as_sf(activate(stream,"nodes"))

  # Get nearest edge
  closest_edge = st_nearest_feature( x=loc, edges )
  #closest_edge_dist = st_distance(sf_samples, edges[closest_edge,], by_element=TRUE)

  # Get distance from associated vertices
  dist_from = drop_units(st_distance(loc, nodes[edges$from[closest_edge],], by_element=TRUE))
  dist_to = drop_units(st_distance(loc, nodes[edges$to[closest_edge],], by_element=TRUE))

  # Make sparse matrix
  A_is = Matrix::sparseMatrix(
    dims = c( length(st_geometry(loc)), length(st_geometry(nodes)) ),
    i = rep(seq_along(st_geometry(loc)),2),
    j = c( edges$from[closest_edge], edges$to[closest_edge] ),
    x = 1 - c( dist_from/(dist_from+dist_to), dist_to/(dist_from+dist_to) )
  )

  # checks
  if( any(abs(Matrix::rowSums(A_is)-1) > tolerance ) ){
    stop("Check rowsums")
  }

  return(A_is)
}

#' @title Make mesh for stream network
#'
#' @description make an object representing spatial information required
#' to specify a stream-network spatial domain, similar in usage to
#' \code{link[fmesher]{fm_mesh_2d}} for a 2-dimensional continuous domain
#'
#' @inheritParams sfnetwork_evaluator
#'
#' @return
#' An object (list) of class `sfnetwork_mesh`. Elements include:
#' \describe{
#' \item{N}{The number of random effects used to represent the network}
#' \item{table}{a table containing a description of parent nodes (from),
#'              childen nodes (to), and the distance separating them}
#' \item{stream}{copy of the stream network object passed as argument}
#' }
#'
#' @export
sfnetwork_mesh <-
function( stream ){

  #require(sfnetworks)
  #require(units)
  #require(Matrix)

  # from, to, dist
  edges = st_as_sf(activate(stream,"edges"))
  nodes = st_as_sf(activate(stream,"nodes"))
  N = nrow(nodes)

  # Old table
  table = data.frame( from = edges$from,
                     to = edges$to,
                     dist = drop_units(st_length(edges)) )

  # Error check
  number_of_parents = table( factor(table$to, levels=seq_len(N)) )
  if( max(number_of_parents) > 1 ){
    stop("Stream network has multiple parents for a node, i.e., isn't ordered upstream as assumed")
  }

  # Sparse matrix form ... doesn't work because sparseMatrix sums distances from multiple cells in M2 prior to nonlinear transformation
  # Have to calculate sparse matricse in CPP using the sparse distance matrix
    # Q = I + Q1 + t(Q1) + Q2 + Q3
    # Q1 = exp(-theta*M1) / (1 - exp(-2*theta*M1)) *
    # Q2 = exp(-2*theta*M2) / (1 - exp(-2*theta*M2))
    # Q3 = exp(-2*theta*M3) / (1 - exp(-2*theta*M3))
  #M1 = sparseMatrix( i=table$from, j=table$to, x=table$dist, dims=c(N,N) )
  #M2 = sparseMatrix( i=table$from, j=table$from, x=table$dist, dims=c(N,N) )
  #M3 = sparseMatrix( i=table$to, j=table$to, x=table$dist, dims=c(N,N) )

  # Sparse distance matrix
  Dist_ss = sparseMatrix( i = table$to,
                          j = table$from,
                          x = table$dist,
                          dims = c(N,N) )

  #
  out = structure( list(
    n = N,
    table = table,
    stream = stream,
    Dist_ss = Dist_ss
  ), class="sfnetwork_mesh" )
  return(out)
}

#' @title Simulate GMRF for stream network
#'
#' @description
#' Simulate values from a GMRF using a tail-down (flow-unconnected) exponential
#' model on a stream network
#'
#' @param sfnetwork_mesh Output from \code{\link{sfnetwork_mesh}}
#' @param theta Decorrelation rate
#' @param n number of simulated GMRFs
#' @param what Whether to return the simulated GMRF or its precision matrix
#'
#' @return a matrix of simulated values for a Gaussian Markov random field
#' arising from a stream-network spatial domain, with row for each spatial random
#' effect and \code{n} columns, using the sparse precision matrix
#' defined in Charsley et al. (2023)
#'
#' @references
#' Charsley, A. R., Gruss, A., Thorson, J. T., Rudd, M. B., Crow, S. K.,
#' David, B., Williams, E. K., & Hoyle, S. D. (2023). Catchment-scale
#' stream network spatio-temporal models, applied to the freshwater stages
#' of a diadromous fish species, longfin eel (Anguilla dieffenbachii).
#' Fisheries Research, 259, 106583. \doi{10.1016/j.fishres.2022.106583}
#'
#' @export
simulate_sfnetwork <-
function( sfnetwork_mesh,
          theta,
          n = 1,
          what = c("samples","Q") ){

  what = match.arg(what)
  N = sfnetwork_mesh$n

  # OLD WAY
  #table = sfnetwork_mesh$table
    # Q1 = exp(-theta*M1) / (1 - exp(-2*theta*M1))
    # Q2 = exp(-2*theta*M2) / (1 - exp(-2*theta*M2))
    # Q3 = exp(-2*theta*M3) / (1 - exp(-2*theta*M3))
  #Q1 = sparseMatrix( i=table$from, j=table$to, x=-exp(-theta*table$dist)/(1-exp(-2*theta*table$dist)), dims=c(N,N) )
  #Q2 = sparseMatrix( i=table$from, j=table$from, x=exp(-2*theta*table$dist)/(1-exp(-2*theta*table$dist)), dims=c(N,N) )
  #Q3 = sparseMatrix( i=table$to, j=table$to, x=exp(-2*theta*table$dist)/(1-exp(-2*theta*table$dist)), dims=c(N,N) )
  #I = Diagonal( n=nrow(Q1) )
  #Q = I + Q1 + Matrix::t(Q1) + Q2 + Q3

  # NEW WAY
  triplet = mat2triplet( sfnetwork_mesh$Dist_ss )
  v1_s = exp(-theta * triplet$x) / (1 - exp(-2 * theta * triplet$x))
  v2_s = exp(-2 * theta * triplet$x) / (1 - exp(-2 * theta * triplet$x))
  # Make components
  P_ss = sparseMatrix( i = triplet$i,
                       j = triplet$j,
                       x = v1_s,
                       dims = c(N,N) )
  D_ss = Diagonal( N ) +
         sparseMatrix( i = triplet$i,
                       j = triplet$i,
                       x = v2_s,
                       dims = c(N,N) )
  # Assemble
  DminusP_ss = D_ss - P_ss
  Q = (t(DminusP_ss) %*% solve(D_ss) %*% DminusP_ss)


  #
  if(what=="samples") out = rmvnorm_prec( n=n, mean=rep(0,nrow(Q)), Q=Q )
  if(what=="Q") out = Q
  return(out)
}
