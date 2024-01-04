
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
#' @inheritParams sfnetwork_evaluator
#'
#' @export
sfnetwork_mesh <-
function( stream ){

  #require(sfnetworks)
  #require(units)
  #require(Matrix)

  # from, to, dist
  edges = st_as_sf(activate(stream,"edges"))
  nodes = st_as_sf(activate(stream,"edges"))
  table = data.frame( from = edges$from,
                     to = edges$to,
                     dist = drop_units(st_length(edges)) )
  N = max( c(table$from, table$to) )

  # Sparse matrix form ... doesn't work because sparseMatrix sums multiple cells in M2 and M3 prior to nonlinear transformation
    # Q = I + Q1 + t(Q1) + Q2 + Q3
    # Q1 = exp(-theta*M1) / (1 - exp(-2*theta*M1)) *
    # Q2 = exp(-2*theta*M2) / (1 - exp(-2*theta*M2))
    # Q3 = exp(-2*theta*M3) / (1 - exp(-2*theta*M3))
  #M1 = sparseMatrix( i=table$from, j=table$to, x=table$dist, dims=c(N,N) )
  #M2 = sparseMatrix( i=table$from, j=table$from, x=table$dist, dims=c(N,N) )
  #M3 = sparseMatrix( i=table$to, j=table$to, x=table$dist, dims=c(N,N) )

  #
  out = structure( list(
    n = N,
    table = table,
    stream = stream
  ), class="sfnetwork_mesh" )
  return(out)
}

#' @title Simulate GMRF for stream network
#'
#' @param sfnetwork_mesh Output from \code{\link{sfnetwork_mesh}}
#' @param theta Decorrelation rate
#' @param n number of simulated GMRFs
#' @param what Whether to return the simulated GMRF or its precision matrix
#'
#' @export
simulate_sfnetwork <-
function( sfnetwork_mesh,
          theta,
          n = 1,
          what = c("samples","Q") ){

  what = match.arg(what)
  table = sfnetwork_mesh$table
  N = sfnetwork_mesh$n

    # Q1 = exp(-theta*M1) / (1 - exp(-2*theta*M1))
    # Q2 = exp(-2*theta*M2) / (1 - exp(-2*theta*M2))
    # Q3 = exp(-2*theta*M3) / (1 - exp(-2*theta*M3))
  Q1 = sparseMatrix( i=table$from, j=table$to, x=-exp(-theta*table$dist)/(1-exp(-2*theta*table$dist)), dims=c(N,N) )
  Q2 = sparseMatrix( i=table$from, j=table$from, x=exp(-2*theta*table$dist)/(1-exp(-2*theta*table$dist)), dims=c(N,N) )
  Q3 = sparseMatrix( i=table$to, j=table$to, x=exp(-2*theta*table$dist)/(1-exp(-2*theta*table$dist)), dims=c(N,N) )
  I = Diagonal( n=nrow(Q1) )
  Q = I + Q1 + Matrix::t(Q1) + Q2 + Q3

  #
  if(what=="samples") out = rmvnorm_prec( n=n, mean=rep(0,nrow(Q)), Q=Q )
  if(what=="Q") out = Q
  return(out)
}
