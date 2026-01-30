
library(igraph)
library(Matrix)

#time <- make_graph( ~ 1 - 2  )
#space <- make_star( 3, mode = "out" )
#V(space)$name = letters[seq_along(space)]
#V(space)$name = seq_along(space)

#
time <- make_graph( ~ 2001 -+ 2002  )
space <- make_star( 3, mode = "out" )
V(space)$name = letters[seq_along(space)]

# Unconnected graphs
empty_time = make_empty_graph( length(V(time)) )
  V(empty_time)$name = V(time)$name
empty_space = make_empty_graph( length(V(space)) )
  V(empty_space)$name = V(space)$name

############
# Bespoke direct ("tensor") product
############

# Identical to netUtils::graph_direct iff all edges are undirected (graph is symmetric)
# Generalizes for directed edges (asymetric graph)
direct_product <- 
function( g1, g2 ){
  P1 = as_adjacency_matrix(g1)
  P2 = as_adjacency_matrix(g2)
  P_joint = kronecker( P1, P2 )
  colnames(P_joint) = rownames(P_joint) = apply( expand.grid(rownames(P1),rownames(P2)), MARGIN = 1, FUN = paste, collapse = "_" )
  g_joint = graph_from_adjacency_matrix( P_joint )
  return(g_joint)
}

graph_both = direct_product( space, time )
plot(graph_both)

gh <- netUtils::graph_direct( time, space )    # tensor/direct product
plot(gh)


############
# Bespoke cartesian product
############

# Identical to netUtils::graph_direct iff all edges are undirected (graph is symmetric) and all names are numeric (has bug otherwise)
# Generalizes for directed edges (asymetric graph)
cartesian_product <- 
function( g1, g2 ){
  P1 = as_adjacency_matrix(g1)
  P2 = as_adjacency_matrix(g2)
  I1 = Diagonal(n=length(V(g1)))
    rownames(I1) = colnames(I1) = V(g1)$names
  I2 = Diagonal(n=length(V(g2)))
    rownames(I2) = colnames(I2) = V(g2)$names

  P_joint = kronecker( P1, I2 ) + kronecker( I1, P2 ) 
  colnames(P_joint) = rownames(P_joint) = apply( expand.grid(rownames(P1),rownames(P2)), MARGIN = 1, FUN = paste, collapse = "_" )
  g_joint = graph_from_adjacency_matrix( P_joint )
  return(g_joint)
}

graph_both = cartesian_product( space, time )
plot(graph_both)

gh = netUtils::graph_cartesian(time, space)
plot(gh)


############
# Superposition of cartesian products
############

superposition <- 
function( g1, g2 ){
  g_joint <- union(g1, g2)
  g_joint <- simplify(g_joint, edge.attr.comb = list(weight = "sum"))
  return( g_joint )
}

G1 = cartesian_product( space, empty_time )
G2 = cartesian_product( empty_space, time )

G = superposition( G1, G2 )

