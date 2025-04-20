
# See: C:\Users\James.Thorson\Desktop\Work files\Collaborations\2025 -- Generalized graphical additive mixed model\mammal_traits

Q_network <-
function( log_theta,
          tree,
          method = 0 ){

  # Locals
  parent_s = tree$edge[,1]  # max(table(tree$edge[,1])) = 2
  child_s = tree$edge[,2]   # max(table(tree$edge[,2])) = 1
  dist_s = tree$edge.length
  n_node = max( tree$edge )
  n_edge = length(parent_s)
  theta = exp( log_theta );
  if(any(table(child_s))>1) stop("Not a tree or stream")

  if( is.na(method) ){
    # Drops edge lengths
    #graph = as.igraph.phylo( tree )
    graph = graph_from_edgelist( tree$edge )
    dist_ss = igraph::distances( graph, weights = graph )
  }

  # Assemble pieces
  if( method == 0 ){
    # NOTE: tmp_s is different than before
    v1_s = exp(-theta * dist_s) / (1 - exp(-2 * theta * dist_s))
    v2_s = exp(-2 * theta * dist_s) / (1 - exp(-2 * theta * dist_s))
    # Make components
    P_ss = sparseMatrix( i = child_s, j = parent_s, x = v1_s, dims = c(n_node,n_node) )
    D_ss = Diagonal( n_node ) + sparseMatrix( i = child_s, j = child_s, x = v2_s, dims = c(n_node,n_node) )
    # Assemble
    DminusP_ss = D_ss - P_ss
    Q = (t(DminusP_ss) %*% solve(D_ss) %*% DminusP_ss)
  }
  if( method == 3 ){
    # Previous method
    tmp_s = -exp(-theta * dist_s) / (1 - exp(-2 * theta * dist_s))
    tmp2_s = exp(-2 * theta * dist_s) / (1 - exp(-2 * theta * dist_s))
    # Start withdiagonal
    Q = sparseMatrix( i = seq_len(n_node), j = seq_len(n_node), x = 1, dims = c(n_node,n_node) )
    # Add extras
    for( s in seq_len(n_edge) ){
      tmp = -exp(-theta*dist_s[s]) / (1-exp(-2*theta*dist_s[s]));
      tmp2 = exp(-2*theta*dist_s[s]) / (1-exp(-2*theta*dist_s[s]));
      Q[parent_s[s], child_s[s]] = tmp;         # Q[parent_s[s], child_s[s]] +
      Q[child_s[s], parent_s[s]] = tmp;         # Q[child_s[s], parent_s[s]] +
      Q[parent_s[s], parent_s[s]] = Q[parent_s[s], parent_s[s]] + tmp2;
      Q[child_s[s], child_s[s]] = Q[child_s[s], child_s[s]] + tmp2;
    }
  }
  #i = c( parent_s, child_s, parent_s, child_s, seq_len(n_node) )
  #j = c( child_s, parent_s, parent_s, child_s, seq_len(n_node) )
  #x = c( tmp_s, tmp_s, tmp2_s, tmp2_s, rep(1,n_node) )
  #Q = sparseMatrix( i = i, j = j, x = x, dims = c(n_node, n_node), repr = "C" )
  return(Q)
}

library(ape)
library(Matrix)
library(igraph)

tree = rtree(2)
theta = 0.2
#tree$edge = tree$edge[,2:1] # data.frame( tree$edge[,2], tree$edge[,1] )

# Existing method
Q3 = Q_network( log_theta = log(theta),
               tree = tree,
               method = 3 )
Q0 = Q_network( log_theta = log(theta),
               tree = tree,
               method = 0 )
range( Q0 - Q3 )
