

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
edges = unique(edge)

# Compute edge lengths
edge_lengths <- sqrt(
  (loc[edges[,1], 1] - loc[edges[,2], 1])^2 +
  (loc[edges[,1], 2] - loc[edges[,2], 2])^2
)

min(edge_lengths)
