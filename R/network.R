check_neighbourhood_consistency <- function(#
  ig,
  vertex_names
) {
  #' For a set of vertices, check how consistent the neighbourhoods are
  #' reporting each v in union(N(i) for i in vertex set) how many of 
  #' vertex_set have this as a neighbour, and pos/neg edge ratio.

  all_neighbour <- Reduce(
    union,
    adjacent_vertices(ig, vertex_names)
  )
  vertex_set <- V(ig)[name %in% vertex_names]

  igp <- subgraph(ig, c(all_neighbour, vertex_set))

  # Convert to IDs in new subgraph
  all_neighbour <- Reduce(
    union,
    adjacent_vertices(igp, vertex_names)
  )
  vertex_set <- V(igp)[name %in% vertex_names]

  
  res_list <- list()
  for(n in all_neighbour) {
    if(n %in% vertex_set) {
      next
    }
    nn <- neighbors(igp, n)

    # Count how many of the vertices in vertex_set are neighbours of this
    # specific node. If the nodes are "consistent", should be |vertex_set|.
    # Can't think of a more elegant way to do this counting today, there
    # must be a nicer way though.
    count_vertices <- 0
    pos <- 0
    neg <- 0
    missing <- 0
    for(m in vertex_set) {
      if(m %in% nn) {
        count_vertices <- count_vertices + 1
        edge <- E(igp)[get_edge_ids(igp, c(m, n), directed = FALSE)]
        if(edge$corr > 0) {
          pos <- pos + 1
        } else {
          neg <- neg + 1
        }
      }
    }

    res_list[[V(igp)[n]$name]] <- list(
      name = V(igp)[n]$name,
      count = count_vertices,
      pos_ratio = pos / (pos + neg),
      type = V(igp)[n]$type
    )
  }

  tbl <- rbindlist(res_list) |>
   as_tibble() |>
   filter(type != "b") |>
   arrange(
    -count
   )
}