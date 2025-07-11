netcomi_graph_analysis <- function(
  hb_se,
  lb_se,
  fun_md,
  metab_md,
  metab_ttest,
  hb_ig,
  lb_ig
) {
  library(igraph)
  library(SpiecEasi)
  library(NetCoMi)
  library(tidyverse)
  library(tools)
  library(readxl)
  library(glue)

  LABELS_NETWORK <- c(
    "High Bioactive",
    "Low Bioactive"
  )
  metab_ttest <- metab_ttest |>
    mutate(metab_id = paste0(Project, `Feature ID`, `Ion mode`),
          enr_in = ifelse(`Individual Estimate` > 0, "High Bioactive", "Low Bioactive"))
  metab_md <- metab_md |>
    mutate(nwk_id =  paste0(project, feature_id, ionmode))
  BIO_COLORS = c("Low Bioactive" = "#619CFF", "High Bioactive" = "#00BA38")

  #### MAKE NETWORKS ####
  # These are preconstructed by SPIEC-EASI. To get to an association matrix, 
  # we want to have the estimated correlation, for any edges which are inferred
  # to be related. From the NetCoMi code, they derive this association matrix
  # as getOptNet * cov2cor - so edges not in the network are 0, others weighted
  # with their correlation. We do the same when converting to igraph, but just
  # use correlation as an edge property.
  association_matrices <- setNames(
    lapply(
      list(hb_se, lb_se),
      function(x) {
        se_net <- x
        cor <- stats::cov2cor(as.matrix(getOptCov(se_net)))
        ass_mat <- as.matrix(cor * SpiecEasi::getRefit(se_net))
        diag(ass_mat) <- 1
        colnames(ass_mat) <- colnames(se_net$est$data)
        rownames(ass_mat) <- colnames(se_net$est$data)
        return(ass_mat)
      }
    ),
    LABELS_NETWORK
  )

  networks <- netConstruct(
    association_matrices[[1]],
    association_matrices[[2]],
    dataType = "partialCorr",
    sparsMethod = "none",
    verbose = 0
  )


  #### TRANSFORM METADATA ####
  # Merge the metadata into a single table
  # read.delim rather than read_delim due to bug in vroom/bioc currently
  metadata <- fun_md %>% as_tibble() %>%
    rename(nwk_vert_id = 1) %>%
    mutate(md_type = "EM.PFAML0")
  if (!is.null(metab_md)) {
    md_b <- as_tibble(metab_md) %>%
      rename(nwk_vert_id = nwk_id) %>%
      mutate(md_type = "um_sig_fecal")
    metadata <- bind_rows(metadata, md_b)
  }
  # Function to get truncated names for long lineages
  trim_lineage <- function(lineage_str, delim = ";", unknown = "?") {
    parts <- unlist(strsplit(lineage_str, delim))
    if (length(parts) == 1) {
      return(lineage_str)
    }
    not_missing <- which(!(parts == "?"))
    last_good <- not_missing[length(not_missing)]
    nice_lineage <- parts[last_good:length(parts)]
    return(paste(nice_lineage, collapse = delim))
  }
  # Restrict to the vertices in the graph and match ordering
  metadata <- metadata %>%
    filter(nwk_vert_id %in% colnames(networks$assoMat1)) %>%
    mutate(short_label = map_chr(nwk_vert_id, trim_lineage)) %>%
    left_join(
      metab_ttest |>
        select(metab_id, enr_in),
      join_by(nwk_vert_id == metab_id)
    ) |>
    mutate(
      color_cat = ifelse(grepl("PFAM", md_type), "PFAM", enr_in)
    )

  ### PAIRED COMPARISON ####
  # Get a paired analysis object
  net_analysis <- netAnalyze(networks,
                            centrLCC = FALSE,
                            avDissIgnoreInf = TRUE,
                            sPathNorm = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            hubPar = c("betweenness"),
                            hubQuant = 0.95,
                            lnormFit = FALSE,
                            normDeg = TRUE,
                            normBetw = TRUE,
                            normClose = TRUE,
                            normEigen = TRUE
  )

  # Paired plot of the network
  pdf(width = 10, height = 5, file = "output/figures/fig_four_paired_network.pdf")
  plot(net_analysis,
      nodeFilter = "names",
      nodeFilterPar = base::union(net_analysis$lccNames1,
                                  net_analysis$lccNames2),
      sameLayout = TRUE,
      labels = setNames(metadata$short_label, metadata$nwk_vert_id),
      layoutGroup = "union",
      rmSingles = "inboth",
      labelScale = FALSE,
      # Node properties
      nodeSize = "degree",
      cexNodes = 1,
      nodeSizeSpread = 1.5,
      nodeColor = "feature",
      colorVec = list(c(BIO_COLORS['High Bioactive'], BIO_COLORS['Low Bioactive'], "#f58231"), 
                      c(BIO_COLORS['High Bioactive'], BIO_COLORS['Low Bioactive'], "#f58231")),
      nodeShape = c("circle", "square"),
      featVecCol = setNames(metadata$color_cat, metadata$nwk_vert_id),
      featVecShape = setNames(metadata$md_type, metadata$nwk_vert_id),
      nodeTransp = 20,
      borderWidth = 0.5,
      borderCol = "lightgrey",
      cexLabels = 0,
      # Hub settings
      highlightHubs = TRUE,
      hubTransp = 0,
      hubBorderWidth = 1,
      hubBorderCol  = "gray30",
      cexHubLabels = 0.3,
      # Edge properties
      edgeWidth = 2,
      edgeTranspLow = 50,
      edgeTranspHigh = 20,
      # Other properties
      cexTitle = 1,
      groupNames = LABELS_NETWORK,
      repulsion = 1,
      mar = c(1,1,1,1))
  
  dev.off()

  #### QUANTITATIVE COMPARISION ####
  comp_network <- netCompare(net_analysis,
                            permTest = FALSE,
                            verbose = TRUE,
                            seed = 123456)

  comp_sum <- summary(comp_network,
                      groupNames = c("Low Bioactive", "High Bioactive"),
                      showCentr = c("degree", "between", "closeness"),
                      numbNodes = 5)

  #### DIFFERENTIAL NETWORK ANALYSIS ####
  diff_net <- diffnet(networks,
                      diffMethod = "fisher",
                      adjust = "lfdr",
                      n1 = 20, n2 = 20)
  
  #### SUBGRAPHS FOR fecal535neg and fecal1214pos ####
  igraphs <- setNames(
    list(hb_ig, lb_ig),
    LABELS_NETWORK
  )
  metab_set <- c('fecal535neg', 'fecal1214pos', 'fecal367neg', 'fecal703neg')
  fecal535neg_subgraphs <- lapply(igraphs, \(x) {
    nb <- neighborhood(x, order = 1, nodes = V(x)[metab_set], mindist = 0)
    union_nb_v <- reduce(nb, union)
    subg <- subgraph(x, union_nb_v)
    # Select only intertype or metab -> metab edges
    subg_mm <- subgraph(subg, V(subg)[type == "b"])
    edges_mm <- E(subg_mm)
    edges_mm_mat <- as.matrix(ends(subg_mm, edges_mm))
    edges_mm_id <- get.edge.ids(subg, c(t(edges_mm_mat)))
    edges_inter <- E(subg)[inter]
    subg_m <- subgraph.edges(subg, union(edges_mm_id, edges_inter))
    # Add a color category attribute
    V(subg_m)$color_cat <- (
      metadata |> 
      group_by(nwk_vert_id) |>
      slice_max(score, with_ties=FALSE) |> filter(nwk_vert_id == "fecal776neg") |>
      ungroup() |>
      column_to_rownames("nwk_vert_id"))[V(subg_m)$name, 'color_cat']
    return(subg_m)
  })

  # Apply colors
  for (name in names(fecal535neg_subgraphs)) {
    arm_graph <- fecal535neg_subgraphs[[name]]
    color_factor  <- as.factor(V(arm_graph)$color_cat)
    node_pal      <- c(BIO_COLORS['High Bioactive'], BIO_COLORS['Low Bioactive'], "#f58231")
    node_col      <- node_pal[as.numeric(color_factor)]
    outpath <- file.path(glue("output/figures/figure_four_inset_{name}.pdf"))
    pdf(outpath, height=6, width=6, file=outpath)
    # Make custom layout
    min_pos <- min(E(arm_graph)$corr[E(arm_graph)$corr > 0]) 
    min_neg <- min(E(arm_graph)$corr[E(arm_graph)$corr < 0]) 
    # Fit all negative correlations into this space, lower being closer to 0
    new_weights <- ifelse(
      E(arm_graph)$corr > 0, 
      E(arm_graph)$corr, 
      (1 - abs(E(arm_graph)$corr / min_neg)) * min_pos 
    ) + 0.01
    # new_weights <- ifelse(E(arm_graph)$corr > 0, E(arm_graph)$corr, min(E(arm_graph)$weight) * 1.2)
    # arm_layout <- layout_with_fr(arm_graph,
    #                              weights=new_weights)
    arm_layout <- layout_with_graphopt(arm_graph)
    plot(arm_graph,
        layout=arm_layout,
        vertex.size=5,
        vertex.color=node_col,
        vertex.shape=V(arm_graph)$shape,
        vertex.label=V(arm_graph)$name |> map_chr(\(x) {
          if(grepl("fecal", x) && !(x %in% metab_set)) {
            return("")
          }
          return(x)
        }),
        vertex.border=NULL,
        vertex.label.color="black",
        vertex.label.family="sans",
        vertex.label.font=ifelse(V(arm_graph)$name %in% metab_set, 2, 1),
        vertex.label.cex=ifelse(V(arm_graph)$name %in% metab_set, .65, .5),
        vertex.label.dist=0.8,
        vertex.label.degree=-pi/2,
        edge.width=log10(abs(E(arm_graph)$corr) * 200)
    )
    dev.off()
  }

  return(
    list(
      net_compare=comp_network,
      net_compare_summary=comp_sum,
      diff_net=diff_net,
      diff_net_summary=summary(diff_net)
    )
  )
}

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