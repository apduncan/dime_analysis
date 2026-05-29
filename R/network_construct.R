se_params <- function() {
  #' SPIEC-EASI parameters to use when learning the paired network
  list(
    nlambda = 100,
    lambda.min.ratio = 1e-1,
    method = "glasso",
    sel.criterion = "bstars",
    pulsar.params = list(rep.num = 100, ncores = 1),
    verbose = TRUE
  )
}

network_prevalence_filter <- function(
  abundance,
  prevalence = 0.2,
  func_filt = TRUE,
  sample_md = NULL
) {
  #' Filter out any features which appear in very few samples, and select
  #' only those which seem to have some relationship to dietary conditions
  #'
  #' @param abundance Abundance as matrix
  #' @param prevalence Minimum permitted prevalence
  #' @param func_filt Apply custom filtering to PFAM function data
  #'
  # Normalise abundance
  norm_abd <- abundance |>
    as.matrix() |>
    prop.table(margin = 2) |>
    as.data.frame()
  stopifnot(all.equal(
    colSums(norm_abd) |> unname(),
    rep(1, dim(norm_abd)[2]))
  )

  prev_obs <- (abundance > 0) |> rowSums()
  thresh <- prevalence * dim(abundance)[2]
  feat_pass <- prev_obs > thresh
  log_info("Input dimensions: {paste0(dim(abundance), collapse = ', ')}")
  log_info("Filter prevalence: {prevalence} ({thresh} samples)")
  log_info("Accept {sum(feat_pass)} features")
  # Use this as indexer to original df
  abundance_filtered <- abundance[feat_pass, ]

  # Custom filtering for high-dimensional function data
  if (func_filt) {
    # Select only features which appear to have some relationship to
    # the dietary conditions. Functional annotation is too high dimensional
    # for network construction
    alpha <- 0.1
    max_features <- 1500

    test_abd <- norm_abd[, sample_md |>
      filter(grepl("after", sample_arm)) |>
      arrange(sample_arm, participant) |>
      pull(1)]
    # Apply paired Wilcoxon tests to identify features potentially related
    # to the dietary conditions
    split_factor <- sample_md |>
      filter(grepl("after", sample_arm)) |>
      arrange(sample_arm, participant) |>
      pull(sample_arm) |>
      as.factor()
    tests <- apply(
      test_abd |> as.matrix(),
      FUN = \(x) {
        splits <- split(x, split_factor)
        return(wilcox.test(splits[[1]], splits[[2]], paired = TRUE)$p.value)
      },
      MARGIN = 1
    )
    # Select all features p <= alpha, limitting to top max_features
    selected_features <- tests[tests <= alpha]
    log_info("Selected {length(selected_features)} features")
    selected_features <- ordered(selected_features)[
      1:min(length(selected_features), max_features)
    ]

    # Select down to significant
    abundance_filtered <- abundance_filtered[
      rownames(abundance_filtered) %in% names(selected_features),
    ]
    norm_filtered <- norm_abd[
      rownames(norm_abd) %in% names(selected_features),
    ]
  }
  abundance_filtered
}

network_split_matrix <- function(
  mat,
  metadata
) {
  #' Split a matrix with samples on columns and features on rows into tables
  #' with just samples from aftr the high and low bioactive interventions.
  #'
  #' @param mat Matrix with features on rows, samples on columns
  #' @param metadata Sample metadata

  #' High
  tbl_meta_high <- metadata |>
    filter(sample_arm == "after_high")
  mat_high <- mat[, tbl_meta_high |> pull(sample_id)]
  #' Low
  tbl_meta_low <- metadata |>
    filter(sample_arm == "after_low")
  mat_low <- mat[, tbl_meta_low |> pull(sample_id)]

  # Remove all features which are 0 in both matrices
  feat_keep <- !((mat_low |> rowSums() == 0) & (mat_high |> rowSums() == 0))
  mat_low <- mat_low[feat_keep, ]
  mat_high <- mat_high[feat_keep, ]

  list(
    high = mat_high,
    low  = mat_low
  )
}

paired_network_se <- function(
  mat_a,
  mat_b,
  se_params
) {
  #' Build paired SPIEC-EASI network using the two datatypes mode, in this
  #' study community KO abundance an untargetted metabolomics peaks.
  #'
  #' @param mat_a Abundance matrix A (here, KO abundance)
  #' @param mat_b Abundance matrix B (here, untargetted metabolomics)
  #' @param se_params List of parameters to passed to spiec.easi() call

  mat_a <- mat_a |> t()
  mat_b <- mat_b |> t()

  log_info("Input dimensions")
  log_info("Matrix A: {paste0(dim(mat_a), collapse=', ')}")
  log_info("Matrix B: {paste0(dim(mat_b), collapse=', ')}")

  # For paired data, it's possible that measurements don't exist for some
  # samples so we should subset to the intersection
  log_info("Drop samples which are not in both data types")
  in_both <- intersect(rownames(mat_a), rownames(mat_b))
  mat_a <- mat_a[in_both, ]
  mat_b <- mat_b[in_both, ]
  log_info("Dimensions after dropping")
  log_info("Matrix A: {paste0(dim(mat_a), collapse=', ')}")
  log_info("Matrix B: {paste0(dim(mat_b), collapse=', ')}")

  # Check samples are using the same names
  if (!setequal(rownames(mat_a), rownames(mat_b))) {
    stop("Input matrices have different sample names")
  }
  # Change to same ordering if needed
  if (!all(rownames(mat_a) == rownames(mat_b))) {
    mat_b <- mat_b[match(rownames(mat_a), rownames(mat_b))]
  }

  # Filter out any completely unknown features (those labelled ?)
  mat_a <- mat_a[, colnames(mat_a) != "?"]
  mat_b <- mat_b[, colnames(mat_b) != "?"]

  #### BUILD NETWORK ####
  # To run on paired data, just provided a length 2 list for data
  se_params$data <- list(mat_a, mat_b)

  # Run with our merged parameters
  log_info("Building network...")
  se_net <- do.call(spiec.easi, args = se_params)
  log_info("Build complete")

  #### MAKE IGRAPH OBJECT ####
  # Create a nicely labelled igraph object while we have all the data loaded
  log_info("Converting to igraph object")
  dtype <- c(rep("a", ncol(mat_a)), rep("b", ncol(mat_b)))

  # Make igraph object and set some properties
  ig <- adj2igraph(make_symmetric(as.matrix(getRefit(se_net))))

  # Vertex properties
  V(ig)$name <- c(colnames(mat_a), colnames(mat_b))
  V(ig)$type <- dtype
  # Assume we want to use shape to indicate data type later
  V(ig)$shape <- ifelse(dtype == "a", "circle", "square")

  # Edge properties
  # Weight
  sr_corr <- cov2cor(as.matrix(getOptCov(se_net)))
  rownames(sr_corr) <- c(colnames(mat_a), colnames(mat_b))
  colnames(sr_corr) <- c(colnames(mat_a), colnames(mat_b))
  edge_ends <- as.matrix(ends(ig, E(ig)))
  edge_weights <- apply(edge_ends[, ], MARGIN = 1, FUN = function(x) {
    sr_corr[x[1], x[2]]
  })
  E(ig)$weight <- abs(edge_weights)
  E(ig)$corr <- edge_weights
  # Boolean to indicate pos/neg coefficient
  E(ig)$pos <- edge_weights > 0
  E(ig)$color <- ifelse(edge_weights > 0, "forestgreen", "brown2")
  # Inter-domain edge (boolean)
  a_cutoff <- ncol(mat_a)
  edge_ends_id <- ends(ig, E(ig), name = FALSE)
  edge_inter <- apply(edge_ends_id[, ], MARGIN = 1, FUN = function(x) {
    (((x[1] > a_cutoff) && (x[2] <= a_cutoff)) ||
      ((x[1] <= a_cutoff) && (x[2] > a_cutoff)))
  })
  E(ig)$inter <- edge_inter

  # #### WRITE RESULTS ####
  cor <- stats::cov2cor(as.matrix(getOptCov(se_net)))
  ass_mat <- as.matrix(cor * SpiecEasi::getRefit(se_net))
  diag(ass_mat) <- 1
  colnames(ass_mat) <- colnames(se_net$est$data)
  rownames(ass_mat) <- colnames(se_net$est$data)
  ass_mat <- as.data.frame(ass_mat) %>% rownames_to_column("id")
  list(
    se_net = se_net,
    igraph = ig,
    adjacency = ass_mat
  )
}

network_write <- function(
  lst_res,
  output_root,
  condition_name
) {
  #' Output network results to file with naming conventions which will match
  #' results distributed in Zenodo data.
  #'
  #' @param lst_res Results from paired_network_se
  #' @param output_root Directory to write results
  #' @param condition_name Descriptive name to use in output filenames

  dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
  ig_out <- file.path(
    output_root,
    glue("{condition_name}.igraph.Rds")
  )
  assoc_out <- file.path(
    output_root,
    glue("{condition_name}.association.Rds")
  )
  se_out <- file.path(
    output_root,
    glue("{condition_name}.se.Rds")
  )
  saveRDS(lst_res$igraph, ig_out)
  saveRDS(lst_res$adjacency, assoc_out)
  saveRDS(lst_res$se_net, se_out)
  return(c(
    ig_out, assoc_out, se_out
  ))
}

# SPIEC-EASI results in a non-symmetric adjacency matrix
# Make symmetric from lower triangular matrix
# This function is taken from the fifer package,
# https://github.com/dustinfife/fifer/

#' Force a matrix to be symmetric
#'
#' @param a A matrix you wish to force to be symmetrical
#' @param lower.tri Should the upper triangle be replaced with the lower
#' triangle?
#' @return a symmetric matrix
#' @author Dustin Fife
#' @export
#' @examples
#' a <- matrix(rnorm(16), ncol = 4)
#' make.symmetric(a, lower.tri = FALSE)
make_symmetric <- function(a, lower.tri = TRUE) {
  if (lower.tri) {
    ind <- upper.tri(a)
    a[ind] <- t(a)[ind]
  } else {
    ind <- lower.tri(a)
    a[ind] <- t(a)[ind]
  }
  a
}
