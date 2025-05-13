taxa_distance <- function(
  tbl_taxa_abundance
) {
  #' Bray-Curtis distance between samples
  #'
  #' Input abundances are total-sum-scaled, then Bray-Curtis distance
  #' calculated. No filtering of taxa is performed.
  #'
  #' @param  tbl_taxa_abundance Taxonomic abundance, samples on columns as
  #' tibble.
  #' @returns Distance matrix produced by vegdist

  tbl_taxa_abundance |>
    column_to_rownames(colnames(tbl_taxa_abundance)[[1]]) |>
    as.matrix() |>
    prop.table(margin = 2) |>
    t() |>
    vegdist(method = "bray")
}

taxa_distance_abs <- function(
  tbl_taxa_abundance
) {
  #' Bray-Curtis distance between samples
  #'
  #' Input abundances have Bray-Curtis distance calculated.
  #' No filtering of taxa is performed.
  #'
  #' @param  tbl_taxa_abundance Taxonomic abundance, samples on columns as
  #' tibble.
  #' @returns Distance matrix produced by vegdist
  #' tbl_taxa_abundance |>
  a <- tbl_taxa_abundance |>
    column_to_rownames(colnames(tbl_taxa_abundance)[[1]]) |>
    as.matrix() |>
    t() |>
    vegdist(method = "bray")
}

taxa_tss <- function(
  tbl_taxa_abundance
) {
  #' Total sum scale taxonominc abundance
  #'
  #' @param tbl_taxa_abundance Taxonominc abundance as tibble, samples on
  #' columns
  #' @returns TSS taxonomic abundance as tibble
  tbl_taxa_abundance |>
  column_to_rownames(colnames(tbl_taxa_abundance)[[1]]) |>
  as.matrix() |>
  prop.table(margin = 2) |>
  data.frame(check.names = FALSE) |>
  rownames_to_column(colnames(tbl_taxa_abundance)[[1]]) |>
  as_tibble()
}

enterosignature_reapply <- function(
  pth_genus
) {
  #' Generate 5ES model weights using cvanmf
  #'
  #' From genus level abundance calculate Enterosignature weights. This is
  #' done using a shell command to the cvanmf CLI. Targets should therefore
  #' be run in an environment where these CLI tools are available.
  #'
  #'
  system(
    glue("reapply -i {pth_genus} -m 5es -o output/models/es/")
  )
  return("output/models/es/")
}

read_es <- function(
  pth_es
) {
  #' Read all tables in the ES output
  #'
  #' @param pth_es Directory containing cvanmf output
  #' @returns List with names matching file names (H, W, model_fit,...) of
  #' tables
  
  # Probably some smarter map way to do this, but I'm slow today...
  olst <- list()
  for (pth in Sys.glob("output/models/es/*.tsv")) {
    name <- basename(pth) |> str_replace(".tsv", "")
    olst[[name]] <- read_delim(pth)
  }
  # Uppercase sample names
  colnames(olst$h) <- toupper(colnames(olst$h))
  olst
}

enterosignature_dbrda <- function(
  tbl_h,
  tbl_w,
  tbl_md
) {
  #' dbRDA for Enterosignature weights
  #'
  #' dbRDA using capscale conditioned on participant and constrained on diet.
  #' Restricted to baseline, HB and LB diets.
  #'
  #' @param tbl_h H matrix from cvanmf as tibble
  #' @param tbl_w W matrix from cvanmf as tibble
  #' @param tbl_md Metadata as tibble
  #' @retruns Results of capscale()

  # Using df as shorthand for "has row names"
  df_md_sub <- tbl_md |>
    remove_before_metadata() |>
    filter(sample_id %in% colnames(tbl_h)) |>
    column_to_rownames("sample_id")
  df_h <- tbl_h |>
    column_to_rownames("...1")
  # Subset H matrix to only BLN, HB, LB
  df_h <- df_h[, rownames(df_md_sub)]
  # Order H matrix and metadata to match
  df_md_sub <- df_md_sub[colnames(df_h), ]
  # es_dist <- vegdist(df_h |> t(), method = "euclidean")
  df_use <- df_h |> t()
  es_dbrda <- capscale(
    df_use ~ Diet + Condition(participant), df_md_sub, distance = "euclidean"
  )
  return(es_dbrda)
}

esdbrda_metab_corr <- function(
  dbrda_es,
  mat_metabolite,
  tbl_md,
  alpha = 0.05,
  ...
) {
  #' Correlation of metabolites to ES dbRDA
  #'
  #' Produce correlation of metabolites to ordination. Identify which of the
  #' species vectors each metabolite vector is closest to and count.
  #'
  #' @returns List with 'vectors' a dataframe of vectors, with column 'closest'
  #' indicating the species the metabolite is closest to, and 'angle' with the 
  #' angle to that centroid.

  # Subset metabolites to samples in dbrda
  mat_metab_prune <- mat_metabolite[
    ,rownames(scores(dbrda_es)$sites)] |> t() |> as.data.frame()
  metab_fit <- envfit(
    dbrda_es, mat_metab_prune, display = "lc", permutations = 1000
  )
  # Calculate angle to each centroid
  df_species <- scores(dbrda_es)$species
  df_vectors <- metab_fit$vectors$arrows |> as.data.frame()
  dim_use <- colnames(df_vectors)
  df_vectors$pval <- metab_fit$vectors$pvals
  df_vectors <- df_vectors |>
    filter(pval <= alpha)
  pairs <- expand.grid(rownames(df_species), rownames(df_vectors))
  colnames(pairs) <- c("species", "metab_vec")
  pairs$cosine <- apply(pairs, MARGIN = 1, FUN = \(x) {
    a <- df_species[,dim_use][x[[1]],] |> as.vector()
    b <- df_vectors[,dim_use][x[[2]],] |> unlist() |> as.vector()
    cosine(a, b)
  })
  # Identify max centroid angle
  summary <- pairs |>
    group_by(metab_vec) |>
    filter(cosine == max(cosine)) |>
    ungroup() |>
    group_by(species) |>
    summarise(max_cos = n())
  summary$max_cos_pc <- summary$max_cos / length(unique(pairs$metab_vec))
  return(
    list(
      full = pairs,
      summary = summary
    )
  )
}

plot_es_dbrda <- function(
  dbrda_es,
  tbl_metadata
) {
  #' Plot ES dbRDA results
}

prevalence_filter <- function(
  tbl,
  prevalence
) {
  #' Remove features with low prevalence
  #'
  #' Any features with non-zero values in less than or equal to 'prevalence'
  #' samples will be dropped.
  #' @param tbl Feature tibble, features on rows, first column feature name
  #' @param prevalence Number of samples (not proportion) as threshold
  fname <- colnames(tbl)[1]
  df <- tbl |> column_to_rownames(fname)
  features_dropped <- rownames(df)[
    rownames(df)[rowSums(df > 0) <= prevalence]
  ]
  log_warn(
    "{length(features_dropped)} dropped by prevalence threshold {prevalence}")
  log_info("Dropped: {features_dropped}")
  return(
    df[setdiff(rownames(df), features_dropped), , drop = FALSE] |>
      rownames_to_column(fname)
  )
}

rarefy_tbl <- function(
  tbl,
  depth = NULL,
  repeats = 100,
  threads = 8,
  seed = 9898
) {
  #' Rarefy a table and produced alpha diversity using rtk
  #'
  #' Perform rarefaction on a count, or count like, tibble. This uses the rtk
  #' package, which will also calculate alpha diversity measures.
  #' @param tbl Feature tibble, samples on columns, first column feature name
  #' @param depth Rarefaction depth. If left NULL will use minimum sample depth.
  #' Currently will raise an error if depth is greater than some existing
  #' samples.
  #' @param repeats Number of rarefactions to compute, which affects mean alpha
  #' div calculations
  #' @param threads Threads for rarefaction
  #' @param seed Seed for random subsampling
  #' @returns

  fname <- colnames(tbl)[1]
  df <- tbl |> column_to_rownames(fname)
  depth <- ifelse(is.null(depth), min(colSums(df)), depth)
  data.rmin <- rtk(
    df,
    ReturnMatrix = 1,
    depth = depth,
    repeats = repeats,
    threads = threads,
    seed = seed
  )
  rare_ko_abd <- data.rmin$raremat[[1]]
  # Get Alpha diversity metrics
  alpha_diversity <- data.frame(
    richness = get.mean.diversity(data.rmin, div = "richness"),
    shannon = get.mean.diversity(data.rmin, div = "shannon"),
    inv_simp = get.mean.diversity(data.rmin, div = "invsimpson"),
    pielou = get.mean.diversity(data.rmin, div = "eveness")
  )
  # Order out from rtk is not the same as order in to rtk when using
  # multithreaded
  rownames(alpha_diversity) <- data.rmin$divvs |>
    map(\(x) {x$samplename}) |>
    unlist()
  return(list(
    matrix = rare_ko_abd,
    alpha_diversity = alpha_diversity
  ))
}

permanova_abundance <- function(tbl_abd, tbl_md, fm_terms, seed, ...) {
  #' permANOVA of abundance
  #'
  #' Perform permANOVA using adonis2 for an abundance table
  #'
  #' @param tbl_abd Abundance tibble
  #' @param tbl_md Metadata table
  #' @param fm_terms Vector of terms to use in formula
  #' @param seed Seed set before adonis call
  #' @param ... Passed to adonis2 amtch match_abundance_metadata
  #' @return List with 'rank_name' and 'adonis2' results object

  rank_name <- colnames(tbl_abd)[[1]]
  # Match metadata and sample ordering
  df_abd <- tbl_abd |> column_to_rownames(rank_name) |> as.data.frame()
  df_md  <- tbl_md |>
    column_to_rownames(colnames(tbl_md)[[1]]) |>
    as.data.frame()
  match_dfs <- match_abundance_metdata(df_abd, df_md, rm_na = fm_terms, ...)
  df_abd <- match_dfs$abd |> t()
  df_md  <- match_dfs$md
  formula <- reformulate(termlabels = fm_terms, response = "df_abd")
  set.seed(seed)
  adonis_res <- adonis2(formula, df_md, ... )
  return(
    list(
      name = rank_name,
      result = adonis_res
    )
  )
}

match_abundance_metdata <- function(
  df_abd,
  df_md,
  rm_na = c(),
  allow_missing = FALSE,
  ...
) {
  #' Match abundance and metadata tables
  #'
  #' Restricts to only the interesction of the samples in both metadata and the
  #' abundance. Will report samples in the abundance missing in the metadata.
  #' If 'allow_missing' is FALSE, will halt. Input is expected as a data.frame,
  #' with abundance having samples on columns, metadata having samples on rows.
  #'
  #' @param df_abd Abundance data.frame, sample on columns
  #' @param df_md Metadata data.frame, samples on rows
  #' @param rm_na Remove samples which have NA values for these columns
  #' @param allow_missing If FALSE, will halt if any have no metadata
  #' @param ... Unused
  #' @returns List with 'abd', 'md'

  sample_intersect <- intersect(colnames(df_abd), rownames(df_md))
  md_missing <- setdiff(colnames(sample), sample_intersect)
  if (length(md_missing) > 0) {
    log_warn("{length(md_missing)} dropped as missing metadata: {md_missing}")
    if(!allow_missing) {
      stop("Metadata for samples missing. To allow, use allow_missing = TRUE")
    }
  }
  df_abd <- df_abd[, sample_intersect, drop = FALSE]
  df_md <- df_md[sample_intersect, , drop = FALSE]

  if (length(rm_na) > 0) {
    # Remove any samples which have NA values in the rm_na columns in metadata
    complete <- apply(is.na(df_md[, rm_na]), MARGIN = 1, FUN = \(x) {!any(x)})
    # print(complete)
    dropped <- rownames(df_md)[!complete]
    if(length(dropped) > 0) {
      log_warn("{length(md_missing)} dropped as NA in columns {dropped}: {md_missing}")
    }
    df_abd <- df_abd[, complete, drop = FALSE]
    df_md <- df_md[complete, , drop = FALSE]
  }
  return(list(abd = df_abd, md = df_md))
}

plot_alpha_diversity <- function(
  tbl_alpha_div,
  tbl_sample_metadata
) {
  #' Plot alpha diversity for LB and HB conditions
  #'
  #' @param tbl_alpha_div Table containing mean alpha diversitites from
  #' rarefaction. Expects a table with rownames being samples, and columns
  #' alpha diversity indices, with accepted column names being
  #' "inv_simp", "pielou", "richness", "shannon".

  index_names <- c(
    inv_simp = "Inv.\nSimspon",
    pielou = "Evenness",
    richness = "Richness",
    shannon = "Shannon"
  )
  # Restict to only those in the table
  index_use <- index_names[
    intersect(names(index_names), colnames(tbl_alpha_div))
  ]

  sample_md <- tbl_sample_metadata
  alpha_test <- tbl_alpha_div |>
    rownames_to_column("sample") %>%
    pivot_longer(-sample,
               names_to = "index",
               values_to = "value") %>%
    left_join(sample_md, by = join_by(sample == sample_id))

  plt_ad <- alpha_test %>%
    filter(sample_arm %in% c("after_high", "after_low")) %>%
    arrange(participant, sample_arm) %>%
    mutate(diet = sample_arm |> map_chr(\(x) {
      switch(x,
            "after_high" = "High Bioactive",
            "after_low" = "Low Bioactive")
    }),
    index = fct_relevel(index, names(index_use))) %>%
    ggplot(aes(x = diet, y = value)) +
    geom_boxplot(aes(fill = diet, color = diet), alpha = 0.3) +
    facet_wrap(
      vars(index),
      ncol = 4,
      scales = "free_y",
      labeller = labeller(index = index_use)
    ) +
    geom_point(aes(color = diet), size = 1) +
    geom_line(aes(group = participant), linewidth = 0.1) +
    stat_compare_means(comparisons = list(c(1, 2)), paired = TRUE, size = 2) +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    scale_x_discrete(limits = c("Low Bioactive", "High Bioactive"),
                    labels = NULL,
                    breaks = NULL) + 
    scale_color_manual(values = BIOACTIVE_COLORS, name = "Diet") +
    scale_fill_manual(values = BIOACTIVE_COLORS, name = "Diet") +
    ylab(NULL) + xlab(NULL) +
    theme(legend.position="bottom") +
    ggtitle("Taxonomic Diversity")
  return(plt_ad)
}

plot_beta_within <- function(
  dst,
  tbl_sample_metadata,
  only_baseline = FALSE
) {
  bc_dis <- dst
  sample_md <- tbl_sample_metadata

  # Convert to long form
  bd_long <- bc_dis %>%
    as.matrix() %>%
    (\(x) {
      x[lower.tri(x, diag = FALSE)] <- NA
      return(x)
      }
    ) %>% as.data.frame() %>%
    rownames_to_column("a") %>%
    pivot_longer(!a, values_to = "bc_dist", names_to = "b") %>%
    filter(!is.na(bc_dist))

  # Select only the samples which are in the same arm
  bd_within_arm <- bd_long %>%
    # Remove the diagonal
    filter(a != b) %>%
    left_join(sample_md %>% dplyr::select(sample_id, sample_arm), 
              by = join_by(a == sample_id)) %>%
    rename(arm_a = sample_arm) %>%
    left_join(sample_md %>% dplyr::select(sample_id, sample_arm), 
              by = join_by(b == sample_id)) %>%
    rename(arm_b = sample_arm) %>%
    # Label any V1 samples as baseline
    mutate(arm_a = ifelse(grepl("V1", a), "baseline", arm_a),
          arm_b = ifelse(grepl("V1", b), "baseline", arm_b)) %>%
    # Only within the same arm, and remove diagonal
    filter(arm_a == arm_b & a != b)

  # Plot the beta diversity distributions for baseline, high, and low diets
  # Remove any before_ samples which are not at timepoint V1
  fig2_data <- bd_within_arm %>%
    mutate(Diet = arm_a |> map_chr(arm_to_diet)) %>%
    filter(!is.na(Diet)) %>%
    mutate(Diet = ifelse(Diet == "Baseline", "Non-Intervention", Diet))
  grey_label <- "Non-Intervention"
  if (only_baseline) {
    fig2_data <- fig2_data |>
      filter(!grepl("before", arm_a))
    grey_label <- "Baseline"
  }
    # Can optionally exclude 20-V1 as it is based on an imputed cell count
    # filter(!(grepl("20-V1", a) | grepl("20-V1", b))) %>%
  fig2_beta_div <- fig2_data |> 
    ggplot(
    aes(x = Diet, fill = Diet, y = bc_dist, color = Diet)) +
    geom_violin(alpha = 0.5, width = 0.5, show.legend = FALSE) +
    # geom_jitter(aes(color = Diet), size = 0.1) +
    geom_boxplot(notch = TRUE, width = 0.1, alpha = 0.5, 
                outlier.size = 2, outlier.stroke = 0, outlier.alpha = 1) +
    scale_fill_manual(
      values = c("High Bioactive" = BIOACTIVE_COLORS[2] |> unname(),
                "Low Bioactive" = BIOACTIVE_COLORS[1] |> unname(),
                grey_label = BIOACTIVE_COLORS[3] |> unname())
    ) +
    scale_color_manual(
      values = c("High Bioactive" = BIOACTIVE_COLORS[2] |> unname(),
                "Low Bioactive" = BIOACTIVE_COLORS[1] |> unname(), 
                grey_label = BIOACTIVE_COLORS[3] |> unname())
    ) +
    stat_compare_means(comparisons = list(
      c(grey_label, "Low Bioactive"),
      c(grey_label, "High Bioactive"),
      c("High Bioactive", "Low Bioactive")),
      method = "wilcox.test",
      size = 2
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(
      limits = c(grey_label, "Low Bioactive", "High Bioactive"),
      labels = NULL,
      breaks = NULL
    ) +
    ggtitle("Beta Diversity") +
    ylab("Bray-Curtis\nDissimilarity") +
    xlab(NULL) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())
  fig2_beta_div
}

scale_taxa <- function(
  tbl_taxa,
  tbl_cell_count
) {
  #' Scale relative abundance to absolute using cell counts
  #'
  #' @param tbl_taxa Taxa as tibble, first column taxa names
  #' @param tbl_cell_count Cell counts, first column sample names, second counts

  ad_fc <- tbl_cell_count
  rare_rel <- prop.table(
    tbl_taxa |>
      column_to_rownames(colnames(tbl_taxa)[1]) |>
      as.matrix(check.names = FALSE),
    margin = 2
  )
  # Put in same order as scaling matrix
  rare_rel <- rare_rel[, ad_fc$sample_id]
  rare_scaled = sweep(
    rare_rel %>% as.matrix(),
    MARGIN = 2,
    STATS = ad_fc$mean_norm_cellcount,
    FUN = "*"
  ) %>%
  as.data.frame()
  return(
    rare_scaled |>
      rownames_to_column(colnames(tbl_taxa)[1])
  )
}

diet_dbrda <- function(
  abundance,
  metadata,
  n_species_vectors = 5,
  centroids = TRUE,
  condition = TRUE
  ) {
  es_dbrda <- list()
  es_dbrda$long_table <- abundance |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "feature", values_to = "rel_weight") |>
    left_join(metadata) |>
    # Restrict to just three timepoints of interest - BLN, HB, and LB
    filter(time_point %in% c("V1", "V2", "V4")) |>
    mutate(Diet = sample_arm |> map_chr(\(x) {
      if(grepl("before", x)) {
        return("Baseline")
      }
      if(grepl("high", x)) {
        return("High Bioactive")
      }
      if(grepl("low", x)) {
        return("Low Bioactive")
      }
      stop("Sample arm not correct, should be baseline, after_high or after_low")
    }))
  
  # Selected down to just points of interest, for capscale need two tables
  # - abundance matrix, - metadata table
  es_dbrda$matrix <- es_dbrda$long_table |>
    select(sample_id, feature, rel_weight) |>
    pivot_wider(names_from = feature, values_from = rel_weight, id_cols = sample_id) |>
    column_to_rownames("sample_id") |>
    as.matrix()
  es_dbrda$metadata <- es_dbrda$long_table |>
    select(sample_id, Diet, participant) |>
    distinct()
  # Ensure ordering is correct
  stopifnot(all(rownames(es_dbrda$matrix) == es_dbrda$metadata$sample_id))
  
  # Run dbrda
  if(condition) {
    es_dbrda$dbrda <- capscale(es_dbrda$matrix ~ Diet + Condition(participant),
                               es_dbrda$metadata,
                               distance = "bray")
  } else {
    es_dbrda$dbrda <- capscale(es_dbrda$matrix ~ Diet,
                               es_dbrda$metadata,
                               distance = "bray")
  }
  es_dbrda$summary <- summary(es_dbrda$dbrda)
  
  # Make a plot of results
  dbrda_dfs <- plot(es_dbrda$dbrda, scaling = 2, display = c("sites", "species", "bp"),
                    tidy=TRUE)
  es_dbrda$tidy <- dbrda_dfs
  
  biplot_mult <- attr(dbrda_dfs$biplot, "arrow.mul")
  
  # Derive percentages of conditioned variance to label axes
  total_inertia <- es_dbrda$summary$tot.chi
  conditioned_inertia <- es_dbrda$summary$partial.chi
  constrained_inertia <- es_dbrda$summary$constr.chi
  unconstrained_interia <- es_dbrda$summary$unconst.chi
  
  es_dbrda$plot2 <- dbrda_dfs$sites |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    left_join(es_dbrda$metadata) %>%
    ggplot() +
    geom_point(aes(x = CAP1, y = CAP2, color = Diet), size = 1) +
    scale_color_manual(values = BIOACTIVE_COLORS) +
    stat_ellipse(aes(x = CAP1, y = CAP2, color = Diet)) +
    # Biplot arrows
    geom_segment(data = dbrda_dfs$biplot |> data.frame() |>
                   rownames_to_column("sample_id") |>
                   left_join(metadata),
                   mapping = aes(
                     xend = CAP1 * biplot_mult, 
                     yend = CAP2 * biplot_mult,
                     color = diet), 
                    x = 0, y = 0)
  
  es_dbrda$plot <- es_dbrda$summary$site |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    left_join(es_dbrda$metadata) %>%
    ggplot() +
    geom_point(aes(x = CAP1, y = CAP2, color = Diet), size = 1) +
    scale_color_manual(values = BIOACTIVE_COLORS) +
    stat_ellipse(aes(x = CAP1, y = CAP2, color = Diet))
  
  # Add any requested vectors - select by vector norm on the constrained axes
  vector_length <- apply(es_dbrda$summary$species[,1:2],
        FUN = \(x) { sum(x^2) },
        MARGIN = 1)
  vector_order <- order(vector_length, decreasing = TRUE)
  vector_length_ordered <- cbind(es_dbrda$features[vector_order],
                                 vector_length[vector_order])
  # Select the top n longest vectors, or number of vectors in frame if less than top_n
  es_dbrda$cap_vector_lengths <- vector_length_ordered
  
  # Find the cosine similarity to each of the centroids
  centroid_cosines <- list()
  for(centroid in rownames(es_dbrda$summary$centroids)) {
    centroid_vector <- es_dbrda$summary$centroids[centroid, 1:2]
    centroid_cosines[[centroid]] <- apply(
      es_dbrda$summary$species[, 1:2],
      FUN = \(x) {cosine(x, centroid_vector)},
      MARGIN = 1
    )
  }
  centroid_cosines$norm <- vector_length
  centroid_cosine_df <- do.call(cbind,centroid_cosines)
  
  # Cosine similarity to each of biplot vectors
  biplot_cosines <- list()
  bp_df <- dbrda_dfs$biplot |> data.frame()
  for(bp_name in rownames(bp_df)) {
    bp_vector <- as.numeric(bp_df[bp_name, 1:2])
    biplot_cosines[[bp_name]] <- apply(
      es_dbrda$summary$species[, 1:2],
      FUN = \(x) {cosine(x, bp_vector)},
      MARGIN = 1
    )
  }
  biplot_cosines$norm <- vector_length
  es_dbrda$biplot_cosines <- do.call(cbind, biplot_cosines) |> as.data.frame()
  
  top_feats <- es_dbrda$summary$species[
    vector_order[1:min(n_species_vectors, dim(es_dbrda$summary$species)[1])],
    1:2]
  es_dbrda$cap_vector_lengths <- centroid_cosine_df
  if(n_species_vectors > 0) {
    es_dbrda$plot <- es_dbrda$plot +
      geom_segment(data = top_feats,
                   mapping = aes(xend = CAP1, yend = CAP2), x = 0, y = 0)
  }
  
  if(centroids) {
    # Add constraint group centroids
    es_dbrda$plot <- es_dbrda$plot +
      geom_text(data = es_dbrda$summary$centroids |> 
                   as.data.frame() |>
                   rownames_to_column("Diet") |>
                   mutate(Diet = gsub("Diet", "", Diet)),
                 mapping = aes(x = CAP1, y = CAP2, label = Diet),
                size = 3,
                color = BIO_TEXT_COLOR
                 )
  }
  return(es_dbrda)
}

decorate_species_dbrda <- function(
  lst_dbrda
) {
  spec_db <- lst_dbrda
  spec_dbrda_fit <- envfit(spec_db$dbrda, spec_db$matrix |> as.data.frame())
  spec_dbrda_fit_sig_taxa <- spec_dbrda_fit$vectors$pvals[spec_dbrda_fit$vectors$pvals <= 0.05]

  # Get and scale arrows
  all_arrows <- scores(spec_dbrda_fit, display="vectors")
  mult <- ordiArrowMul(spec_dbrda_fit, fill=2)
  scaled_arrows <- all_arrows * mult
  trim_tax <- function(x) {
    lineage <- strsplit(x, split=";", fixed=TRUE) |> unlist()
    # Find last non-? and concatenate
    last_known <- (lineage != "?") |> sum()
    txt <- paste0(lineage[last_known:length(lineage)], collapse=";")
    return(txt)
  }
  sig_arrows <- scaled_arrows[
    (spec_dbrda_fit$vectors$pvals[spec_dbrda_fit$vectors$pvals <= 0.05]) |> names(),
    ]
  rownames(sig_arrows) <- rownames(sig_arrows) |> map_chr(trim_tax)
  library(ggrepel)
  spec_db_w_arrows <- spec_db$plot +
    geom_segment(aes(xend=CAP1, yend=CAP2),
                x=0, y=0,
                data = sig_arrows |> as.data.frame() |> rownames_to_column("taxon"),
                linewith=0.3) +
    geom_text_repel(aes(x=CAP1, y=CAP2, label=taxon),
              data = sig_arrows |> as.data.frame() |> rownames_to_column("taxon"),
              size=2)
  return(spec_db_w_arrows)
}

pcoa <- function(
  tbl,
  method = "bray"
) {
  #' Run standard PCoA using cmdscale and vegdist
  #' 
  #' @param tbl Abundances as tibble, first column feature names, samples
  #' on columns
  #' @param method Distance method

  d <- vegdist(
    tbl |>
      column_to_rownames(colnames(tbl)[[1]]) |>
      t(),
    method = method
  )
  return(cmdscale(d, eig=TRUE))
}

plot_pcoa <- function(
  pcoa_res,
  tbl_sample_metadata
) {
  #' Plot PCoA ordination
  #'
  #' @param pcoa_res Result from cmdscale
  #' @param tbl_sample_metadata Sample metadata

  tbl_coords <- pcoa_res$points |>
    as.data.frame(check.names = FALSE) |>
    rownames_to_column("sample_id") |>
    rename(PC1 = 2, PC2 = 3) |>
    left_join(tbl_sample_metadata |> remove_before_metadata())
  
  tbl_coords |>
    ggplot(
      aes(
        x = PC1,
        y = PC2,
        color = Diet
      )
    ) +
    geom_line(
      data = tbl_coords %>% filter(Diet %in% c("High Bioactive", "Baseline")),
      mapping = aes(
        x = PC1,
        y = PC2,
        color = Diet,
        group = participant
      ),
      color = BIOACTIVE_COLORS[1],
      linewidth = 0.15
    ) +
    geom_line(
      data = tbl_coords %>% filter(Diet %in% c("Low Bioactive", "Baseline")),
      mapping = aes(
        x = PC1,
        y = PC2,
        color = Diet,
        group = participant
      ),
      color = BIOACTIVE_COLORS[2],
      linewidth = 0.15
    ) +
    geom_point() +
    scale_color_manual(
      values = BIOACTIVE_COLORS3
    )
}

remove_before_samples <- function(
  tbl,
  tbl_sample_metadata
) {
  after_only <- remove_before_metadata(tbl_sample_metadata)
  keep <- colnames(tbl) %in% after_only$sample_id
  keep[[1]] <- TRUE
  return(
    tbl[, keep]
  )
}

paper_figure_two <- function(
  plt_species_pcoa,
  plt_species_dbrda,
  plt_tax_shannon,
  plt_fun_shannon,
  plt_tax_bray,
  plt_fun_bray,
  plt_mgs_cazymes
) {
  layout <- "
  AABBCDEF
  GGGGGGGG
  "
  remove_guides <- guides(
    color = "none",
    fill = "none",
    shape = "none",
    size = "none"
  )

  # Modify a few plots to look better when patchworked
  plt_pcoa_mod <- plt_species_pcoa +
    ggtitle("Species PCoA") +
    theme(legend.position = "bottom")
  plt_dbrda_mod <- plt_species_dbrda +
    labs(
      title = "Species dbRDA",
      subtitle = "Formula: Condition(individual) + Diet"
    ) +
    remove_guides
  plt_taxad_mod <- plt_tax_shannon +
    labs(
      title = "Alpha Diversity",
      subtitle = "Taxonomy"
    ) +
    ylab("Shannon") +
    theme(strip.text = element_blank()) +
    remove_guides
  plt_funad_mod <- plt_fun_shannon +
    labs(
      title = "",
      subtitle = "Function"
    ) +
    ylab("") +
    theme(strip.text = element_blank()) +
    remove_guides
  plt_taxbd_mod <- plt_tax_bray +
    labs(
      title = "Differences with Diet",
      subtitle = "Taoxonomy"
    ) +
    remove_guides
  plt_funbd_mod <- plt_fun_bray +
    labs(
      title = "",
      subtitle = "Function"
    ) +
    ylab("") +
    remove_guides
  plt_cazy_mod <- plt_mgs_cazymes +
    labs(
      title = "CAZyme Substrate Usage in Diet Associated MGS"
    )
  return(
    wrap_plots(
      plt_pcoa_mod,
      plt_dbrda_mod,
      plt_taxad_mod,
      plt_funad_mod,
      plt_taxbd_mod,
      plt_funbd_mod,
      plt_mgs_cazymes,
      design = layout
    ) +
    plot_layout(
      guide = "collect",
    ) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "", "D", "", "E"))
    ) &
    theme(legend.position = "bottom") &
    THEME_DIME
  )
}

paper_figure_two_alternate <- function(
  plt_species_pcoa,
  plt_species_dbrda,
  plt_tax_shannon,
  plt_fun_shannon,
  plt_tax_bray,
  plt_fun_bray,
  plt_mgs_cazymes
) {
  layout <- "
  AABBCCDD
  EFGGGGGG
  "
  remove_guides <- guides(
    color = "none",
    fill = "none",
    shape = "none",
    size = "none"
  )

  # Modify a few plots to look better when patchworked
  plt_pcoa_mod <- plt_species_pcoa +
    ggtitle("Species PCoA") +
    theme(legend.position = "bottom")
  plt_dbrda_mod <- plt_species_dbrda +
    labs(
      title = "Species dbRDA",
      subtitle = "Formula: Condition(individual) + Diet"
    ) +
    remove_guides
  plt_taxad_mod <- plt_tax_shannon +
    labs(
      title = "Alpha Diversity",
      subtitle = "Taxonomy"
    ) +
    theme(
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 8)
    ) +
    remove_guides
  plt_funad_mod <- plt_fun_shannon +
    labs(
      title = "",
      subtitle = "Function"
    ) +
    theme(
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 8)
    )
  plt_taxbd_mod <- plt_tax_bray +
    labs(
      title = "Differences with Diet",
      subtitle = "Taoxonomy"
    ) +
    remove_guides
  plt_funbd_mod <- plt_fun_bray +
    labs(
      title = "",
      subtitle = "Function"
    ) +
    ylab("") +
    remove_guides
  plt_cazy_mod <- plt_mgs_cazymes +
    labs(
      title = "CAZyme Substrate Usage in Diet Associated MGS"
    )
  return(
    wrap_plots(
      plt_pcoa_mod,
      plt_dbrda_mod,
      plt_taxad_mod,
      plt_funad_mod,
      plt_taxbd_mod,
      plt_funbd_mod,
      plt_mgs_cazymes,
      design = layout
    ) +
    plot_layout(
      guide = "collect",
    ) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "", "D", "", "E"))
    ) &
    theme(legend.position = "bottom") &
    THEME_DIME
  )
}