#' Functions to handle load, filtering and processing data

read_metadata <- function(pth) {
  #' Load DIME metadata from delimited text format
  #' 
  #' Loads metadata and determines which samples are baseline samples
  #' 
  #' @param pth Path to file
  read_delim(pth) |>
  rename(sample_id = 1) |>
  mutate(
    baseline = (
      (sample_group == "before") &
      (ifelse(sequence == "low_high", "before_low", "before_high") == sample_arm)
    )
  )
}

remove_before_metadata <- function(sample_md) {
  #' Select only post-diet and baseline metadata
  #' 
  #' Discard the pre-diet samples from metadata, and assign a long form 
  #' descriptive name in the Diet column.
  #' 
  #' @param sample_md Tibble of sample metadata
  sample_md |>
  filter(
    baseline | (sample_arm == "after_high") | (sample_arm == "after_low")
  ) |>
  mutate(
    Diet = map_chr(sample_arm, \(x) {
      if(grepl("before", x)) {
        return("Baseline")
      }
      if(x == "after_high") {
        return("High Bioactive")
      }
      return("Low Bioactive")
    })
  )
}

read_metabolites <- function(pth) {
  #' Read metabolite peak values in matrix form
  #' 
  #' @param pth Path to metabolite matrix
  read_delim(pth) |>
    column_to_rownames("gid") |>
    as.matrix()
}

read_food_composition <- function(
  pth_adj_nutrients,
  pth_nutrients
) {
  #' Residual adjusted food composition
  #'
  #' Returns three dataframes in a list:
  #' * 'all' - Nutrients & Bioactives
  #' * 'nutrients' - Nutrients only
  #' * 'bioactives' - Bioactives only
  #' @param pth_adjust_nutrients Residual adjusted nutrient & bioactive intake
  #' @param pth_nutrients Table for only nutrients. Used to identify nutrients
  #' and bioactives in full table.
  #' @returns List with tibbles of nutrient and/or bioactive intake

  tbl_res_adj <- read_delim(pth_adj_nutrients)
  tbl_nutrients <- read_delim(pth_nutrients)

  # Create sets of nutrient & bioactive column names to subset data
  nutrient_names <- tbl_nutrients |> colnames() |> 
    setdiff(c("Participants", "time_point", "sample_group", "sample_arm",
    "SAM_ID", "sample_id", "diet"))
  alldiet_names <- colnames(tbl_res_adj) |> setdiff(c("sample_id"))
  bioactive_names <- setdiff(alldiet_names, nutrient_names)

  df_all <- tbl_res_adj |>
    column_to_rownames("sample_id") |> t() |>
    order_table(rownames = TRUE)
  df_nutrient <- df_all[nutrient_names,]
  df_bioactive <- df_all[bioactive_names,]
  return(list(
    all = df_all,
    nutrient = df_nutrient,
    bioactive = df_bioactive
  ))
}

read_taxa <- function(
  pth_taxa
) {
  #' Taxonomic abundance
  #'
  #' Taxonomic abundance as a tibble (first col is taxon name)
  #' @param pth_taxa Path to MATAFILER format taxonomic abundance
  #' @returns Taxonomic abundance as tibble, samples on columns
  read_delim(pth_taxa) |>
    filter(1 != "-1") |>
    rename_all(toupper) |>
    # Filter negative control
    select(-`PLATE-1-NC`)
}

read_ko <- function(
  pth_ko
) {
  #' KEGG Ortholog abundance
  #'
  #' KO abundance as a tibble (first col is KO)
  #' @param pth_ko Path to MG-TK format taxonomic abundance
  #' @returns KO abundance as a tibble, samples on columns
  read_delim(pth_ko) |>
    filter(1 != "-1") |>
    rename_all(toupper) |>
    select(-`PLATE-1-NC`)
}

read_pfam <- function(
  pth_pfam
) {
  #' PFAM abundance
  #'
  #' PFAM abundance as a tibble (first col is PFAM)
  #' @param pth_pfam Path to MG-TK format abundance
  #' @returns PFAM abundance as a tibble, samples on columns
  read_delim(pth_pfam) |>
    filter(1 != "-1") |>
    rename_all(toupper)
}

filter_unknown <- function(
  tbl
) {
  #' Remove counts for unknown features
  #'
  #' These are collected in the -1 row in MG-TK output. Simple function
  #' but created to give some semantic meaning to what is being done.
  #'
  #' @param tbl Count tibble, samples on columns, first row labels
  #' @returns Tibble with -1 row removed
  fname <- colnames(tbl)
  return(
    tbl |>
      filter(if_any(1, ~ . != -1))
  )
}

order_table <- function(tbl, rownames = FALSE) {
  #' Arrange columns into fixed order
  #'
  #' For tables which have samples on columns, arrange the columns into a
  #' fixed order, omitting and missing columns. This is sorted 01-V1, 01-V2, ...
  #' 20-V3, 20-V4. This method assumes that the first column is a label,
  #' and all other columns correspond to samples.
  #' 
  #' @param tbl Tibble with first column label, and samples on all other columns
  #' @param rownames Boolean indicating if feature ids are in rownames. If FALSE
  #' assumes first column contains feature IDs
  #' @returns Tibble with sorted columns

  dime_order <- expand.grid(c(seq(1, 40)), c("V1", "V2", "MP", "V3", "V4")) |>
    apply(MARGIN = 1, FUN = paste0, collapse = "-") |>
    str_replace(" ", "0")
  order_intersect <- intersect(dime_order, colnames(tbl))
  if (!rownames) {
    order_intersect <- c(colnames(tbl)[[1]], order_intersect)
  }
  return(tbl[, order_intersect])
}

arm_to_diet <- function(x) {
  return(
    switch(x,
           "after_high" = "High Bioactive",
           "after_low" = "Low Bioactive",
           "before_low" = "Baseline",
           "before_high" = "Baseline",
           "baseline" = "Baseline",
           NA
    ))
}

read_cellcounts <- function(
  pth
) {

  raw <- read_delim(pth)
  # Loading flow cytometry cell counts for scaling
  ad_fc <- raw %>%
    # Participant identifiers use O rather than 0 as in microbiome
    mutate(Participants = Participants |> map_chr(function (x) {
      str_pad(gsub("O", "", x), width = 2, side = "left", pad = "0")
    })) %>%
    group_by(Participants, Visit) %>%
    summarise(mean_norm_cellcount = mean(Normalized))

  # 20-V1 is missing due to insufficient stool volume
  # Impute a value using mean of all participant 20 values
  p20_fc <- raw |> 
    mutate(Participants = Participants |> map_chr(function (x) {
      str_pad(gsub("O", "", x), width = 2, side = "left", pad = "0")
    })) |>
    group_by(Participants) |>
    summarise(total_mean = mean(Normalized)) |>
    filter(Participants == "20") |>
    pull(total_mean)

  # Add this to the ad_fc dataframe
  ad_fc <- rbind(ad_fc |> as.data.frame(), c("20", "V1", p20_fc)) |> 
    as_tibble() |> 
    mutate(mean_norm_cellcount = as.numeric(mean_norm_cellcount)) |>
    mutate(sample_id = paste0(Participants, "-", Visit))
  return(ad_fc)
}

read_dbcan <- function(
  dir_dbcan
) {
  #' Read dbcan results and metadata
  #'
  #' Return a list including counts of dbcan enzymes in each MGS, and
  #' metadata about the dbcan enzymes.
  #' @param dir_dbcan Directory with dbcan results
  #' @returns List with 'counts_raw', 'annotation_raw', 'annotation_tbl',
  #' 'counts_long' 

  dbsub_count <- read_delim(file.path(dir_dbcan, "dbsub_gene_count.csv"))
  dbsub_annot <- read_delim(file.path(dir_dbcan, "dbsub_substrate_lookup.tsv"))
  # Tidy the dbsub annotations
  # Make each entry a non-redundant list
  dbsub_annot$nr <- dbsub_annot |>
    pull("Substrate") |>
    map(strsplit, split = ", ", fixed = TRUE) |>
    map(unlist) |>
    map(unique) |>
    map(paste, collapse=",")
  expanded_annot <- map2(
    dbsub_annot$`dbCAN subfam`,
    dbsub_annot$nr,
    \(x, y) {
      substrates <- strsplit(y, ",", fixed = TRUE)
      rlist <- substrates |> unlist() |> map(\(s) {c(x, s)})
    }
  )
  expanded_annot_tbl <- do.call(rbind, expanded_annot |> flatten()) |>
    as.data.frame()
  colnames(expanded_annot_tbl) <- c("dbsub", "Substrate")
  dbsub_count_long <- dbsub_count |>
    pivot_longer(-Gene) |>
    filter(value > 0) |>
    rename(MGS = name) |>
    left_join(expanded_annot_tbl, join_by(Gene == `dbsub`), multiple = "all") |>
    filter(!is.na(Substrate))

  return(list(
    counts_raw = dbsub_count,
    annotation_raw = dbsub_annot,
    counts_long = dbsub_count_long,
    annotation_tbl = expanded_annot_tbl
  ))
}

substrate_difference <- function(
  tbl_dbsub_long,
  tbl_traits
) {
  dbsub_count_long <- tbl_dbsub_long
  traits_table <- tbl_traits
  substrate_counts <- dbsub_count_long |>
    group_by(MGS, Substrate) |>
    summarise(substrate_enzymes = sum(value)) |>
    left_join(
      traits_table |>
      mutate(
        association = ifelse(
          hb_associated == 1,
          "High Bioactive",
          ifelse(
            lb_associated == 1,
            "Low Bioactive",
            "None"
          )
        )
      ),
      join_by(MGS == Name)
    )
  
  filter_out <- c("and  beta-glucan", "glucomannan", "raffinose",
  "human milk polysaccharide")
  noassoc_color <- "#C19A6B"

  # Drop some irrelevant substrates - those which only exist in one condition
  substrate_exists <- substrate_counts |>
    group_by(association, Substrate) |>
    summarise(anyzymes = sum(substrate_enzymes) > 0) |>
    pivot_wider(
      id_cols = Substrate,
      names_from = association,
      values_from = anyzymes,
      values_fill = 0) |>
      column_to_rownames("Substrate")
  keep_substrates <- rownames(
    substrate_exists[rowSums(substrate_exists) > 2, ]
  )
  plt_substrate <- substrate_counts |>
    filter(Substrate %in% keep_substrates) |>
    # For now manually filter those with low copy number
    filter(!(Substrate %in% filter_out)) |>
    # Reorder bioactive conditions to match other plots
    ungroup() |>
    mutate(
      association = fct_relevel(
        factor(association),
        "None",
        "Low Bioactive",
        "High Bioactive"
      )
    ) |>
    ggplot(
      aes(
        x = factor(association),
        y = substrate_enzymes,
        fill = association,
        color = association
      )
    ) +
    geom_boxplot(alpha = 0.65, linewidth = 1) +
    facet_wrap(
      ~Substrate,
      scale="free_y",
      ncol = 8,
      labeller = labeller(Substrate = label_wrap_gen(18))
    ) +
    ggpubr::stat_compare_means(
      label = "p.format",
      label.x.npc = "centre",
      label.y.npc = "top",
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    scale_x_discrete(limits = c("None", "Low Bioactive", "High Bioactive")) +
    scale_fill_manual(
      values = c(
        "High Bioactive" = unname(BIOACTIVE_COLORS[2]),
        "Low Bioactive" = unname(BIOACTIVE_COLORS[1]),
        "None" = noassoc_color
      ),
      name = "Diet Association"
    ) +
    scale_color_manual(
      values = c(
        "High Bioactive" = unname(BIOACTIVE_COLORS[2]),
        "Low Bioactive" = unname(BIOACTIVE_COLORS[1]),
        "None" = noassoc_color
      ),
      name = "Diet Association"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab("") +
    ylab("Number of CAZymes in MGS") +
    ggtitle(
      "Number of CAZymes in MGS acting on substrates by diet association"
    ) +
    THEME_DIME
  
  # Plot with only significant substrates
  tests <- ggplot_build(plt_substrate)$data[[2]]
  tests$substrate <- substrate_counts |>
    filter(Substrate %in% keep_substrates) |>
    # For now manually filter those with all single copy
    filter(!(Substrate %in% filter_out)) |>
    pull(Substrate) |>
    unique() |>
    sort()
  # Manually perform adjustment at stat_compare_means does not do this
  # across facets
  tests$p.adj.manual <- p.adjust(tests$p, method = "BH")
  sig_tests <- tests |>
    filter(p <= 0.05)
  
  filtered_substrate_counts <- substrate_counts |>
    filter(Substrate %in% keep_substrates) |>
    # For now manually filter those with all single copy
    filter(!(Substrate %in%
    c("and  beta-glucan", "glucomannan", "raffinose"))) |>
    filter(Substrate %in% sig_tests$substrate) |>
    # Reorder bioactive conditions to match other plots
    ungroup() |>
    mutate(
      association = fct_relevel(
        factor(association),
        "None",
        "Low Bioactive",
        "High Bioactive"
      )
    )
  
  sig_tests_use <- sig_tests |>
    left_join(
      filtered_substrate_counts |>
      group_by(Substrate) |>
      summarise(label_pos = max(substrate_enzymes) * 1.05),
      join_by(substrate == Substrate),
    ) |>
    mutate(
      Substrate = substrate,
      p.adj.format = glue("q={round(p.adj.manual, 2)}, p={round(p, 2)}")
    )
  
  plt_sig_substrate <- filtered_substrate_counts |>
    ggplot(
      aes(
        x = factor(association),
        y = substrate_enzymes,
        fill = association,
        color = association
      )
    ) +
    geom_boxplot(alpha = 0.65, linewidth = 1) +
    facet_wrap(
      ~Substrate,
      scale="free_y",
      ncol = 8,
      labeller = labeller(Substrate = label_wrap_gen(18))
    ) +
    # ggpubr::stat_compare_means(
    #   label = "p.format",
    #   label.x.npc = "centre",
    #   label.y.npc = "top",
    #   show.legend = FALSE
    # ) +
    geom_text(
      data = sig_tests_use,
      aes(
        x = "Low Bioactive",
        y = label_pos,
        label = p.adj.format,
        color = "black",
        fill = "black"
      ),
      color = "black"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    scale_x_discrete(limits = c("None", "Low Bioactive", "High Bioactive")) +
    scale_fill_manual(
      values = c(
        "High Bioactive" = unname(BIOACTIVE_COLORS[2]),
        "Low Bioactive" = unname(BIOACTIVE_COLORS[1]),
        "None" = noassoc_color
      ),
      name = "Diet Association"
    ) +
    scale_color_manual(
      values = c(
        "High Bioactive" = unname(BIOACTIVE_COLORS[2]),
        "Low Bioactive" = unname(BIOACTIVE_COLORS[1]),
        "None" = noassoc_color
      ),
      name = "Diet Association"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab("") +
    ylab("Number of CAZymes in MGS") +
    ggtitle(
      "Number of CAZymes in MGS acting on substrates by diet association"
    ) +
    THEME_DIME

  return(list(all = plt_substrate, sig = plt_sig_substrate))
}

#' Clean and add matching sample identifiers for each sample
#'
#' @param tbl_bioactives_unadjusted Bioactives before adjustment as tibble
#' @returns Tibble with bioactive columns, and sample_id
clean_unadjusted_bioactives <- function(
  tbl_bioactives_unadjusted
) {
  tbl_cleaned <- tbl_bioactives_unadjusted |>
    mutate(
      participant_id = Participants |>
        str_sub(-2, -1) |>
        str_replace("O", "0"),
      sample_id = paste0(participant_id, "-", time_point)
    ) |>
    relocate(sample_id, participant_id) |>
    # Select down to only bioactive columns
    select(sample_id, Lignans:NEPP) |>
    # Rename a typoed bioactive
    rename(Capsaicinoids = Capasicosoids)
  return(tbl_cleaned)
}