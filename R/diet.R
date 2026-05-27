metabolite_food_correlation <- function(
  tbl_metab_peaks,
  tbl_metab_ttest,
  tbl_metab_annotations,
  lst_food_freqs,
  diet_condition
) {
  #' Looks at only those metabolites which are different, and increased in
  #' the diet of interest
  mb_xls <- tbl_metab_ttest
  peaks <- tbl_metab_peaks
  ident <- tbl_metab_annotations

  # Select only metabolites higher in this condition
  sig_metab_long <- mb_xls |>
    mutate(
      qval = `Qvalue ttest`,
      log2fc = `Individual Estimate`,
      enr_in = ifelse(qval < 0.05,
                          ifelse(log2fc > 0, "High Bioactive", "Low Bioactive"),
                          "None")
    ) |>
    filter(qval <= 0.05)
  sig_metab_diet <- sig_metab_long |>
    filter(enr_in == diet_condition)
  sig_metab_diet_wide <- peaks |>
    filter(full_id %in% sig_metab_diet$full_id) |>
    filter(!is.na(dime_name)) |>
    # Want just participant name and HB/LB
    mutate(
      participant = map(sample_name, \(x) {
        unlist(strsplit(x, "-", fixed = TRUE))[[1]]
      }),
      diet = map(sample_arm, \(x) {
        if(x == "after_high") {
          return("High Bioactive")
        }
        if(x == "after_low") {
          return("Low Bioactive")
        }
        return("")
      }),
      simple_id = paste0(participant, "_", diet)
    ) |>
    filter(diet == diet_condition) |>
    pivot_wider(id_cols = simple_id,
                names_from = full_id,
                values_from = into_driftcor_batchcor_pqn,
                values_fill = 0)
  grp_data <- lst_food_freqs[[diet_condition]] |>
    mutate(simple_id = paste0(Participant, "_", diet_condition)) |>
    select(-Participant) |>
    column_to_rownames("simple_id")

  # Map over combinations of metabolite / food group
  diet_combinations <- expand.grid(
    food_group = grp_data |> colnames(),
    metabolites = colnames(sig_metab_diet_wide |> select(-simple_id))
  )
  diet_tests <- map2(diet_combinations$food_group, diet_combinations$metabolites, \(fdgrp, metb) {
    # Get data and ensure in correct order
    met_vec <- sig_metab_diet_wide[, metb |> as.vector()] |> as.vector() |> unlist()
    food_vec <- grp_data[sig_metab_diet_wide$simple_id, fdgrp |> as.vector()] |> as.vector()
    test_res <- cor.test(x=met_vec, y=food_vec, 
                        alternative = "greater",
                        method = "spearman")
    test_res$food_group <- fdgrp
    test_res$metabolite <- metb
    return(test_res)
  })
  diet_test_df <- do.call(rbind, diet_tests) |> as.data.frame()
  library(mutoss)
  step_down_q <- multiple.down(diet_test_df$p.value |> unlist(), alpha = 0.05)
  diet_test_df$q_val_sd <- step_down_q$adjPValues
  # Compare to a standard BH
  diet_test_df$q_val_bh <- p.adjust(
    diet_test_df$p.value|> unlist(), method = "BH")
  diet_test_df$sig_char <- map2_chr(diet_test_df$p.value, diet_test_df$q_val_bh,
                                \(p, q) {
                                  if(q <= 0.01) {
                                    return ("***")
                                  }
                                  if(q <= 0.05) {
                                    return ("**")
                                  }
                                  if(q <= 0.1) {
                                    return ("*")
                                  }
                                  if(p <= 0.01) {
                                    return("..")
                                  }
                                  if(p <= 0.05) {
                                    return(".")
                                  }
                                  return("")
                                })
  return(diet_test_df)
}

plot_metabolite_food_correlation <- function(
  tbl_fm_corr,
  alpha,
  ...
) {
  library(pheatmap)
  filt <- tbl_fm_corr |>
    mutate(estimate = estimate |> unlist(),
          food_group = food_group |> unlist(),
          metabolite = metabolite |> unlist()) |>
    # Filter to only metabolites with a p =< 0.05
    dplyr::group_by(metabolite) |>
    filter(min(p.value |> unlist()) <= alpha) |>
    ungroup()
  filt_r <- filt |>
    pivot_wider(id_cols = food_group,
                names_from = metabolite,
                values_from = estimate) |> 
    column_to_rownames("food_group")
  filt_pchar <- filt |>
    mutate(estimate = estimate |> unlist(),
            food_group = food_group |> unlist(),
            metabolite = metabolite |> unlist()) |>
    pivot_wider(id_cols = food_group,
                names_from = metabolite,
                values_from = sig_char) |>
    column_to_rownames("food_group")
  # Annotation for whether metabolites are urine or faeces
  col_annotation <- data.frame(
    row.names = colnames(filt_r),
    Source = ifelse(grepl("urine", colnames(filt_r)), "Urine", "Faeces")
  )
  filt_r |>
    pheatmap(cluster_cols = TRUE, cluster_rows = FALSE,
            display_numbers = filt_pchar,
            width = 10,
            height = 4,
            treeheight_col = 0,
            breaks = seq(-1, 1, length.out=101),
            annotation_col = col_annotation,
            ...
            )
}

distance_nutrient <- function(
  df_nutrient
) {
  #' Euclidean distance between nutrient intakes
  #'
  #' Data is centred and scaled. No filtering of nutrient is performed.
  #' @param df_nutrient Nutrient intake as a dataframe, samples on columns
  #' @returns Distance matrix from vegdist
  df_nutrient |>
    as.matrix() |>
    t() |>
    scale() |>
    vegdist(method = "euclidean")
}

bioactive_summary <- function(
  mat_bioactive,
  bioactive_unadjusted,
  tbl_sample_metadata
) {
  #' Calculate some summary statistics for bioactives

  # Mean intake of every bioactive in each arm - residual adjusted
  long_bioactive <- mat_bioactive |>
    as.data.frame() |>
    rownames_to_column("bioactive") |>
    pivot_longer(
      -bioactive,
      names_to = "sample",
      values_to = "value"
    ) |>
    dplyr::left_join(
      tbl_sample_metadata,
      join_by(sample == sample_id)
    ) |>
    # Remove participant 09 due to unusual dietary reports
    filter(
      participant != "09"
    )
  bioactive_arm_mean <- long_bioactive |>
    group_by(
      bioactive, sample_arm
    ) |>
    summarise(
      mean = mean(value)
    )

  # Count number of participants whose intake is higher in HB
  hb_inc <- long_bioactive |>
    filter(time_point %in% c("V1", "V2", "V4")) |>
    mutate(diet = map_chr(sample_arm, \(x) arm_to_diet(x))) |>
    select(bioactive, participant, diet, value) |>
    pivot_wider(
      id_cols = c(participant, bioactive),
      names_from = diet,
      values_from = value
    ) |>
    mutate(
      HB_LB = `High Bioactive` - `Low Bioactive`,
      HB_increase = HB_LB > 0
    )
  hb_inc_bioactive <- hb_inc |>
    group_by(bioactive) |>
    summarise(
      n_increased = sum(HB_increase),
      mean_increase = mean(HB_LB)
    ) |>
    arrange(mean_increase)
  hb_inc_participant <- hb_inc |>
    group_by(participant) |>
    summarise(n_increased = sum(HB_increase)) |>
    arrange(n_increased)

  # Summarise change per individual LB-HB, BLN-HB, BLN-LB
  # To make change ratio calculation simpler, use unadjusted values
  long_unadjusted <- bioactive_unadjusted |>
    pivot_longer(
      -sample_id,
      names_to = "bioactive",
      values_to = "value"
    ) |>
    dplyr::left_join(
      tbl_sample_metadata,
      join_by(sample_id == sample_id)
    ) |>
    # Filter participant 09 due to unusual dietary reports
    filter(
      participant != "09"
    )
  log2fc_bioactives <- long_unadjusted |>
    filter(time_point %in% c("V1", "V2", "V4")) |>
    mutate(diet = map_chr(sample_arm, \(x) arm_to_diet(x))) |>
    # Where residual adjusted value is negative, treat this as 0 for LFC
    filter(diet %in% c("High Bioactive", "Low Bioactive", "Baseline")) |>
    pivot_wider(
      names_from = diet,
      values_from = value,
      id_cols = c(bioactive, participant)
    ) |>
    mutate(
      HB_LB = log2(`High Bioactive` / `Low Bioactive`),
      HB_LB = replace(HB_LB, is.infinite(HB_LB), NA),
      HB_BLN = log2(`High Bioactive` / `Baseline`),
      HB_BLN = replace(HB_BLN, is.infinite(HB_BLN), NA),
      LB_BLN = log2(`Low Bioactive` / Baseline),
      LB_BLN = replace(LB_BLN, is.infinite(LB_BLN), NA)
    )

  limit <- max(abs(log2fc_bioactives$HB_LB), na.rm = TRUE) * c(-1, 1)
  # Order bioactives by mean Log2FC
  hblb_order <- log2fc_bioactives |>
    group_by(bioactive) |>
    summarise(mean = mean(HB_LB, na.rm = TRUE)) |>
    arrange(mean) |>
    pull(bioactive)
  plt_hblb_log2fc <- log2fc_bioactives |>
    filter(!is.na(HB_LB)) |>
    mutate(
      raised_in = ifelse(HB_LB > 0, "High Bioactive", "Low Bioactive"),
      bioactive = fct_relevel(bioactive, hblb_order)
    ) |>
    ggplot() +
    geom_line(
      aes(x=bioactive, y=HB_LB, group=participant),
      color="grey", size=0.1
    ) +
    geom_point(
      aes(y=HB_LB, x=bioactive, color=raised_in)
    ) +
    scale_color_manual(
        name="Increased In",
        values=BIOACTIVE_COLORS
    ) +
    coord_flip() +
    ylab("Log2 Fold Change (HB/LB)") +
    xlab("Bioactive") +
    THEME_DIME +
    theme(legend.position = "bottom")

  plt_hbbln_log2fc <- log2fc_bioactives |>
    filter(!is.na(HB_BLN)) |>
    mutate(
      raised_in = ifelse(HB_BLN > 0, "High Bioactive", "Baseline"),
      bioactive = fct_relevel(bioactive, hblb_order)
    ) |>
    ggplot() +
    geom_line(
      aes(x=bioactive, y=HB_BLN, group=participant),
      color="grey", size=0.1
    ) +
    geom_point(
      aes(y=HB_BLN, x=bioactive, color=raised_in)
    ) +
    scale_color_manual(
        name="Increased In",
        values=BIOACTIVE_COLORS3
    ) +
    coord_flip() +
    ylab("Log2 Fold Change (HB/BLN)") +
    xlab("Bioactive") +
    THEME_DIME +
    theme(legend.position = "bottom")

  plt_lbbln_log2fc <- log2fc_bioactives |>
    filter(!is.na(LB_BLN)) |>
    mutate(
      raised_in = ifelse(LB_BLN > 0, "Low Bioactive", "Baseline"),
      bioactive = fct_relevel(bioactive, hblb_order)
    ) |>
    ggplot() +
    geom_line(
      aes(x=bioactive, y=LB_BLN, group=participant),
      color="grey", size=0.1
    ) +
    geom_point(
      aes(y=LB_BLN, x=bioactive, color=raised_in)
    ) +
    scale_color_manual(
        name="Increased In",
        values=BIOACTIVE_COLORS3
    ) +
    coord_flip() +
    ylab("Log2 Fold Change (LB/BLN)") +
    xlab("Bioactive") +
    THEME_DIME +
    theme(legend.position = "bottom")

  return(list(
    plots=list(
      log2fc_hb_lb = plt_hblb_log2fc,
      log2fc_hb_bln = plt_hbbln_log2fc,
      log2fc_lb_bln = plt_lbbln_log2fc,
      log2fc_combined = plt_hblb_log2fc + plt_hbbln_log2fc + plt_lbbln_log2fc
    ),
    adjusted_arm_mean = bioactive_arm_mean
  ))
}

#' Perform PCA on nutrient or bioactive adjusted compositions.
#' The purpose of this is to provide a multivariate look at separation of
#' baseline, low, and high bioactive diets.
diet_ordination <- function(
  diet_composition_adj,
  tbl_sample_metadata
) {
  #' Participant 09 is removed from comparison of dietary data
  #' The diet based on their diaries is highly dissimilar to others, but their
  #' microbiome and other values appear normal - we suspect this is due to
  #' issues with dietary recording.
  tbl_filtered <- diet_composition_adj[
    , !grepl("09", colnames(diet_composition_adj))] |>
    t()

  #' Version of the table with metadata attached
  tbl_filtered_md <- tbl_filtered |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    left_join(tbl_sample_metadata, by = "sample_id")

  # PCA and extract results
  pca_res <- prcomp(tbl_filtered, scale. = TRUE, center = TRUE)
  pca_sum <- summary(pca_res)
  pca_vexp <- scales::percent(pca_sum$importance[2, 1:2], accuracy = 0.1)
  pca_df <- data.frame(pca_res$x[ ,1:3]) |>
    rownames_to_column("sample_id") |>
    left_join(tbl_sample_metadata, by = "sample_id") |>
    mutate(Diet = sample_arm |> map_chr(arm_to_diet))

  # Make styled figure
  fig <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diet)) +
    scale_color_manual(
      values = c("High Bioactive" = BIOACTIVE_COLORS3[1] |> unname(),
                "Low Bioactive" = BIOACTIVE_COLORS3[2] |> unname(),
                "Baseline" = BIOACTIVE_COLORS3[3] |> unname())
    ) +
    stat_ellipse() + 
    geom_line(data = pca_df %>% filter(Diet %in% c("High Bioactive", "Baseline")), 
              mapping = aes(x = PC1,
                            y = PC2,
                            color = Diet,
                            group = participant),
              color = BIOACTIVE_COLORS3[1],
              linewidth = 0.15) +
    geom_line(data = pca_df %>% filter(Diet %in% c("Low Bioactive", "Baseline")), 
              mapping = aes(x = PC1,
                            y = PC2,
                            color = Diet,
                            group = participant),
              color = BIOACTIVE_COLORS3[2],
              linewidth = 0.15) +
    geom_point() +
    xlab(glue("PC1 ({pca_vexp[1]})")) +
    ylab(glue("PC2 ({pca_vexp[2]})")) +
    # ggtitle("PCA of Bioactive Intake") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.text.y = element_blank())

  return(
    list(
      figure=fig,
      pca=pca_res,
      pca_summary=pca_sum,
      data=tbl_filtered
    )
  )
}

# Stacked bar of mean biaoctive intake between arms, and intake per participant
bioactive_bar <- function(
  bioactive_unadjusted,
  tbl_sample_metadata
) {
  long_bioactive <- bioactive_unadjusted |>
    left_join(tbl_sample_metadata, by = "sample_id") |>
    mutate(diet = sample_arm |> map_chr(arm_to_diet)) |>
    select(sample_id, Lignans:NEPP, diet, participant) |>
    # Exclude participant 09
    filter(participant != "09") |>
    # Fix some typos in table
    rename(`Cinnamic acid` = `Cinnaimc acid`) |>
    pivot_longer(Lignans:NEPP, names_to = "bioactive", values_to = "intake")

  # Mean intake for all arm (HB, LB, BLN)
  mean_per_bioactive <- long_bioactive |>
  group_by(diet, bioactive) |>
  summarise(mean = mean(intake),
            sd = sd(intake))

  # Mean total intake
  mean_total_intake <- long_bioactive |>
    group_by(diet, participant) |>
    summarise(sum = sum(intake)) |>
    group_by(diet) |>
    summarise(mean = mean(sum),
              sd = sd(sum))

  # Plot intake
  # Pick out the 8 most abundant categories, and group the rest under other
  n <- 8
  top_n <- (long_bioactive |>
    group_by(bioactive) |>
    summarise(mean = mean(intake)) |>
    arrange(desc(mean)))$bioactive[1:8]

  fig_bioactive_stack <- long_bioactive |>
    mutate(
      bioactive_limit = ifelse(bioactive %in% top_n, bioactive, "Other")) |>
    mutate(
      bioactive_limit = fct_relevel(bioactive_limit, c(top_n, "Other"))) |>
    group_by(diet, bioactive_limit) |>
    summarise(intake = mean(intake)) |>
    ggplot() +
    geom_col(aes(x = diet, y = intake, fill = bioactive_limit)) +
    scale_x_discrete(
      labels = c(
        low = "Low",
        high = "High",
        baseline = "Baseline"
      ),
      limits = c("Baseline", "Low Bioactive", "High Bioactive")
    ) +
    scale_fill_manual(
      values = as.vector(alphabet(n + 1)),
      limits = c(top_n, "Other")
    ) +
    ylab("Mean daily intake (mg)") +
    xlab(NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    ggtitle("Bioactive Intake") +
    guides(fill = guide_legend(nrow = 3))

  # Individual plots for supplementary
  # Include participant 09 in these
  fig_ind_stack <- bioactive_unadjusted |>
    left_join(tbl_sample_metadata, by = "sample_id") |>
    mutate(diet = sample_arm |> map_chr(arm_to_diet)) |>
    select(sample_id, Lignans:NEPP, diet, participant) |>
    # Fix some typos in table
    rename(`Cinnamic acid` = `Cinnaimc acid`) |>
    pivot_longer(Lignans:NEPP, names_to = "bioactive", values_to = "intake") |>
    mutate(
      bioactive_limit = ifelse(bioactive %in% top_n, bioactive, "Other")) |>
    mutate(
      bioactive_limit = fct_relevel(bioactive_limit, c(top_n, "Other"))) |>
    group_by(diet, bioactive_limit, participant) |>
    # summarise(intake = mean(intake)) |>
    ggplot() +
    geom_col(aes(x = diet, y = intake, fill = bioactive_limit)) +
    facet_wrap(~participant) +
    scale_x_discrete(
      labels = c(
        low = "Low",
        high = "High",
        baseline = "Baseline"
      ),
      limits = c("Baseline", "Low Bioactive", "High Bioactive")
    ) +
    scale_fill_manual(
      values = as.vector(alphabet(n + 1)),
      limits = c(top_n, "Other")
    ) +
    ylab("Mean daily intake (mg)") +
    xlab(NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_text(angle=90)
          ) +
    ggtitle("Bioactive Intake") +
    guides(fill = guide_legend(nrow = 3))
  
  return(list(
    plots = list(
      stack = fig_bioactive_stack,
      individual_stack = fig_ind_stack
    ),
    data = long_bioactive,
    mean_per_bioactive = mean_per_bioactive,
    mean_total_intake = mean_total_intake
  ))
}

paper_figure_one <- function(
  plt_bioactive_mean_intake,
  plt_bioactive_log2fc,
  plt_bioactive_pca,
  plt_nutrient_pca,
  loc_study_design
) {
  layout <- "ABCD"
  fig1_diet <- (plt_bioactive_mean_intake + 
     guides(fill = guide_legend(
      #  label.theme = element_text(size = rel(6.5)),
       nrow = 3,
       title = NULL,
       keywidth = unit(6, "pt"),
       keyheight = unit(3, "pt")
      )) +
     theme(legend.box.margin = margin(0, 0, 0, 0),
           legend.margin = margin(0, 0, 0, 0)
          #  axis.text.x = element_text(size = rel(0.6))
      )
  ) +
  # Per participant Log2FC
  (
    plt_bioactive_log2fc +
      theme(
        panel.spacing.y = unit(0, "cm"),
        axis.text.y = element_text(size = rel(0.6)),
        strip.text = element_text(size = rel(0.5))
      ) +
      guides(
        color = guide_legend(nrow = 2)
      )
  ) +
  # PCAs
  (
    plt_bioactive_pca +
      ggtitle("Bioactives")
  ) +
  (
    plt_nutrient_pca +
      ggtitle("Nutrients") +
      guides(colour = FALSE)
  ) +
  # Layout and annotations
  plot_layout(guides = "keep", design = layout, widths = c(1, 1, 2, 2)) + 
  plot_annotation(tag_levels = list(c("B", "C", "D", "E"))) &
  theme(legend.position = "bottom",
        legend.text = element_text(size = rel(0.7)),
        title = element_text(size = rel(0.7))
        )
  return(fig1_diet)
}