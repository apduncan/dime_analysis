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