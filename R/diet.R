FOOD_GROUPS = list(
  coffee = list(
    keep = c("coffee"),
    discard = c("cake", "biscuit", "processed milks")
  ),
  tea = list(
    keep = c("tea"),
    discard = c("burger", "raw", "rich tea", "cake")
  ),
  wine = list(
    keep = c("wine"),
    discard = c("vinegar", "Wine Gums")
  ),
  alcohol = list(
    keep = c("beer", "wine", "alcohol", "spirits"),
    discard = c("vinegar", "Wine Gums")
  ),
  fruit = list(
    keep = c("fruit"),
    discard = c("juice", "alcohol")
  ),
  chocolate = list(
    keep = c("chocolate", "cococa"),
    discard = c("non-chocolate")
  ),
  berry = list(
    keep = c("berry", "berries"),
    discard = c("cake", "ice", "Chilled Desserts")
  ),
  bread = list(
    keep = c("bread"),
    discard = c("Doritos", "cool original", "cakes", "meat")
  ),
  pasta = list(
    keep = c("spaghetti", "pasta"),
    discard = c("NO DISCARD")
  )
  # vegetables_general = list(
  #   keep = c("Vegetables-General"),
  #   discard = c("NO DISCARD")
  # ),
  # vegetables_roots = list(
  #   keep = c("Roots | Tubers | Bulbs"),
  #   discard = c("NO DISCARD")
  # ),
  # vegetables_leafy = list(
  #   keep = c("Leafy Vegetables"),
  #   discard = c("NO DISCARD")
  # ),
  # vegetables_peas = list(
  #   keep = c("Beans | Peas | Lentils"),
  #   discard = c("NO DISCARD")
  # ),
  # nuts_seeds = list(
  #   keep = c("Nuts & Seeds"),
  #   discard = c("NO DISCARD")
  # )
)

food_group_frequency <- function(
  tbl_food_groups,
  lst_groups
) {
  #' Frequency of food group intake
  #'
  #' Each group is a list with 'keep' as a vector of terms to keep, and
  #' 'discard' a vector of terms to discard. Returns a list with frequency
  #' during HB, LB, and BLN dietary periods.
  #'
  #' @param tbl_food_groups Itemised food diaries as tibble
  #' @param lst_groups List of food groups, with name being name of group
  #' @returns List with `High Bioactive`, `Low Bioactive`, and `Baseline`
  #' tibbles of food intake frequency.

  full_arms <- tbl_food_groups
  groups <- FOOD_GROUPS

  plt_list <- list()
  i <- 1
  for(name in names(groups)) {
    print(name)
    group <- groups[[name]]
    # Entries to keep
    fd_grps <- do.call(rbind,
                      group$keep |> map(\(x) {
                        full_arms %>%filter(
                          grepl(x, SRC_DB_Category, ignore.case = TRUE) |
                          grepl(x, Food_Category, ignore.case = TRUE) |
                          grepl(x, Food_Description, ignore.case = TRUE) |
                          grepl(x, Food_Name, ignore.case = TRUE)
                        )
                      }))
    # Discard things we don't want
    for(x in group$discard) {
      fd_grps <- fd_grps |> filter(!(
        grepl(x, SRC_DB_Category, ignore.case = TRUE) |
          grepl(x, Food_Category, ignore.case = TRUE) |
          grepl(x, Food_Description, ignore.case = TRUE) |
          grepl(x, Food_Name, ignore.case = TRUE)
      ))
    }
    # Count
    fd_sum <- fd_grps %>%
      group_by(Participant, diet) %>%
      summarise(count = n()) %>%
      # Make complete
      ungroup() %>%
      complete(Participant, diet, fill=list(count = 0))
    # Divide by number of days in log period to get a crude fruit / day measure
    log_days <- full_arms %>%
      select(Participant, diet, Day) %>%
      distinct() %>%
      group_by(Participant, diet) %>%
      summarise(log_duration = n())
    fd_prop <- fd_sum %>%
      left_join(log_days, by = join_by(Participant, diet)) %>%
      mutate(count_per_day = count / log_duration)
    # Store the results for this food group
    fd_prop$group <- name
    plt_list[[name]] <- list()
    plt_list[[name]]$name <- name
    plt_list[[name]]$frequency_df <- fd_prop
    for (field in c("SRC_DB_Category", "Food_Category", "Food_Description",
                   "Food_Name")) {
      plt_list[[name]][[field]] <- table(fd_grps[[field]])
    }
  }

  group_freq_list <- list()
  # Combine to dataframes for LB, HB, BLN
  for (arm in c("High Bioactive", "Low Bioactive", "Baseline")) {
    all_grp_df <- plt_list |> map(\(x) {
      x$frequency_df |>
        filter(diet == arm)
    })
      all_grp_df <- do.call(rbind, all_grp_df)
    # Pivot to make participant row, and group the column
      all_grp_df <- all_grp_df |>
        pivot_wider(
          id_cols = Participant,
          names_from = group,
          values_from = count_per_day,
          values_fill = 0
        )
      group_freq_list[[arm]] <- all_grp_df
    }
  return(group_freq_list)
}

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