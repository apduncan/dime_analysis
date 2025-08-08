#' Functions for metabolite analysis

metabolite_dbrda <- function(
  mat_metab,
  tbl_md,
  participant_rm = c()
) {
  #' Performed unconstrained dbRDA analysis of metabolite data
  #'
  #' Metabolite abundances are centred and scaled, and dbRDA carried out using
  #' capscale. dbRDA is conditioned by participant, and uses euclidean distance.
  #' Sequential perMANOVA is carried out (by = "margin") using participant as
  #' strata.
  #'
  #' @param mat_metab Metabolite matrix
  #' @param tbl_md Sample metadata as a tibble
  #' @param participant_rm Participants to remove (outliers etc)
  #' @return List containg dbRDA, perMANOVA results, and plot. Data used 
  #' as inputs is in the 'data' item; results in 'plot', 'dbrda', and 
  #' 'permanova'
  # Load the metabolite data
  set.seed(4298)
  metabo_peaks <- mat_metab |>
    t() |>
    scale()
  # Filter to remove any participants requested
  sample_filt <- tbl_md |>
    filter(!(participant %in% participant_rm))

  # Restrict to those samples which are in the passed metadata
  metabo_peaks <- metabo_peaks[sample_filt$sample_id, ]

  # Perform dbRDA and perMANOVA test
  metabo_rda <- capscale(
    metabo_peaks ~ Condition(participant),
    sample_filt,
    dist = "euclidean"
  )
  metabo_summary <- summary(metabo_rda)
  metab_permanova <- adonis2(
    metabo_peaks ~ participant + diet,
    sample_filt,
    method = "euclidean",
    by = "margin",
    strata = sample_filt$participant
  )

  # Fit features to ordination
  # fit <- envfit(metabo_rda, metabo_peaks)

  # Calculate some summaries of variance
  conditioned_inertia <- metabo_summary$partial.chi / metabo_summary$tot.chi
  postcond_inertia <- metabo_summary$unconst.chi
  mds_props <- (
    metabo_summary$cont$importance |>
    as.data.frame()
  )["Proportion Explained", ] |> unlist() |> as.vector()
  mds_props_pc <- mds_props * 100
  diet_pval <- metab_permanova$`Pr(>F)`[[2]]

  # Plotting
  df_met_dbrda <- metabo_summary$sites |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    select(sample_id, MDS1, MDS2) |>
    left_join(sample_filt) 
  # Determine x and y positions for labels
  x_lab <- max(df_met_dbrda$MDS1) * 1.05
  y_lab <- max(df_met_dbrda$MDS2) * 1.05
  plt_met_dbrda <- df_met_dbrda |>
    ggplot(aes(x = MDS1, y = MDS2, color = Diet)) +
    scale_color_manual(values = c(BIOACTIVE_COLORS, c(Baseline = "darkgrey"))) +
    xlab(glue("MDS1 ({round(mds_props_pc[1], 2)}%)")) +
    ylab(glue("MDS2 ({round(mds_props_pc[2], 2)}%)")) +
    stat_ellipse() +
    annotate("text", 
      label = glue(
        "Interindividual variance = {round(conditioned_inertia*100, 2)}%\nDiet P(>F) = {diet_pval}"
      ), 
      x = x_lab, y = y_lab,
      size = 2,
      vjust = "top",
      hjust = "right"
    ) +
    geom_point() +
    guides(color = guide_legend(position = "bottom"))
  return(list(
    data_used = list(
      metabolites = metabo_peaks,
      metadata = sample_filt,
      plot = df_met_dbrda
    ),
    plot = plt_met_dbrda,
    dbrda = metabo_rda,
    permanova = metab_permanova
  ))
}

metabolite_volcano <- function(
  tbl_metab_ttest,
  tbl_metab_peaks,
  tbl_metab_annotations
) {
  #' Volcano plot summarising post/post diet t-tests
  #'
  #' Analysis of metabolite levels carried out by George Savva. This function
  #' visualises results and displays using the same visual style as other
  #' plots in the paper.
  #'
  #' @param tbl_metab_ttest T-test results as tibble
  #' @param tbl_metab_peaks Peak values as a tibble
  #' @param tbl_metab_annotations Peak labels as a tibble
  #' @return Volcano plot as ggplot2 object

  mb_xls <- tbl_metab_ttest
  peaks <- tbl_metab_peaks
  ident <- tbl_metab_annotations

  ident <- ident |>
    mutate(feature_id = paste0(project, feature_id, ionmode)) |>
    group_by(feature_id) |>
    slice_max(score, n = 1) |>
    select(feature_id, plot_manual_label)

  plt_mtb_volcano <- mb_xls |>
    mutate(feature_id = paste0(Project, `Feature ID`, `Ion mode`)) |>
    select(feature_id, `Individual Estimate`, `Qvalue ttest`, Project,
           `Pvalue ttest`) |>
    left_join(ident) |>
    rename(log2fc = `Individual Estimate`, qval = `Qvalue ttest`,
           pval = `Pvalue ttest`) |>
    mutate(
      enr_in = ifelse(qval < 0.05,
                          ifelse(log2fc > 0, "High Bioactive", "Low Bioactive"),
                          "None"),
      plot_manual_label = ifelse(qval < 0.05, plot_manual_label, NA),
      Project = str_to_title(Project),
      Project = ifelse(Project == "Fecal", "Faecal", "Urine")
    ) |>
    ggplot(aes(x = (log2fc), y = -(log(pval)))) +
    geom_point(aes(color = enr_in), size = 1) +
    scale_color_manual(
      values = c(
        "High Bioactive" = "#00BA38",
        "Low Bioactive" = "#619CFF",
        "q < 0.05" = "black"
      ),
      name = "Diet") +
    geom_text_repel(
      aes(label = plot_manual_label),
      size = 2.5,
      min.segment.length = 0.2
    ) +
    facet_wrap(~Project) +
    xlab("Log2 Fold Change") +
    ylab("-log(p value)") +
    theme_pubclean()
  return(plt_mtb_volcano)
}

metabolite_figure <- function(
  plt_volcano,
  plt_dbrda_fecal,
  plt_dbrda_urine,
  plt_correlations
) {
  #' Combine metabolite subfigures
  #'
  #' Figure is vertically stacked, with A) Volcano, B) Ordination,
  #' C) Correlations
  #'
  #' @param plt_volcano Volcano plot
  #' @param plt_dbrda_fecal Ordination for faeces
  #' @param plt_dbrda_urine Ordination for urine
  #' @param plt_correlations Correlation of metabolites and food groups

  plt <- plt_volcano +
    (
      plt_dbrda_fecal +
      plt_dbrda_urine +
      plot_layout(guides = "collect") &
      THEME_DIME &
      theme(legend.position = "top")
    ) +
    as.ggplot(plt_correlations) +
    plot_layout(
      design = "A
      B
      C",
      heights = c(3, 3, 4)
    ) +
    plot_annotation(tag_levels = c("A", "B", "", "C"))
  return(plt)
}

metabolite_distance <- function(
  mat_metabolite
) {
  #' Euclidean distance based on centred, scaled data
  #' 
  #' @param tbl_metabolite Metabolite peaks as matrix, samples on columns
  #' @returns Distance matrix produced by vegdist
  mat_metabolite |>
    t() |>
    scale() |>
    vegdist(method = "euclidean")
}

metabolite_labelled_supplementary <- function(
  tbl_metab_annotations,
  output_path
) {
  #' Produce a supplementary metabolite annotation table
  #'
  #' Saves the metabolite annotations in a format suitable for inclusion as a
  #' supplementary table. Outputs as an Excel sheet.
  #'
  #' @param tbl_metab_annoations
  #' @param output_path Location to save in Excel format
  #' @returns Path to Excel sheet
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  tbl_metab_annotations |>
    filter(!is.na(isotopes) | !is.na(adduct) | !is.na(name)) |>
    select(project, ionmode, pcgroup, feature_id, isotopes, adduct, name) |>
    rename(
      Material = project,
      IonMode = ionmode,
      Group = pcgroup,
      FeatureID = feature_id,
      # Isotopes = isotopes,
      Adduct = adduct,
      Name = name
    ) |>
    # Standardise annotations where identity is unclear
    mutate(
      Name = map_chr(Name, \(x) {
        name_bits <- strsplit(x, "???", fixed = TRUE) |> unlist()
        name <- name_bits[[1]]
        if (length(name_bits) > 1) {
          name <- glue("{name}(?)")
        }
        print(name)
        return(name)
      })
    ) |>
    write.xlsx(file = output_path)
  return(output_path)
}

scfa_boxplot <- function(
  tbl_scfa,
  tbl_sample_metadata
) {
  #' Compare SCFA concentration between diets
  #'
  #' @param tbl_scfa SCFA concentration, sample on rows
  #' @param tbl_sample_metdata Metadata table
  sample_md <- tbl_sample_metadata
  fig_scfa <- tbl_scfa %>%
    rename(sample_id = ID) %>%
    left_join(sample_md) %>%
    # Restrict to only Baseline, Low, High
    filter(
      sample_arm %in% c("after_high", "after_low") | time_point == "V1") %>%
    mutate(Diet = map_chr(sample_arm, arm_to_diet)) %>%
    select(participant, Diet, Acetate:Butyrate) %>%
    pivot_longer(Acetate:Butyrate, names_to = "scfa", values_to = "amount") %>%
    arrange(participant, Diet, scfa) %>%
    ggplot(aes(x = Diet, y = amount, fill = Diet, color = Diet)) +
    geom_boxplot(alpha = 0.6, show.legend = FALSE) +
    facet_wrap(~scfa, ncol = 1, scales = "free_y") +
    scale_x_discrete(
      limits = c("Baseline", "Low Bioactive", "High Bioactive")
    ) +
    scale_y_continuous(
      expand = expand_scale(mult = c(0.1, 0.1))
    ) +
    scale_fill_manual(
      values = c("High Bioactive" = BIOACTIVE_COLORS3[2] |> unname(),
                "Low Bioactive" = BIOACTIVE_COLORS3[1] |> unname(),
                "Baseline" = BIOACTIVE_COLORS3[3] |> unname()),
      limits = c("Baseline", "Low Bioactive", "High Bioactive")
    ) +
    scale_color_manual(
      values = c("High Bioactive" = BIOACTIVE_COLORS3[2] |> unname(),
                "Low Bioactive" = BIOACTIVE_COLORS3[1] |> unname(),
                "Baseline" = BIOACTIVE_COLORS3[3] |> unname()),
      limits = c("Baseline", "Low Bioactive", "High Bioactive")
    ) +
    stat_compare_means(
      comparisons = list(c("Low Bioactive", "High Bioactive")),
      method = "wilcox.test",
      paired = TRUE,
      size = 2,
      symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
        symbols = c("****", "***", "**", "*", "ns")
      )
    ) +
    ylab(NULL) +
    xlab(NULL) +
    theme_minimal() +
    ggtitle("SCFAs") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(fig_scfa)
}

scfa_bioactive_correlation <- function(
  tbl_scfa,
  tbl_bioactive,
  tbl_sample_metadata
) {
  #' Correlation of SCFAs to bioactives
  #'
  #' @param tbl_scfa SCFA concentration, sample on rows
  #' @param tbl_bioactive Bioactive intake residual adjusted, samples on columns
  #' @param tbl_sample_metdata Sample metadata

  scfa_df <- tbl_scfa |> column_to_rownames("ID")
  bioactive_df <- tbl_bioactive |> t()
  sample_md <- tbl_sample_metadata

  within_arm_correlate <- function(arm) {
    # Reduce to only only requested diet
    sample_arm_to_correlate <- arm
    diet_samples <- sample_md |> filter(sample_arm == sample_arm_to_correlate) |> pull(sample_id)
    
    scfa_restricted <- scfa_df[diet_samples, ] |> t()
    diet_restricted <- bioactive_df[diet_samples, ] |> t()
    
    # Correlate each of the SCFAs and bioactives
    # Make a matrix with pairs to be tested
    met_rows <- length(rownames(diet_restricted))
    test_pair <- rownames(scfa_restricted) %>%
      map(\(x) {
        cbind(rep(x, met_rows), rownames(diet_restricted))
      }) %>%
      reduce(\(x, y) rbind(x, y))
    
    # Test each pair
    tests <- mapply(FUN = \(x, y) {
      # browser()
      scfa_vec   <- unlist(as.vector(scfa_restricted[x, ]))
      diet_vec   <- unlist(as.vector(diet_restricted[y, ]))
      pearson   <- cor.test(scfa_vec, diet_vec, method = "pearson")
      spearman  <- cor.test(scfa_vec, diet_vec, method = "spearman")
      c(x, y, pearson$estimate, pearson$p.value, spearman$estimate, spearman$p.value)
    }, test_pair[, 1], test_pair[, 2])
    tests <- as.data.frame(t(tests))
    colnames(tests) <- c("scfa", "bioactive", "pearson_r", "pearson_p", 
                        "spearman_rho", "spearman_p")
    tests <- tests %>% mutate(
      spearman_p   = as.numeric(spearman_p),
      spearman_rho = as.numeric(spearman_rho),
      pearson_r    = as.numeric(pearson_r),
      pearson_p    = as.numeric(pearson_p)
    )
    
    # Multiple test correction
    tests$pearson_q <- p.adjust(tests$pearson_p, method = "BH")
    tests$spearman_q <- p.adjust(tests$spearman_p, method = "BH")
    return(tests)
  }

  high_bioactive_scfa <- within_arm_correlate("after_high")
  low_bioactive_scfa <- within_arm_correlate("after_low")
  return(list(
    high_bioactive = high_bioactive_scfa,
    low_bioactive = low_bioactive_scfa
  ))
}

# scfa_bioactive_heatmap <- function(
#   lst_scfa_correlation
# ) {
#   # Heatmap of the results
#   sig_labeller <- function(pval) {
#     if(pval <= 0.01) {
#       return("**")
#     }
#     if(pval <= 0.05) {
#       return("*")
#     }
#     return("")
#   }

#   corr_hmap <- function(tests, corrtype = "pearson") {
#     corr_field <- ifelse(corrtype == "pearson", "pearson_r", "spearman_rho")
#     p_field <- ifelse(corrtype == "pearson", "pearson_p", "spearman_p")
#     corr_label <- ifelse(corrtype == "pearson", "Pearson's r", "Spearman's rho")
#     limit <- max(abs(tests[, corr_field]))
#     plot <- tests |>
#       mutate(sig_label = map_chr(!!sym(p_field), sig_labeller)) |>
#       mutate(scfa = fct_relevel(scfa, colnames(scfa_df)) |> fct_rev()) |>
#       mutate(bioactive = fct_relevel(bioactive, colnames(bioactive_df))) |>
#       ggplot(aes(y = scfa, x = bioactive, fill = !!sym(corr_field))) +
#       geom_tile() +
#       geom_text(aes(label = sig_label), color = "white") +
#       scale_fill_distiller(
#         palette = "RdYlBu",
#         limits = c(-limit, limit),
#         name = corr_label
#       ) +
#       theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#       ) +
#       xlab("Bioactive") +
#       ylab("Short Chain Fatty Acid")
#     return(plot)
#   }

#   corr_hmap(high_bioactive_scfa, corrtype = "pearson")
#   return(corr_hmap)
# }

plot_scfa_bioactive_correlation <- function(
  tbl_scfa,
  ...
) {
  library(pheatmap)

  # Significance character
  sig_labeller <- function(pval) {
    if(pval <= 0.01) {
      return("**")
    }
    if(pval <= 0.05) {
      return("*")
    }
    return("")
  }
  sig_q_labeller <- function(p, q) {
    if (q <= 0.01) {
      return("***")
    }
    if (q <= 0.05) {
      return("**")
    }
    if (q <= 0.1) {
      return("*")
    }
    if (p <= 0.01) {
      return("..")
    }
    if (p <= 0.05) {
      return(".")
    }
    return("")
  }
  tbl_scfa <- tbl_scfa |>
    mutate(
      sig_p = map_chr(pearson_p, sig_labeller),
      sig_char = map2_chr(pearson_p, pearson_q, sig_q_labeller)
    )

  # Convert long to matrix of r and p vals
  filt_r <- tbl_scfa |>
    pivot_wider(id_cols = scfa,
                names_from = bioactive,
                values_from = pearson_r) |>
    column_to_rownames("scfa")

  filt_pchar <- tbl_scfa |>
    pivot_wider(id_cols = scfa,
                names_from = bioactive,
                values_from = sig_char) |>
    column_to_rownames("scfa")

  # Annotation for whether metabolites are urine or faeces
  row_annotation <- data.frame(
    row.names =rownames(filt_r),
    Source = rep("Faeces", dim(filt_r)[[1]])
  )
  # Custom palette from -1 to 1 to be consistent
  hmap <- filt_r |>
    pheatmap(cluster_cols = FALSE,
            cluster_rows = TRUE,
            display_numbers = filt_pchar,
            breaks = seq(-1, 1, length.out = 101),
            width = 10,
            height = 4,
            treeheight_col = 0,
            treeheight_row = 0,
            # annotation_row = row_annotation,
            legend = TRUE,
            annotation_legend = FALSE,
            ...
            )
  return(hmap)
}

plot_dbrda_correlation_heatmap <- function(
  lst_fit_dbrda,
  tbl_metab_annotations
) {
  #' Heatmap of named metabolites and their correlations to dbRDA centroids
  #'
  #' @param lst_fit_dbrda Results from esdbrd_metab_corr, a list with
  #' angles between metabolite vectors fit envfit and species centroids.
  #' @param tbl_metab_annotations Table of annotations for untargetted
  #' metabolite peaks.
  #' @return ggplot tile

  # Convert cosine into long format and join to metadata
  tbl_ann <- tbl_metab_annotations |>
    mutate(id = paste0(project, feature_id, ionmode)) |>
    distinct(id, .keep_all = TRUE) |>
    filter(!is.name(name))
  plt_subset <- lst_fit_dbrda$full |>
    left_join(tbl_ann, join_by(metab_vec == id)) |>
    filter(!is.na(plot_manual_label))
  plt_tile <- plt_subset |>
    mutate(
      plotname = glue("({metab_vec}) {plot_manual_label}"),
      cos_rnd = round(cosine, 3)
    ) |>
    ggplot(aes(x = species, y = plotname, fill = cosine, label = cos_rnd)) +
    geom_tile() +
    geom_text() +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0) +
    ggtitle("Angle between ES centroid and labelled metabolite vectors") +
    xlab("Enterosignature") +
    ylab("Metabolite Peak")
  return(plt_tile)
}

paper_figure_three <- function(
  plt_volcano,
  plt_scfa_box,
  plt_scfa_correlation
) {
  #' Combine metabolite subfigures
  #'
  #' Figure is vertically stacked, with A) Volcano, B) SCFA box,
  #' C) SCFA correlations, D) UM correlations to food groups
  #'
  #' @param plt_volcano Volcano plot
  #' @param plt_scfa_box Difference in SCFA concentration HB/LB
  #' @param plt_scfa_correlation Correlation betwen nutrients and SCFA
  #' @param plt_um_correlations Correlation of metabolites and food groups

  layout <- "AAAA
  BCCC"

  plt <- wrap_plots(
    plt_volcano +
      THEME_DIME +
      labs(title = "Metabolite Differential Abundance"),
    plt_scfa_box +
      theme_minimal() +
      THEME_DIME +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
      ) +
      labs(
        title = "SCFA",
        subtitle = "Between Diets"
      ),
    as.ggplot(plt_scfa_correlation) +
      labs(
        title = "SCFA",
        subtitle = "Correlation to Bioactives"
      ),
    design = layout,
    heights = c(2, 3),
    widths = c(1, 2, 2, 1)
  ) +
    plot_annotation(tag_levels = "A")
  return(plt)
}

paper_figure_three_alternate <- function(
  plt_volcano,
  plt_scfa_box,
  plt_scfa_correlation
) {
  #' Combine metabolite subfigures, omitting food groups
  #'
  #' Figure is vertically stacked, with A) Volcano, B) SCFA box,
  #' C) SCFA correlations
  #'
  #' @param plt_volcano Volcano plot
  #' @param plt_scfa_box Difference in SCFA concentration HB/LB
  #' @param plt_scfa_correlation Correlation betwen nutrients and SCFA

  layout <- "AAAA
  BCCC"

  plt <- wrap_plots(
    plt_volcano +
      THEME_DIME +
      labs(title = "Metabolite Differential Abundance"),
    plt_scfa_box +
      theme_minimal() +
      THEME_DIME +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
      ) +
      labs(
        title = "SCFA",
        subtitle = "Between Diets"
      ),
    as.ggplot(plt_scfa_correlation) +
      labs(
        title = "SCFA",
        subtitle = "Correlation to Bioactives"
      ),
    design = layout,
    heights = c(2, 3, 3),
    widths = c(1, 2, 2, 1)
  ) +
    plot_annotation(tag_levels = "A")
  return(plt)
}

metabolite_review_list <- function(
  tbl_metab_ttest,
  tbl_metab_peaks,
  tbl_metab_annotations,
  tbl_taxa_corr
) {
  #' Select metabolites which are of interest for further identification

  # These are metabolites which form 'hubs' in the network. They are the
  # most interesting to me
  network <- c('fecal535neg', 'fecal367neg', 'fecal703neg', 'fecal1214pos')

  # This is an externally produced list of fecal metabolites which correlate
  # to microbial taxa. These are also of interest due to their relationship to
  # the microbiome
  tax_corr <- tbl_taxa_corr |>
    mutate(network_hub = metabolite %in% network) |>
    rename(estimated_taxa_correlation = corr_est) |>
    arrange(-network_hub, -abs(estimated_taxa_correlation)) |>
    # Select these out of the ttest dataframe to get LFC, Lars/Jan are
    # interested in that also
    left_join(
      tbl_metab_ttest |> select(full_id, `Individual Estimate`),
      join_by(metabolite == full_id)
    ) |>
    rename(log2fc = `Individual Estimate`)

  # Select highest LFCs which are signficant
  tbl_top_lfc <- tbl_metab_ttest |>
    filter(`Qvalue ttest` <= 0.05) |>
    slice_max(`Individual Estimate`, n = 50) |>
    arrange(-abs(`Individual Estimate`)) |>
    rename(
      log2fc = `Individual Estimate`,
      metabolite = full_id
    ) |>
    select(metabolite, log2fc) |>
    mutate(
      top_50 = TRUE
    )
  
  return(list(
    taxa = tax_corr,
    lfx = tbl_top_lfc
  ))
}

metabolite_highest_labelled <- function(
  tbl_metab_ttest,
  tbl_metab_annotations,
  material = "fecal"
) {
  ident_lbl <- tbl_metab_annotations |>
    mutate(full_id = paste0(project, feature_id, ionmode))
  tbl_top_lfc <- tbl_metab_ttest |>
    filter(Project == material) |>
    left_join(ident_lbl, join_by(full_id == full_id)) |>
    filter(`Qvalue ttest` <= 0.05) |>
    filter(!is.na(plot_manual_label)) |>
    slice_max(`Individual Estimate`, n = 50) |>
    arrange(-abs(`Individual Estimate`)) |>
    rename(
      log2fc = `Individual Estimate`,
      metabolite = full_id
    ) |>
    select(metabolite, log2fc, plot_manual_label) |>
    mutate(
      top_50 = TRUE
    )
}

format_additional_idents <- function(tbl, confidence_max) {
  #' Convert additional identifications into a format compatible with
  #' existing volcano plot code

  # Excepts to have columns score, project(f/u), ionmode (neg/pos), name,
  # plot_manual_label
  parts <- str_match(tbl$metabolite, "([a-z]*)(\\d*)([a-z]*)")
  tbl <- tbl |>
    mutate(
      project = parts[, 2],
      ionmode = parts[, 4],
      feature_id = parts[, 3],
      score = 1.0
    ) |>
    filter(plot_confidence <= confidence_max)
  return(tbl)
}

combine_additional_idents <- function(
  tbl_metab_annotations,
  tbl_metab_revised,
  confidence_max_new = 2
) {
  #' Collaborators performed additional identification for fecal metabolites.
  #' Want to keep the existing annotation of urinary metabolites.
  tbl_new_fecal <- tbl_metab_revised |>
    format_additional_idents(confidence_max = confidence_max_new) |>
    select(project, feature_id, ionmode, plot_manual_label, score)
  tbl_exi_urine <- tbl_metab_annotations |>
    filter(project == "urine") |>
    select(project, feature_id, ionmode, plot_manual_label, score) |>
    filter(!grepl("?", plot_manual_label, fixed = TRUE))
  return(rbind(tbl_new_fecal, tbl_exi_urine))
}
