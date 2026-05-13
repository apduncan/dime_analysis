# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse", "vegan", "glue", "patchwork", "ggplot2", "readxl",
    "ggrepel", "ggpubr", "pheatmap", "ggplotify", "patchwork",
    "usedist", "reticulate", "lsa", "rtk", "logger", "openxlsx",
    "colorspace", "igraph", "data.table", "igraph", "SpiecEasi", "NetCoMi",
    "tidyverse", "tools"
  )
)

tar_source()

# Target list
list(
  # ==== READ METADATA ====
  tar_target(
    pth_sample_metadata,
    "data/source/sample_metadata.csv",
    format = "file"
  ),
  tar_target(
    tbl_sample_metadata,
    read_metadata(pth_sample_metadata)
  ),

  # ==== READ METABOLITES ====
  # These are matrices of peak values in each sample, extracted from
  # peaklist_long_norm.rds
  tar_target(
    pth_metab_urine,
    "data/derived/untargetted_metabolomics/um_full_urine.tsv",
    format = "file"
  ),
  tar_target(
    mat_metab_urine,
    read_metabolites(pth_metab_urine) |>
    order_table(rownames = TRUE)
  ),
  tar_target(
    pth_metab_faeces,
    "data/derived/untargetted_metabolomics/um_full_fecal.tsv",
    format = "file"
  ),
  tar_target(
    mat_metab_faeces,
    read_metabolites(pth_metab_faeces) |>
    order_table(rownames = TRUE)
  ),
  tar_target(
    pth_metab_faeces_sig,
    "data/derived/untargetted_metabolomics/um_sig_fecal.tsv",
    format = "file"
  ),
  tar_target(
    mat_metab_faeces_sig,
    read_metabolites(pth_metab_faeces_sig) |>
      order_table(rownames = TRUE)
  ),
  # These are the results of statistical analysis carried out by George Savva
  # Code for this analysis is available in data/source/untargetted_metabolomics
  # as quarto markdown and knit versions.
  tar_target(
    pth_metab_ttest,
    "data/source/untargetted_metabolomics/test_results.xlsx"
  ),
  tar_target(
    tbl_metab_ttest,
    read_xlsx(pth_metab_ttest) |>
      mutate(full_id = paste0(Project, `Feature ID`, `Ion mode`))
  ),
  # This is output of untargetted metabolomics, which has been further processed
  # in several ways
  tar_target(
    pth_metab_peaks,
    paste0(c("data/source/untargetted_metabolomics/analysis/raw/",
      "DIME_results/save/peaklist_long_norm.rds"),
      collapse = ""
    )
  ),
  tar_target(
    tbl_metab_peaks,
    readRDS(pth_metab_peaks) |>
      mutate(full_id = paste0(gsub("DIME", "", project), feature_id, ionmode))
  ),
  # Initial tentative annotations of some urine & fecal peaks
  tar_target(
    pth_metab_annotations,
    paste0(c("data/source/untargetted_metabolomics/",
             "initial_identifications.xlsx"),
      collapse = ""
    )
  ),
  tar_target(
    tbl_metab_annotations,
    read_xlsx(pth_metab_annotations)
  ),
  # Make an output with metabolite labels
  tar_target(
    pth_metab_label_supplementary,
    metabolite_labelled_supplementary(
      tbl_metab_annotations = tbl_metab_annotations,
      output_path = "output/tables/supplementary/metabolite_labels.xlsx"
    )
  ),
  # Additional identifications for faecal metabolites of interest
  tar_target(
    pth_metab_additional_annotations,
    paste0(c("data/source/untargetted_metabolomics/",
             "additional_fecal_identifications.xlsx"),
      collapse = ""
    )
  ),
  tar_target(
    tbl_metab_additional_annotations,
    read_xlsx(pth_metab_additional_annotations) |>
      select(metabolite, estimated_taxa_correlation, network_hub, log2fc,
      `most likely (or comment)`, `certainty`, plot_manual_label,
      plot_confidence) |>
      rename(identity_or_comment = `most likely (or comment)`)
  ),

  # ==== READ SCFA ====
  tar_target(
    pth_scfa,
    "data/source/scfa/SCFA_DIME.csv",
    format = "file"
  ),
  tar_target(
    tbl_scfa,
    read_delim(pth_scfa)
  ),

  # ==== READ FUNCTIONAL ANNOTATIONS ====
  # ==== READ DBCAN ====
  tar_target(
    pth_dbcan,
    "data/source/microbiome/dbcan",
    format = "file"
  ),
  tar_target(
    lst_dbcan,
    read_dbcan(pth_dbcan)
  ),
  tar_target(
    pth_mgs_association,
    "data/source/microbiome/dbcan/s1_traits.tsv",
    format = "file"
  ),
  tar_target(
    tbl_mgs_association,
    read_delim(pth_mgs_association)
  ),

  # ==== READ DIET ====
  tar_target(
    pth_diet_composition_adj,
    "data/source/diet/residual_adjusted.tsv",
    format = "file"
  ),
  tar_target(
    pth_diet_nutrients,
    "data/source/diet/nutrient_diet.tsv",
    format = "file"
  ),
  tar_target(
    lst_diet_composition,
    read_food_composition(
      pth_diet_composition_adj,
      pth_diet_nutrients
    )
  ),
  tar_file_read(
    tbl_bioactives_unadjusted,
    file.path("data", "source", "diet", "bioactive_diet.tsv"),
    read_delim(file = !!.x)
  ),
  tar_target(
    tbl_bioactives_unadjusted_clean,
    clean_unadjusted_bioactives(tbl_bioactives_unadjusted)
  ),

  # ==== READ MICROBIOME ====
  # Microbiome
  tar_target(
    pth_genus,
    "data/source/microbiome/taxa/MGS.matL5.txt",
    format = "file"
  ),
  tar_target(
    pth_species,
    "data/source/microbiome/taxa/MGS.matL6.txt",
    format = "file"
  ),
  tar_file_read(
    tbl_mgs,
    "data/source/microbiome/taxa/MGS.matL7.txt",
    read_taxa(!!.x) |>
      order_table()
  ),
  tar_file_read(
    tbl_cell_count,
    "data/source/flow_cytometry/DIME_samples_all_1in100.csv",
    read_cellcounts(!!.x)
  ),
  tar_target(
    tbl_species,
    read_taxa(pth_species) |>
    order_table()
  ),
  tar_target(
    tbl_genus,
    read_taxa(pth_genus) |>
    order_table()
  ),
  tar_target(
    tbl_species_tss,
    taxa_tss(tbl_species) |>
    order_table()
  ),
  # Absolute microbiome abundances
  tar_target(
    tbl_rare_scaled_species,
    scale_taxa(
      lst_species_alpha$matrix |>
        data.frame(check.names = FALSE) |>
        rownames_to_column("species"),
      tbl_cell_count
    )
  ),

  # ==== CALCULATE MICROBIOME DISSIMILARITIES ====
  tar_target(
    dst_species,
    taxa_distance(tbl_species_tss)
  ),
  tar_target(
    dst_rare_species,
    taxa_distance(
      lst_species_alpha$matrix |>
        data.frame(check.names = FALSE) |>
        rownames_to_column("species")
    )
  ),
  tar_target(
    dst_rare_scaled_species,
    taxa_distance_abs(
      tbl_rare_scaled_species
    )
  ),

  # ==== CALCULATE ENTEROSIGNATURES ====
  tar_target(
    pth_es,
    enterosignature_reapply(
      pth_genus
    )
  ),
  tar_target(
    lst_es,
    read_es(pth_es)
  ),

  # ==== ANALYSIS OF ENTEROSIGNATURES ====
  tar_target(
    dbrda_es,
    enterosignature_dbrda(
      lst_es$h,
      lst_es$w,
      tbl_sample_metadata
    )
  ),
  tar_target(
    lst_metab_es_dbrda_fit,
    esdbrda_metab_corr(
      dbrda_es,
      mat_metab_faeces,
      tbl_sample_metadata
    )
  ),

  # ==== ANALYSIS OF TAXONOMIC DIVERSITY ====
  # Alpha diversity
  tar_target(
    lst_species_alpha,
    rarefy_tbl(
      tbl_species,
      repeats = 100,
      threads = 6,
      seed = 7754322456789
    )
  ),
  tar_target(
    plt_species_shannon,
    plot_alpha_diversity(
      lst_species_alpha$alpha_diversity |> select(shannon),
      tbl_sample_metadata = tbl_sample_metadata
    )
  ),
  tar_target(
    plt_species_ad_all,
    plot_alpha_diversity(
      lst_species_alpha$alpha_diversity,
      tbl_sample_metadata = tbl_sample_metadata
    )
  ),
  tar_target(
    pth_species_ad_all,
    write_figure(
      plt = plt_species_ad_all,
      pth = "supp_alpha_div_all",
      width = 5,
      height = 3,
      scale = 1
    )
  ),
  tar_target(
    tbl_summarise_ad,
    summarise_alpha_div(
      lst_species_alpha,
      tbl_sample_metadata
    )
  ),
  tar_target(
    pth_summarise_ad,
    write_table(
      tbl_summarise_ad,
      "summarise_alpha_diversity"
    )
  ),

  # Beta diversity
  tar_target(
    plt_rare_scaled_species_bc_within,
    plot_beta_within(
      dst_rare_scaled_species,
      tbl_sample_metadata,
      only_baseline = FALSE
    )
  ),
  tar_target(
    plt_mgs_substrate,
    substrate_difference(
      lst_dbcan$counts_long,
      tbl_mgs_association
    )
  ),
  tar_target(
    tbl_summarise_bd,
    summarise_beta_within(
      dst_rare_scaled_species,
      tbl_sample_metadata
    )
  ),
  tar_target(
    pth_summarise_bd,
    write_table(
      tbl_summarise_bd,
      "summarise_beta_diversity"
    )
  ),

  # ==== READ FUNCTION ====
  tar_target(
    pth_ko,
    "data/source/microbiome/function/KGML0.txt",
    format = "file"
  ),
  tar_target(
    tbl_ko,
    read_ko(pth_ko)
  ),
  tar_target(
    lst_ko_alpha,
    rarefy_tbl(
      tbl_ko |> filter_unknown() |> prevalence_filter(5),
      repeats = 100,
      threads = 6,
      seed = 9898
    )
  ),
  # PFAM annotation (from eggnogMapper)
  tar_target(
    pth_pfam,
    "data/source/microbiome/function/EM.PFAML0.tsv",
    format = "file"
  ),
  tar_target(
    tbl_pfam,
    read_pfam(pth_pfam)
  ),

  # ==== ANALYSIS OF FUNCTION DIVERSITY ====
  tar_target(
    tbl_tss_ko,
    taxa_tss(
      lst_ko_alpha$matrix |>
        as_tibble() |>
        rownames_to_column("KO")
    )
  ),
  tar_target(
    lst_ko_permanova,
    permanova_abundance(
      tbl_abd = tbl_tss_ko,
      tbl_md = tbl_sample_metadata |>
        remove_before_metadata(),
      fm_terms = c("participant", "sample_arm"),
      by = "terms",
      seed = 9898,
      method = "bray",
      permutations = 10000
    )
  ),
  tar_target(
    plt_ko_shannon,
    plot_alpha_diversity(
      lst_ko_alpha$alpha_diversity |> select(shannon),
      tbl_sample_metadata
    )
  ),
  tar_target(
    plt_ko_ad_all,
    plot_alpha_diversity(
      lst_ko_alpha$alpha_diversity,
      tbl_sample_metadata
    )
  ),
  tar_target(
    tbl_scale_ko,
    scale_taxa(
      lst_ko_alpha$matrix |> as.data.frame() |> rownames_to_column("KO"),
      tbl_cell_count
    )
  ),
  tar_target(
    dst_scale_ko,
    taxa_distance_abs(tbl_scale_ko)
  ),
  tar_target(
    plt_ko_bray,
    plot_beta_within(
      dst_scale_ko,
      tbl_sample_metadata
    )
  ),
  tar_target(
    lst_species_dbrda,
    diet_dbrda(
      tbl_species_tss |> column_to_rownames("L6") |> as.matrix() |> t() |>
        as.data.frame(check.names = FALSE),
      metadata = tbl_sample_metadata,
      centroids = TRUE,
      condition = TRUE,
      n_species_vectors = 0
    )
  ),
  tar_target(
    plt_spec_dbrda,
    decorate_species_dbrda(
      lst_species_dbrda
    )
  ),
  tar_target(
    pcoa_species_tss,
    pcoa(tbl_species_tss |> remove_before_samples(tbl_sample_metadata))
  ),
  tar_target(
    plt_species_pcoa,
    plot_pcoa(
      pcoa_species_tss,
      tbl_sample_metadata
    )
  ),

  # ==== ANALYSIS OF METABOLITES ====
  tar_target(
    dbrda_urine,
    metabolite_dbrda(
      mat_metab = mat_metab_urine,
      tbl_md = remove_before_metadata(tbl_sample_metadata),
      participant_rm = c("07") # Outlier
    )
  ),
  tar_target(
    dbrda_faeces,
    metabolite_dbrda(
      mat_metab = mat_metab_faeces,
      tbl_md = remove_before_metadata(tbl_sample_metadata),
      participant_rm = c()
    )
  ),
  tar_target(
    plt_metab_volcano,
    metabolite_volcano(
      tbl_metab_ttest,
      tbl_metab_peaks,
      tbl_metab_annotations
    )
  ),
  tar_target(
    plt_metab_volcano_revised_ident,
    metabolite_volcano(
      tbl_metab_ttest,
      tbl_metab_peaks,
      combine_additional_idents(
        tbl_metab_annotations,
        tbl_metab_additional_annotations
      )
    )
  ),
  tar_target(
    plt_metab_figure,
    metabolite_figure(
      plt_metab_volcano,
      dbrda_faeces$plot,
      dbrda_urine$plot,
      plt_food_metabolite_corr_hb
    )
  ),
  tar_target(
    pth_metab_figure,
    write_figure(
      plt_metab_figure,
      "metabolites/metabolite_figure",
      height = 8.5,
      width = 7
    )
  ),
  tar_target(
    plt_scfa_box,
    scfa_boxplot(
      tbl_scfa,
      tbl_sample_metadata |> remove_before_metadata()
    )
  ),
  tar_target(
    plt_dbrda_metab_assoc,
    plot_dbrda_correlation_heatmap(
      lst_metab_es_dbrda_fit,
      tbl_metab_annotations
    )
  ),
  tar_target(
    pth_dbrda_metab_assoc,
    write_figure(
      plt_dbrda_metab_assoc,
      "metabolites/supp_es_metab_assoc_hmap",
      height = 3,
      width = 7,
      scale = 2
    )
  ),

  # ==== CALCULATE FOOD AND METABOLITE DISSIMILARITIES ====
  tar_target(
    dst_urine,
    metabolite_distance(mat_metab_urine)
  ),
  tar_target(
    dst_faeces,
    metabolite_distance(mat_metab_faeces)
  ),
  tar_target(
    dst_diet_all,
    distance_nutrient(lst_diet_composition$all)
  ),
  tar_target(
    dst_nutrient,
    distance_nutrient(lst_diet_composition$nutrient)
  ),
  tar_target(
    dst_bioactive,
    distance_nutrient(lst_diet_composition$bioactive)
  ),


  # ==== ANALYSIS FOOD AND METABOLITES ====
  tar_target(
    lst_scfa_correlation,
    scfa_bioactive_correlation(
      tbl_scfa,
      lst_diet_composition$bioactive,
      tbl_sample_metadata
    )
  ),
  tar_target(
    plt_scfa_correlation,
    plot_scfa_bioactive_correlation(
      lst_scfa_correlation$high_bioactive,
      fontsize_row = 6,
      fontsize_col = 6,
      fontsize = 6,
      fontsize_number = 6,
      angle_col = 45
    )
  ),

  # ==== ANALYSIS FOOD METABOLITE TAXA ====
  tar_target(
    df_mantel_tests,
    mantel_tests(
      dst_diet_all,
      dst_nutrient,
      dst_bioactive,
      dst_species,
      dst_urine,
      dst_faeces
    )
  ),

  
  # ==== PLOT MANUSCRIPT FIGURES ====
  tar_target(
    fig_two,
    paper_figure_two(
      plt_species_pcoa  = plt_species_pcoa,
      plt_species_dbrda = plt_spec_dbrda,
      plt_tax_shannon   = plt_species_shannon,
      plt_tax_bray      = plt_rare_scaled_species_bc_within,
      plt_fun_shannon   = plt_ko_shannon,
      plt_fun_bray      = plt_ko_bray,
      plt_mgs_cazymes   = plt_mgs_substrate$sig
    )
  ),
  tar_target(
    pth_fig_two,
    write_figure(
      plt = fig_two,
      pth = "figure_two",
      width = 8,
      height = 4,
      scale = 2
    )
  ),
  tar_target(
    fig_three,
    paper_figure_three(
      plt_volcano = plt_metab_volcano_revised_ident,
      plt_scfa_box = plt_scfa_box,
      plt_scfa_correlation = plt_scfa_correlation
    )
  ),
  tar_target(
    pth_fig_three,
    write_figure(
      plt = fig_three,
      pth = "figure_three",
      width = 6,
      height = 5,
      scale = 1.5
    )
  ),
  tar_target(
    fig_two_alternate,
    paper_figure_two_alternate(
      plt_species_pcoa  = plt_species_pcoa,
      plt_species_dbrda = plt_spec_dbrda,
      plt_tax_shannon   = plt_species_ad_all,
      plt_tax_bray      = plt_rare_scaled_species_bc_within,
      plt_fun_shannon   = plt_ko_ad_all,
      plt_fun_bray      = plt_ko_bray,
      plt_mgs_cazymes   = plt_mgs_substrate$sig
    )
  ),
  tar_target(
    pth_fig_two_alternate,
    write_figure(
      plt = fig_two_alternate,
      pth = "figure_two_alternate",
      width = 8,
      height = 4,
      scale = 2
    )
  ),
  tar_target(
    fig_three_alternate,
    paper_figure_three_alternate(
      plt_volcano = plt_metab_volcano,
      plt_scfa_box = plt_scfa_box,
      plt_scfa_correlation = plt_scfa_correlation
    )
  ),
  tar_target(
    pth_fig_three_alternate,
    write_figure(
      plt = fig_three_alternate,
      pth = "figure_three_alternate",
      width = 6,
      height = 2 * (7 / 3),
      scale = 1.5
    )
  ),

  # === GENERATE NETWORKS ===
  #' Data is split into two dataframe, one for samples after the high bioactive
  #' intervention, and another after the low bioactive intervention. Separate
  #' networks are then learnt for each of these.
  #'
  #' Code for network construction is provided, but making the figures does
  #' not depend on the outputs and it will not be run if you make any of the
  #' figures. This is as the network construction is quite time consuming.
  #' Instead, outputs saved as Rds are provided and used for the analysis.
  #' However, you can use `tar_make(lst_se_outputs)` to rerun the process,
  #' and `tar_make(pth_se_output)` to save these to disk. To use these new
  #' results in downstream analyses, you would need to replace the distributed
  #' results with these generated ones.

  # Split and tidy data
  tar_target(
    lst_metab_split,
    network_split_matrix(
      mat = mat_metab_faeces_sig,
      metadata = tbl_sample_metadata
    )
  ),
  tar_target(
    lst_pfam_split,
    network_split_matrix(
      mat = tbl_pfam |>
        column_to_rownames("L0") |>
        network_prevalence_filter(
          prevalence = 0.2,
          func_filt = TRUE,
          sample_md = tbl_sample_metadata
        ),
      metadata = tbl_sample_metadata
    )
  ),
  tar_target(
    lst_se_params,
    se_params()
  ),
  # Learn two networks
  tar_target(
    lst_net_high,
    paired_network_se(
      mat_a = lst_pfam_split$high,
      mat_b = lst_metab_split$high,
      se_params = lst_se_params
    )
  ),
  tar_target(
    lst_net_low,
    paired_network_se(
      mat_a = lst_pfam_split$low,
      mat_b = lst_metab_split$low,
      se_params = lst_se_params
    )
  ),
  # Output to files
  tar_target(
    vct_net_high_pths,
    network_write(
      lst_net_high,
      file.path("output/network_rerun"),
      condition_name = "after_high"
    )
  ),
  tar_target(
    vct_net_low_pths,
    network_write(
      lst_net_low,
      file.path("output/network_rerun"),
      condition_name = "after_low"
    )
  ),

  # ==== ANALYSIS OF NETWORKS ====
  # Association matrices
  tar_file_read(
    net_hb_ass, "data/derived/networks/after_high.association.Rds", readRDS(!!.x)
  ),
  tar_file_read(
    net_lb_ass, "data/derived/networks/after_low.association.Rds", readRDS(!!.x)
  ),
  # Igraph objects
  tar_file_read(
    net_lb_ig, "data/derived/networks/after_low.igraph.Rds", readRDS(!!.x)
  ),
  tar_file_read(
    net_hb_ig, "data/derived/networks/after_high.igraph.Rds", readRDS(!!.x)
  ),
  tar_file_read(
    net_tbl_pfam_md, "data/derived/networks/feature_md.EM.PFAML0.tsv",
    read_delim(!!.x)
  ),
  # Do global network property analysis
  tar_target(
    obj_net_an,
    netcomi_graph_analysis(
      hb_ass = net_hb_ass,
      lb_ass = net_lb_ass,
      fun_md = net_tbl_pfam_md,
      metab_md = tbl_metab_annotations,
      metab_ttest = tbl_metab_ttest,
      hb_ig = net_hb_ig,
      lb_ig = net_lb_ig
    )
  ),

  # List of metabolites for further investigation
  tar_target(
    tbl_taxa_corr,
    read_delim("data/derived/taxa/metabolite_taxa_weights.csv")
  ),
  tar_target(
    lst_interesting_metab,
    metabolite_review_list(
      tbl_metab_ttest,
      tbl_metab_peaks,
      tbl_metab_annotations,
      tbl_taxa_corr
    )
  ),

  # Check consistency of edges for HBcom-down metabolites
  tar_target(
    ig_hb,
    readRDS("data/derived/networks/after_high.igraph.Rds")
  ),
  tar_target(
    tbl_ig_hb_consistency,
    check_neighbourhood_consistency(
      ig_hb,
      vertex_names = c("fecal703neg", "fecal1214pos")
    )
  )
)
