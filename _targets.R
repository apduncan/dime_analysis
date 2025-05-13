# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse", "vegan", "glue", "patchwork", "ggplot2", "readxl",
    "ggrepel", "ggpubr", "pheatmap", "ggplotify", "patchwork",
    "usedist", "reticulate", "lsa", "rtk", "logger", "openxlsx",
    "colorspace"
  )
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  # Reading data
  # Metdata
  tar_target(
    pth_sample_metadata,
    "data/source/sample_metadata.csv",
    format = "file"
  ),
  tar_target(
    tbl_sample_metadata,
    read_metadata(pth_sample_metadata)
  ),
  # Metabolites
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
    pth_metab_ttest,
    "data/source/untargetted_metabolomics/allFeaturesPvalues2.xlsx"
  ),
  tar_target(
    tbl_metab_ttest,
    read_xlsx(pth_metab_ttest) |>
      mutate(full_id = paste0(Project, `Feature ID`, `Ion mode`))
  ),
  tar_target(
    pth_metab_peaks,
    paste0(c("data/source/untargetted_metabolomics/raw/DIME_results/save/",
             "peaklist_long_norm.rds"),
      collapse = ""
    )
  ),
  tar_target(
    tbl_metab_peaks,
    readRDS(pth_metab_peaks) |>
      mutate(full_id = paste0(gsub("DIME", "", project), feature_id, ionmode))
  ),
  tar_target(
    pth_metab_annotations,
    paste0(c("data/source/untargetted_metabolomics/",
             "selected_features_qib_annotated_manual_additions.xlsx"),
      collapse = ""
    )
  ),
  tar_target(
    tbl_metab_annotations,
    read_xlsx(pth_metab_annotations)
  ),
  tar_target(
    pth_metab_label_supplementary,
    metabolite_labelled_supplementary(
      tbl_metab_annotations = tbl_metab_annotations,
      output_path = "output/tables/supplementary/metabolite_labels.xlsx"
    )
  ),
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
  tar_target(
    pth_scfa,
    "data/source/scfa/SCFA_DIME.csv",
    format = "file"
  ),
  tar_target(
    tbl_scfa,
    read_delim(pth_scfa)
  ),
  # Diet
  tar_target(
    pth_food_diary,
    paste0(c("data/source/food_groups/",
             "cleaned_Nutritics_for_Fred_with_bioactives_170223.csv"),
      collapse = ""),
    format = "file"
  ),
  tar_target(
    pth_food_labelled,
    paste0(c("data/source/food_groups/",
             "Nutrtics_raw_Data_modified_wth_biosample_IDs_FB_10.08.23.xlsx"),
      collapse = ""),
    format = "file"
  ),
  tar_target(
    tbl_food_groups,
    read_food_groups(pth_food_diary, pth_food_labelled)
  ),
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
  tar_target(
    tbl_rare_scaled_species,
    scale_taxa(
      lst_species_alpha$matrix |>
        data.frame(check.names = FALSE) |>
        rownames_to_column("species"),
      tbl_cell_count
    )
  ),
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
  # Function
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
  # Analysis
  # Metabolite analysis
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

  # Food & nutrient
  tar_target(
    lst_food_freqs,
    food_group_frequency(tbl_food_groups, FOOD_GROUPS)
  ),

  # Distance matrices
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


  # Food & metabolite
  tar_target(
    tbl_food_metab_corr_hb,
    metabolite_food_correlation(
      tbl_metab_peaks,
      tbl_metab_ttest,
      tbl_metab_annotations,
      lst_food_freqs,
      "High Bioactive"
    )
  ),
  tar_target(
    plt_food_metabolite_corr_hb,
    plot_metabolite_food_correlation(
      tbl_food_metab_corr_hb,
      0.05,
      fontsize_row = 6,
      fontsize_col = 6,
      fontsize = 6,
      fontsize_number = 4
    )
  ),
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
  # Food, metabolite, taxa
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
  # Compose manuscript figures
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
      plt_volcano = plt_metab_volcano,
      plt_scfa_box = plt_scfa_box,
      plt_scfa_correlation = plt_scfa_correlation,
      plt_um_correlations = plt_food_metabolite_corr_hb
    )
  ),
  tar_target(
    pth_fig_three,
    write_figure(
      plt = fig_three,
      pth = "figure_three",
      width = 6,
      height = 7,
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

  # List of metabolites for review by Lars/Jan
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
  )
)
