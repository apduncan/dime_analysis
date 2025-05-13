mantel_tests <- function(
  dst_diet,
  dst_nutrient,
  dst_bioactive,
  dst_taxa,
  dst_urine,
  dst_faeces
) {
  #' Mantel tests for relationships between diet and metabolome/microbiome
  #'
  #' @returns Tibbles with columns for p-value, correlation, and the
  #' types of data

  mantel_to_vec <- function(mantel_res, a, b) {
    return(setNames(
      c(mantel_res$signif, mantel_res$statistic, a, b),
      c("pval", "corr", "other_type", "diet_type")
    ))
  }
  mantel_to_diet <- function(dist_mat, label) {
    # Subset distances to only those shared by dietary and other distances
    shared_samples <- intersect(
      attr(dist_mat, "Labels"),
      attr(dst_diet, "Labels")
    )
    dst_diet_sub <- dist_subset(dst_diet, shared_samples)
    dst_nutrient_sub <- dist_subset(dst_nutrient, shared_samples)
    dst_bioactive_sub <- dist_subset(dst_bioactive, shared_samples)
    dst_mat_sub <- dist_subset(dist_mat, shared_samples)
    res_all <- vegan::mantel(dst_mat_sub, dst_diet)
    res_nutrient <- vegan::mantel(dst_mat_sub, dst_nutrient)
    res_bioactive <- vegan::mantel(dst_mat_sub, dst_bioactive)
    do.call(rbind,
      list(
        mantel_to_vec(res_all, label, "all"),
        mantel_to_vec(res_nutrient, label, "nutrient"),
        mantel_to_vec(res_bioactive, label, "bioactive")
      )
    )
  }
  # Diet to microbiome composition
  species_dist_subset <- dist_subset(dst_taxa, attr(dst_diet, "Labels"))
  all_mantel_results <- do.call(
    rbind,
    list(
      mantel_to_diet(species_dist_subset, "taxa"),
      mantel_to_diet(dst_faeces, "faecal"),
      mantel_to_diet(dst_urine, "urine")
    )
  ) |> as.data.frame()
  all_mantel_results$pval <- as.numeric(all_mantel_results$pval)
  all_mantel_results$corr <- as.numeric(all_mantel_results$corr)
  # Adjust p values
  all_mantel_results$qval <- p.adjust(all_mantel_results$pval, method = "BH")
  all_mantel_results
}