# PCA for BTA

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(dplyr)

# Prepping df_surveydata for PCA
sample_lookup <- df_meta$meta %>%
  dplyr::distinct(sample_id, GrabSite, GrabSite_base, GrabNumber, year) # pulling sampleID - grabsite and year 

df_surveydata <- df_surveydata %>%
  dplyr::select(-dplyr::any_of(c("sample_id", "sample_id.x", "sample_id.y"))) %>%   # <-- ADD THIS LINE
  dplyr::left_join(sample_lookup, by = c("GrabSite", "GrabSite_base", "GrabNumber", "year"))

# checking these match
sum(is.na(df_surveydata$sample_id)) # if didn't have a match which ones
df_surveydata %>% dplyr::filter(is.na(sample_id)) %>% dplyr::distinct(GrabSite, GrabSite_base, GrabNumber, year)

# Creating weighted values function
standardise_traits <- function(filled_matrix, trait_groups) {
  
  out <- filled_matrix
  
  for (grp in names(trait_groups)) {
    cols    <- trait_groups[[grp]]
    grp_sum <- rowSums(out[cols], na.rm = TRUE)
    
    scaled <- out[cols] / grp_sum
    scaled[grp_sum == 0, ] <- 0   # keep "no information" taxa as 0 instead of getting an NA value if no proportion
    
    out[cols] <- scaled
  }
  
  out
}

# Building weighted matrix 
  # converting raw abundance to proportion of that sample's total
build_sample_trait_matrix <- function(abundance_data,
                                      trait_matrix_std,
                                      taxon_col    = "bta_name",
                                      sample_col   = "GrabSite",
                                      abund_col    = "value",
                                      trait_cols   = NULL,
                                      trait_groups = NULL) {
  
  if (is.null(trait_cols)) {
    if (is.null(trait_groups)) {
      stop("Supply either trait_cols or trait_groups so build_sample_trait_matrix() ",
           "knows which columns are real traits (not taxonomy/metadata columns).")
    }
    trait_cols <- unlist(trait_groups, use.names = FALSE)
  }
  
  rel_abund <- abundance_data %>%
    select(all_of(c(taxon_col, sample_col, abund_col))) %>% # standardising column names
    rename(taxon = all_of(taxon_col),
           sample = all_of(sample_col),
           abund  = all_of(abund_col)) %>%
    mutate(abund = suppressWarnings(as.numeric(abund))) %>%
    filter(!is.na(abund), abund > 0) %>% # dropping 0 abundance
    group_by(sample) %>%
    mutate(rel_abund = abund / sum(abund)) %>%
    ungroup()
  
  trait_lookup <- trait_matrix_std %>% # joining trait profile onto abundance counts
    rename(taxon = all_of(taxon_col)) %>%
    select(taxon, all_of(trait_cols))
  
  joined <- rel_abund %>%
    left_join(trait_lookup, by = "taxon")
  
  missing_taxa <- joined %>%
    filter(if_any(all_of(trait_cols), is.na)) %>%
    distinct(taxon)
  
  if (nrow(missing_taxa) > 0) {
    warning(nrow(missing_taxa),
            " taxa present in abundance data have no match in the trait matrix: ",
            paste(missing_taxa$taxon, collapse = ", "),
            "\n  -> these will contribute NA and be dropped from weighting.") # These are the species that didn't make the inclusion criteria for building BTA matrix
  }
  
  sample_trait <- joined %>%
    filter(!if_any(all_of(trait_cols), is.na)) %>%
    mutate(across(all_of(trait_cols), ~ .x * rel_abund)) %>% # Getting community weighted mean 
    group_by(sample) %>%
    summarise(across(all_of(trait_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename(!!sample_col := sample)
  
  sample_trait
}

#  Run the PCA
run_bta_pca <- function(sample_trait_matrix, sample_col = "GrabSite") {
  mat <- sample_trait_matrix %>%
    column_to_rownames(sample_col)
  
  zero_var <- names(mat)[sapply(mat, function(x) sd(x, na.rm = TRUE) == 0)]
  if (length(zero_var) > 0) {
    message("Removing zero-variance trait column(s) before PCA: ",
            paste(zero_var, collapse = ", "))
    mat <- mat[, setdiff(names(mat), zero_var), drop = FALSE]
  }
  
  pca_res <- FactoMineR::PCA(mat, scale.unit = TRUE, ncp = 5, graph = FALSE)
  
  scores <- as_tibble(pca_res$ind$coord, rownames = sample_col) %>%
    rename_with(~ str_replace(.x, "^Dim\\.", "PC"), starts_with("Dim"))
  
  var_contrib <- as_tibble(pca_res$var$contrib, rownames = "trait") %>%
    rename_with(~ str_replace(.x, "^Dim\\.", "PC"), starts_with("Dim"))
  var_coord <- as_tibble(pca_res$var$coord, rownames = "trait") %>%
    rename_with(~ str_replace(.x, "^Dim\\.", "PC"), starts_with("Dim"))
  
  list(
    pca         = pca_res,
    scores      = scores,
    var_contrib = var_contrib,
    var_coord   = var_coord
  )
  

}


# Fit a linear model --> PC1/PC2 ~ management period + phi + depth

run_pca_lm <- function(pca_scores, metadata,
                       sample_col = "GrabSite",
                       period_col = "management_period",
                       phi_col    = "Mean_phi",
                       depth_col  = "Depth") {
  
  df <- pca_scores %>%
    left_join(
      metadata %>% select(all_of(c(sample_col, period_col, phi_col, depth_col))),
      by = sample_col
    ) %>%
    rename(period = all_of(period_col),
           phi    = all_of(phi_col),
           depth  = all_of(depth_col)) %>%
    mutate(period = factor(period))  
  
  lm_PC1 <- lm(PC1 ~ period + phi + depth, data = df)
  lm_PC2 <- lm(PC2 ~ period + phi + depth, data = df)
  
  message("--- PC1 ~ period + phi + depth ---")
  print(summary(lm_PC1))
  message("\n--- PC2 ~ period + phi + depth ---")
  print(summary(lm_PC2))
  
  list(
    data        = df,
    lm_PC1      = lm_PC1,
    lm_PC2      = lm_PC2,
    summary_PC1 = summary(lm_PC1),
    summary_PC2 = summary(lm_PC2)
  )
}


# Putting it all together and visualising
run_bta_pca_workflow <- function(filled_matrix,
                                  trait_groups,
                                  abundance_data,
                                  metadata_data  = NULL,
                                  taxon_col      = "bta_name",
                                  sample_col     = "GrabSite",
                                  abund_col      = "value",
                                  station_col    = NULL,
                                  station_filter = NULL,
                                  period_col     = "management_period",
                                  phi_col        = "Mean_phi",
                                  depth_col      = "Depth",
                                  top_n_contrib  = 15,
                                  run_diagnostics = TRUE) {
  
    # restrict to GrabSites where station starts with a given letter looking for mearle = D or sediment = T
    if (!is.null(station_filter)) {
      if (is.null(station_col)) {
        stop("station_col must be supplied (e.g. 'GrabSite_station') when using station_filter.")
      }
      
      message("Filtering to stations starting with: '", station_filter, "'")
      
      abundance_data <- abundance_data %>%
        dplyr::filter(stringr::str_starts(.data[[station_col]], station_filter))
      
      if (!is.null(metadata_data)) {
        metadata_data <- metadata_data %>%
          dplyr::filter(stringr::str_starts(.data[[station_col]], station_filter))
      }
      
      n_sites <- abundance_data %>% dplyr::distinct(.data[[sample_col]]) %>% nrow()
      message("  -> ", n_sites, " GrabSite(s) retained.")
      if (n_sites == 0) {
        stop("No GrabSites matched station_filter = '", station_filter,
             "' in column '", station_col, "'. Check spelling/column name.")
      }
    }

  message("Step 1: Standardising trait matrix...")
  trait_matrix_std <- standardise_traits(filled_matrix, trait_groups)
  
  message("Step 2: Building sample x trait (community-weighted) matrix...")
  sample_trait_matrix <- build_sample_trait_matrix(
    abundance_data   = abundance_data,
    trait_matrix_std = trait_matrix_std,
    taxon_col        = taxon_col,
    sample_col       = sample_col,
    abund_col        = abund_col,
    trait_groups     = trait_groups
  )
  
  message("Step 3: Running PCA...")
  pca_out <- run_bta_pca(sample_trait_matrix, sample_col = sample_col)
  
  message("Step 4: Building metadata for plots and models...")
  if (is.null(metadata_data)) {
    metadata_data <- abundance_data %>%
      distinct(across(all_of(c(sample_col, period_col, phi_col, depth_col))))
  }
  
  habillage_vec <- metadata_data[[period_col]][
    match(rownames(pca_out$pca$ind$coord), metadata_data[[sample_col]])
  ]
  
  message("Step 5: Generating plots...")
  plot_scree <- fviz_eig(pca_out$pca, addlabels = TRUE)
  
  # Individuals plot (samples only)
  plot_ind <- fviz_pca_ind(
    pca_out$pca,
    habillage   = habillage_vec,
    addEllipses = TRUE,
    ellipse.type = "confidence",
    label       = "none", 
    pointsize   = 3,
    repel       = TRUE
  ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = NULL)
  
  # Getting all traits 
  all_traits <- pca_out$var_contrib$trait 
  
  # Top traits (by cos2 = quality of representation on the PC1-PC2 plane)
  cos2_df <- as.data.frame(pca_out$pca$var$cos2[, 1:2]) %>%
    tibble::rownames_to_column("trait") %>%
    mutate(cos2_sum = Dim.1 + Dim.2) %>%
    arrange(desc(cos2_sum))
  
  top_traits <- cos2_df %>%
    slice_head(n = top_n_contrib) %>%
    pull(trait)
  
  # Full plot for variables (all traits, arrows removed)
  plot_var_arrows <- fviz_pca_var(
    pca_out$pca,
    select.var = list(name = all_traits),
    col.var    = "black",
    repel      = TRUE,
    arrowsize  = 0.5,
    labelsize  = 3
  ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position   = "none"
    ) +
    ggplot2::labs(
      title = NULL,
      x = paste0("PC1 (", round(pca_out$pca$eig[1, 2], 1), "%)"),
      y = paste0("PC2 (", round(pca_out$pca$eig[2, 2], 1), "%)")
    )
  
  # Heatmap for traits 
  var_cor <- pca_out$pca$var$cor[all_traits, 1:2]
  
  var_cor_df <- as.data.frame(var_cor) %>%
    tibble::rownames_to_column("trait") %>%
    tidyr::pivot_longer(cols = c(Dim.1, Dim.2), names_to = "PC", values_to = "corr") %>%
    mutate(
      PC    = recode(PC, Dim.1 = "PC1", Dim.2 = "PC2"),
      trait = factor(trait, levels = rev(all_traits))
    )
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    
    # narrower tiles + narrower panel width in the layout below
    plot_corr_heat <- ggplot2::ggplot(var_cor_df, ggplot2::aes(x = PC, y = trait)) +
      ggplot2::geom_tile(ggplot2::aes(fill = corr), color = "white",
                         linewidth = 0.2, width = 0.95, height = 0.95) +
      ggplot2::scale_fill_gradient2(
        low = "#B2182B", mid = "white", high = "#2166AC",
        midpoint = 0, limits = c(-1, 1), name = "Corr"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title  = ggplot2::element_blank(),
        panel.grid  = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
        legend.key.width = grid::unit(0.4, "cm")
      )
    
    plot_var <- plot_var_arrows + plot_corr_heat +
      patchwork::plot_layout(widths = c(2.2, 0.4))  # heatmap panel made narrower
    
  } else {
    message("Install 'patchwork' for the Figure 8-style combined panel.")
    plot_var <- plot_var_arrows
  }
  
  # Curated top-N plot, now colour-coded by cos2 (representation quality)
  plot_var_top <- fviz_pca_var(
    pca_out$pca,
    select.var = list(name = top_traits),
    col.var    = "cos2",
    gradient.cols = c("#B2182B", "white", "#2166AC"),
    repel      = TRUE,
    arrowsize  = 0.6,
    labelsize  = 3.5
  ) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Top contributing / best-represented traits", color = "cos2")
  
  plot_contrib_PC1 <- fviz_contrib(pca_out$pca, choice = "var", axes = 1, top = top_n_contrib)
  plot_contrib_PC2 <- fviz_contrib(pca_out$pca, choice = "var", axes = 2, top = top_n_contrib)
  
  print(plot_ind)
  print(plot_var)       # full-trait biplot + heatmap (Figure 8 style)
  print(plot_var_top)  
  print(plot_contrib_PC1)
  print(plot_contrib_PC2)
  
  message("Step 6: Fitting PC1/PC2 ~ period + phi + depth linear models...")
  lm_out <- run_pca_lm(
    pca_scores = pca_out$scores,
    metadata   = metadata_data,
    sample_col = sample_col,
    period_col = period_col,
    phi_col    = phi_col,
    depth_col  = depth_col
  )
  
  if (run_diagnostics) {
    message("Step 7: Model diagnostics...")
    par(mfrow = c(2, 2))
    plot(lm_out$lm_PC1); mtext("PC1 model diagnostics", side = 3, line = -1, outer = TRUE)
    plot(lm_out$lm_PC2); mtext("PC2 model diagnostics", side = 3, line = -1, outer = TRUE)
    par(mfrow = c(1, 1))
    
    if (requireNamespace("car", quietly = TRUE)) {
      message("VIF - PC1 model:")
      print(car::vif(lm_out$lm_PC1))
      message("VIF - PC2 model:")
      print(car::vif(lm_out$lm_PC2))
    } else {
      message("Install the 'car' package to check VIF for phi/depth collinearity.")
    }
  }
  invisible(list(
    sample_trait_matrix = sample_trait_matrix,
    pca_out  = pca_out,
    plots = list(
      scree       = plot_scree,
      ind         = plot_ind,
      var         = plot_var,
      var_top     = plot_var_top,
      contrib_PC1 = plot_contrib_PC1,
      contrib_PC2 = plot_contrib_PC2
    ),
    lm_out   = lm_out
  ))
}
