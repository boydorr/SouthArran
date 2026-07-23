# Analysis

# workflow using adjusted survey data 
  # pivot wider 
  # remove data of species list from df_survey data that is not confirmed to a species level
  # combine any duplicated data 
  # add environmental grab data 
  # mutate column to add period of before and after (before = pre 2016) and (after = post 2016)

# Data preparation for analysis 
prepare_community_data <- function(df_survey,
                                   df_env,
                                   species_col = "accepted_name",
                                   abundance_col = "value",
                                   status_col = "worms_status",
                                   survey_site_col = "GrabSite",
                                   survey_year_col = "year",
                                   env_site_col = "GrabSite",
                                   env_year_col = "Year",
                                   cutoff_year = 2016,
                                   exclude_samples = NULL) {
  df_survey <- df_survey %>%
    mutate(sample_id = paste(.data[[survey_site_col]], .data[[survey_year_col]], sep = "_"))
  
  df_env <- df_env %>%
    mutate(sample_id = paste(.data[[env_site_col]], .data[[env_year_col]], sep = "_"))
  df_env <- df_env %>%
    mutate(sample_id = paste(.data[[env_site_col]], .data[[env_year_col]], sep = "_")) %>%
    distinct(sample_id, .keep_all = TRUE)
  
  # Exclude outlier samples 
  if (!is.null(exclude_samples)) {
    df_survey <- df_survey %>% filter(!sample_id %in% exclude_samples)
    df_env    <- df_env    %>% filter(!sample_id %in% exclude_samples)
  }
  
  # Cleaning the data
  df_clean <- df_survey %>%
    filter(.data[[status_col]] == "accepted",
           str_count(.data[[species_col]], "\\S+") >= 2,  # keep species-level records only
           .data[[abundance_col]] != "P") %>% 
    mutate(!!abundance_col := as.numeric(.data[[abundance_col]])) %>% 
    group_by(sample_id, .data[[species_col]]) %>%
    summarise(abundance = sum(.data[[abundance_col]], na.rm = TRUE),
              .groups = "drop")                        # combine duplicated records
  
  # pivot to appropriate form for analysis 
  community_wide <- df_clean %>%
    pivot_wider(names_from = all_of(species_col),
                values_from = abundance,
                values_fill = 0) %>%
    arrange(sample_id)
  
  # Joining grab data and species count data and adding period information 
  meta <- community_wide %>%
    select(sample_id) %>%
    left_join(df_env, by = "sample_id") %>%
    left_join(df_survey %>% distinct(sample_id, .data[[survey_year_col]]), by = "sample_id") %>%
    mutate(period = if_else(.data[[survey_year_col]] < cutoff_year, "Before", "After"),
           period = factor(period, levels = c("Before", "After")))
  stopifnot(nrow(meta) == nrow(community_wide))
  
  # Combine metadata/environmental data with the species abundance matrix
  # into a single wide table, matched on sample_id
  combined <- meta %>%
    left_join(community_wide, by = "sample_id")
  
  list(
    matrix     = community_wide %>% select(-sample_id) %>% as.data.frame(),
    meta       = meta,
    sample_ids = community_wide$sample_id,
    combined   = combined
  )
}


#NMDS ANALYSIS
run_nmds_analysis <- function(prep,
                              group_col = "period",
                              covariates = c("Mean_phi", "Depth"),
                              distance = "bray",
                              k = 2,
                              trymax = 100,
                              standardize = TRUE,
                              seed = 123,
                              run_permanova = TRUE,
                              run_betadisper = TRUE,
                              run_simper = FALSE,
                              plot_title = NULL) {
  
  comm_matrix <- prep$matrix
  meta <- prep$meta
  
  stopifnot(nrow(comm_matrix) == nrow(meta))
  if (!group_col %in% names(meta)) {
    stop(paste0("group_col '", group_col, "' not found in meta"))
  }
  
  # Track sample_id explicitly, not via rownames
  sample_ids <- meta$sample_id
  
  # Drop samples with missing covariate data (e.g. no sediment data collected)
  if (!is.null(covariates)) {
    complete_rows <- stats::complete.cases(meta[, covariates, drop = FALSE])
    if (any(!complete_rows)) {
      dropped_ids <- sample_ids[!complete_rows]
      warning(paste0(length(dropped_ids), " sample(s) with missing covariate data removed: ",
                     paste(dropped_ids, collapse = ", ")))
      comm_matrix <- comm_matrix[complete_rows, , drop = FALSE]
      meta <- meta[complete_rows, , drop = FALSE]
      sample_ids <- sample_ids[complete_rows]
    }
  }
  
  # Drop empty rows (samples with zero total abundance) to avoid metaMDS/vegdist errors
  row_sums <- rowSums(comm_matrix)
  if (any(row_sums == 0)) {
    empty_ids <- sample_ids[row_sums == 0]
    warning(paste0(length(empty_ids), " sample(s) with zero abundance removed: ",
                   paste(empty_ids, collapse = ", ")))
    comm_matrix <- comm_matrix[row_sums > 0, , drop = FALSE]
    meta <- meta[row_sums > 0, , drop = FALSE]
    sample_ids <- sample_ids[row_sums > 0]
  }
  
  # Drop species columns that are all-zero after filtering
  comm_matrix <- comm_matrix[, colSums(comm_matrix) > 0, drop = FALSE]
  
  # Standardization
  comm_std <- if (standardize) vegan::wisconsin(comm_matrix) else comm_matrix
  rownames(comm_std) <- sample_ids
  
  # NMDS
  set.seed(seed)
  nmds_result <- vegan::metaMDS(comm_std,
                                distance = distance,
                                k = k,
                                trymax = trymax,
                                autotransform = FALSE)
  
  # Extract site scores
  site_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
  site_scores$sample_id <- sample_ids   # always correct now, no rownames guessing
  
  plot_df <- site_scores %>%
    dplyr::left_join(meta, by = "sample_id")
  
  # Plot
  if (is.null(plot_title)) {
    plot_title <- paste0("NMDS of community structure (stress = ",
                         round(nmds_result$stress, 3), ")")
  }
  
  nmds_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = NMDS1, y = NMDS2,
                                                     color = .data[[group_col]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = "t", linewidth = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = plot_title, color = group_col)
  
  # PERMANOVA
  permanova_result <- NULL
  if (run_permanova) {
    explanatory_var <- if (!is.null(covariates)) {
      paste(c(covariates, group_col), collapse = " + ")
    } else {
      group_col
    }
    formula <- stats::as.formula(paste("comm_std ~", explanatory_var))
    permanova_result <- vegan::adonis2(formula, data = meta,
                                       method = distance, permutations = 999,
                                       by = "margin") 
  }
  
  # Betadisper (test of within-group dispersion homogeneity)
  betadisper_result <- NULL
  betadisper_anova <- NULL
  if (run_betadisper) {
    dist_matrix <- vegan::vegdist(comm_std, method = distance)
    betadisper_result <- vegan::betadisper(dist_matrix, meta[[group_col]])
    betadisper_anova <- stats::anova(betadisper_result)
  }
  
  # SIMPER
  if (run_simper) {
    simper_result  <- vegan::simper(comm_std, meta[[group_col]])
    simper_summary <- summary(simper_result, ordered = TRUE)
    
    comp_name <- names(simper_summary)[1]
    
    simper_top10 <- simper_summary[[comp_name]] %>%
      tibble::rownames_to_column("species") %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::rename(mean_before = av_before, mean_after = av_after) %>%  
      dplyr::mutate(
        comparison = comp_name,
        direction = dplyr::case_when(
          mean_after > mean_before ~ "Higher After",
          mean_after < mean_before ~ "Higher Before",
          TRUE ~ "No change"
        )
      )
  } else {
    simper_top10 <- NULL
  }
  
  # Return everything as a list
  list(
    nmds        = nmds_result,
    stress      = nmds_result$stress,
    scores      = plot_df,
    plot        = nmds_plot,
    permanova   = permanova_result,
    betadisper  = betadisper_result,
    betadisper_anova = betadisper_anova,
    simper      = simper_top10,
    comm_used   = comm_std,
    meta_used   = meta
  )
}

# Function for biological pattern measurements
analyses_bio_patterns <- function(comm, 
                                   meta,
                                   period_col = "period",
                                   covariates = c("Mean_phi", "Depth"),
                                   plot = TRUE) {
  
  # fail safes?
  stopifnot(nrow(comm) == nrow(meta))
  stopifnot(period_col %in% names(meta))
  stopifnot(all(covariates %in% names(meta)))
  
  # Diversity metrics
  richness <- specnumber(comm) # S
  shannon  <- diversity(comm, index = "shannon")  # H'
  pielou   <- shannon / log(richness) # J = H' / ln(S)
  
  div_df <- meta %>%
    mutate(richness = richness,
           shannon  = shannon,
           pielou   = pielou)
  
  div_df[[period_col]] <- factor(div_df[[period_col]],
                                 levels = c("Before", "After"))
  
  # Model formulas 
  explanatory_var <- paste(c(covariates, period_col), collapse = " + ") # management period and environmental covariates
  # Then using the metrics as response variables 
  f_richness <- as.formula(paste("richness ~", explanatory_var))
  f_shannon  <- as.formula(paste("shannon ~", explanatory_var))
  f_pielou   <- as.formula(paste("pielou ~", explanatory_var))
  
  # Fit models
  lm_richness <- lm(f_richness, data = div_df)
  lm_shannon  <- lm(f_shannon,  data = div_df)
  lm_pielou   <- lm(f_pielou,   data = div_df)
  
  models <- list(richness = lm_richness,
                 shannon  = lm_shannon,
                 pielou   = lm_pielou)
  
  # Put results in a table 
  results <- bind_rows(
    tidy(lm_richness) %>% mutate(response = "Richness"),
    tidy(lm_shannon)  %>% mutate(response = "Shannon-Wiener"),
    tidy(lm_pielou)   %>% mutate(response = "Pielou's evenness")
  )
  
  list(data = div_df, models = models, summary = results) 
}