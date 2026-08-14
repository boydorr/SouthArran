# Libraries loading 
library(dplyr) 
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(lmerTest)
library(car)
library(broom)
library(purrr)
library(tibble)
library(stringr)

# Workflow
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
                                   station_col = "GrabSite_station",
                                   cutoff_year = 2016,
                                   exclude_samples = NULL) {
  #' Data preparation pipeline 
  #' @description A functional that takes raw survey data (species and grabs) ad creates a clean abundance matrix for future analysis
  #' @param df_survey species abundance survey data
  #' @param df_env environmental data (grab location, and sediment analysis) 
  #' 
  #' @return Four objects in a list (matrix, meta, sample_ids, and combined which is metadata + species matrix combined)
  
  df_survey <- df_survey %>%
    mutate(sample_id = paste(.data[[survey_site_col]], .data[[survey_year_col]], sep = "_"))
  
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
  
  # Pivot to appropriate form for analysis (wide)
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
           period = factor(period, levels = c("Before", "After")),
           GrabSite_station = factor(GrabSite_station)) 
  
  stopifnot(nrow(meta) == nrow(community_wide))
  
  # Combine metadata/environmental data with the species abundance matrix into a single wide table
  combined <- meta %>%
    left_join(community_wide, by = "sample_id") #  matched on sample_id
  
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
                              ref_level = NULL,
                              run_permanova = TRUE,
                              run_betadisper = TRUE,
                              run_simper = FALSE,
                              station_col    = NULL,
                              station_filter = NULL,
                              exclude_samples = NULL,
                              plot_title = NULL,
                              colors = NULL) {
  #' Runs a multivariate community ecology workflow
  #' 
  #' @description  This function creates NMDS ordination and additional standard supporting tests such as PERMANOVA, betadisper, SIMPER
  #' @param prep The cleaned data prepped from the previous function
  #' 
  #' @return A list of raw meta object, NMDS plot, PERMANOVA, betadisper anova, and SIMPER tables
  
  comm_matrix <- prep$matrix
  meta <- prep$meta
  
  if (!is.null(ref_level)) {
    meta[[group_col]] <- relevel(factor(meta[[group_col]]), ref = ref_level)
  }
  
  stopifnot(nrow(comm_matrix) == nrow(meta))
  if (!group_col %in% names(meta)) {
    stop(paste0("group_col '", group_col, "' not found in meta"))
  }
  
  if (!is.null(station_filter)) {
    if (is.null(station_col) || !station_col %in% names(meta)) {
      stop("station_col must be a valid column in meta when station_filter is set")
    }
    keep_station <- !is.na(meta[[station_col]]) & 
      stringr::str_starts(meta[[station_col]], station_filter)
    comm_matrix <- comm_matrix[keep_station, , drop = FALSE]
    meta        <- meta[keep_station, , drop = FALSE]
  }
  
  # Track sample_id, not via rownames
  sample_ids <- meta$sample_id
  
  # Exclude outlier samples 
  if (!is.null(exclude_samples)) {
    
    keep <- !(sample_ids %in% exclude_samples)
    
    comm_matrix <- comm_matrix[keep, , drop = FALSE]
    meta        <- meta[keep, , drop = FALSE]
    sample_ids  <- sample_ids[keep]
  }
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
  
  # Drop empty rows (samples with zero total abundance) to avoid errors
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
  
  rownames(meta) <- meta$sample_id
  
  
  # NMDS
  set.seed(seed)
  nmds_result <- vegan::metaMDS(comm_std,
                                distance = distance,
                                k = k,
                                trymax = trymax,
                                autotransform = FALSE)
  
  # Extract site scores
  site_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
  site_scores$sample_id <- sample_ids   # no rownames guessing
  
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
  
  if (!is.null(colors)) {
    nmds_plot <- nmds_plot + ggplot2::scale_color_manual(values = colors)
  }
  
  # PERMANOVA
  if (run_permanova) {
    explanatory_var <- if (!is.null(covariates)) {
      paste(c(covariates, group_col), collapse = " + ")
    } else {
      group_col
    }
    formula <- stats::as.formula(paste("comm_std ~", explanatory_var))
    
    use_strata <- "GrabSite_station" %in% names(meta) && !anyNA(meta$GrabSite_station)
    
    meta <- meta[rownames(comm_std), , drop = FALSE]
    
    permanova_result <- vegan::adonis2(
      formula,
      data = meta,
      method = distance,
      permutations = 999,
      by = "margin",
      strata = if (use_strata) droplevels(factor(meta$GrabSite_station)) else NULL)
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
    
    # Pull every comparison into one long table
    simper_all <- purrr::imap_dfr(simper_summary, function(comp_df, comp_name) {
      comp_df %>%
        tibble::rownames_to_column("species") %>%
        dplyr::mutate(comparison = comp_name)
    })
    
    # Per-comparison top 10
    simper_top10 <- simper_all %>%
      dplyr::group_by(comparison) %>%
      dplyr::group_modify(~ {
        grp_labels <- strsplit(.y$comparison[1], "_")[[1]]
        avg_cols   <- grep("^av", names(.x), value = TRUE)
        .x %>%
          dplyr::slice_head(n = 10) %>%
          dplyr::mutate(
            direction = dplyr::case_when(
              .data[[avg_cols[2]]] > .data[[avg_cols[1]]] ~ paste("Higher in", grp_labels[2]),
              .data[[avg_cols[2]]] < .data[[avg_cols[1]]] ~ paste("Higher in", grp_labels[1]),
              TRUE ~ "No change"
            )
          )
      }) %>%
      dplyr::ungroup()
    
    # Overall top 10 unique species, ranked on significance
    simper_top10_overall <- simper_top10 %>%
      dplyr::filter(p < 0.05) %>%
      dplyr::group_by(species) %>%
      dplyr::slice_max(average, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(average)) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::select(species, average, sd, ratio, comparison, direction, p)
    
  } else {
    simper_all <- NULL
    simper_top10 <- NULL
    simper_top10_overall <- NULL
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
    simper_all  = simper_all,
    simper_top10_overall = simper_top10_overall,
    comm_used   = comm_std,
    meta_used   = meta)
}



# Function for biological pattern measurements
analyses_bio_patterns <- function(comm, 
                                  meta,
                                  period_col = "period",
                                  ref_level = NULL,
                                  covariates = c("Mean_phi", "Depth"),
                                  station_col = "GrabSite_station",
                                  station_filter = NULL,
                                  plot = TRUE) {
  #'@description A function that computes biodiversity measurements (Richness, Shannon Diversity, Pielou's evenness) from a community matrix and test them against group covariates (period/ protection levels) using linear models
  #'@param comm species abundance matrix (species x samples)
  #'@param meta metadata frame
  #'
  #'@return A list with data, models, summary 
  
  # fail safes
  stopifnot(nrow(comm) == nrow(meta))
  stopifnot(period_col %in% names(meta))
  stopifnot(all(covariates %in% names(meta)))
  if (!is.null(station_col)) stopifnot(station_col %in% names(meta))
  
  # Diversity metrics
  richness <- specnumber(comm) # S
  H        <- diversity(comm, index = "shannon")
  shannon  <- exp(H) # exponential Shannon (effective species)
  pielou   <- ifelse(richness > 1, H / log(richness), NA_real_) 
  
  div_df <- meta %>%
    mutate(richness = richness,
           shannon  = shannon,
           pielou   = pielou)
  
  div_df[[period_col]] <- as.factor(div_df[[period_col]])
  
  # Setting reference level 
  if (!is.null(ref_level)) {
    div_df[[period_col]] <- relevel(div_df[[period_col]], ref = ref_level)
  }
  
  explanatory_var <- paste(c(covariates, period_col), collapse = " + ")
  
  use_mixed <- !is.null(station_col)
  
  if (use_mixed) {
    div_df[[station_col]] <- as.factor(div_df[[station_col]])
    rhs <- paste0(explanatory_var, " + (1 | ", station_col, ")")
  } else {
    rhs <- explanatory_var
  }
  f_richness <- as.formula(paste0("richness ~ ", rhs))
  f_shannon  <- as.formula(paste0("shannon ~ ",  rhs))
  f_pielou   <- as.formula(paste0("pielou ~ ",   rhs))
  
  if (use_mixed) {
    lm_richness <- lmerTest::lmer(f_richness, data = div_df, REML = TRUE)
    lm_shannon  <- lmerTest::lmer(f_shannon,  data = div_df, REML = TRUE)
    lm_pielou   <- lmerTest::lmer(f_pielou,   data = div_df, REML = TRUE)
  } else {
    lm_richness <- lm(f_richness, data = div_df)
    lm_shannon  <- lm(f_shannon,  data = div_df)
    lm_pielou   <- lm(f_pielou,   data = div_df)
  }
  
  models <- list(richness = lm_richness,
                 shannon  = lm_shannon,
                 pielou   = lm_pielou)
  
  results <- bind_rows(
    tidy(lm_richness) %>% mutate(response = "Richness"),
    tidy(lm_shannon) %>% mutate(response = "Exponential Shannon (Hill N1)"),
    tidy(lm_pielou)   %>% mutate(response = "Pielou's evenness")
  )
  
  list(data = div_df, models = models, summary = results)
}