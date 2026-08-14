# Functional Redundancy 
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(broom)
library(vegan)
library(cluster)
library(ape)
library(picante)
library(car)
library(lmerTest)


# 1. Scatter + fitted line + R^2 annotation + slope 
plot_fd_h_scatter <- function(fd_fr, group_col) {
  #' Scatter plot visualisation 
  #' 
  #' @description Function that will plot FD/H' to represent functional redundancy 
  #' @param fd_fr Dataframe with FD values
  #' @param group_col Name string in fd_fr
  #' 
  #' @return A plot of the FD/H' slope
  
  fits <- fd_fr %>%
    dplyr::group_by(.data[[group_col]]) %>% # grouping by "group_col" (period or protection)
    dplyr::summarise(
      beta = coef(lm(FD ~ H))[["H"]],
      r2   = summary(lm(FD ~ H))$r.squared,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = paste0("\u03b2 = ", round(beta, 3),# building a text string for Beta
                     ", R\u00b2 = ", round(r2, 3))
    )
  
  ggplot2::ggplot(fd_fr, ggplot2::aes(x = H, y = FD)) + # plotting the relationship FD/H'
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black") +
    ggplot2::geom_text(data = fits, ggplot2::aes(x = -Inf, y = Inf, label = label),
                       hjust = -0.1, vjust = 1.5, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ .data[[group_col]]) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "H' (Shannon diversity)", y = "FD (Functional diversity)")
}

#2.  Interaction test
test_slope_difference <- function(fd_fr, 
                                  group_col,
                                  station_col = NULL,
                                  ref_level = NULL) {
  #' @description Function that tests whether the FD~H' slope differ significantly between groups
  #' 
  #' @return A statistical output 
  
  df <- fd_fr
  df[[group_col]] <- factor(df[[group_col]])
  
  if (!is.null(ref_level)) {
    df[[group_col]] <- stats::relevel(df[[group_col]], ref = ref_level)
  }
  
  use_mixed <- !is.null(station_col) # deciding if to fit a mixed effect model
  
  if (use_mixed) {
    df[[station_col]] <- factor(df[[station_col]])
    f <- as.formula(
      paste0("FD ~ H * ", group_col, " + (1 | ", station_col, ")")
    )
    mod <- lmerTest::lmer(f, data = df, REML = TRUE)
    result <- list(
      model   = mod,
      anova   = anova(mod),
      summary = summary(mod)
    )
  } else {
    f <- as.formula(paste("FD ~ H *", group_col))
    mod <- lm(f, data = df)
    result <- list(
      model   = mod,
      anova   = car::Anova(mod, type = 3),
      summary = summary(mod)
    )
  }
  result
}

# 3. Models
run_fr_lm <- function(fd_fr,
                      period_col = "period",
                      phi_col = "Mean_phi",
                      depth_col = "Depth",
                      station_col = NULL,
                      ref_level = NULL) {
  #' @description A functional that performs linear model with FR itself (response) against period/protection + phi/depth (explanatory variable)
  #' 
  #' @return Linear model
  
  df <- fd_fr
  df[[period_col]] <- factor(df[[period_col]])
  
  # setting reference level
  if (!is.null(ref_level)) {                    
    df[[period_col]] <- stats::relevel(df[[period_col]], ref = ref_level)
  }
  
  if (is.null(station_col)) {
    
    f <- reformulate(
      c(period_col, phi_col, depth_col),
      response = "FR") 
    
    mod <- lm(f, data = df)
    
  } else {
    
    df[[station_col]] <- factor(df[[station_col]])
    
    f <- as.formula(
      paste(
        "FR ~",
        period_col, "+",
        phi_col, "+",
        depth_col,
        "+ (1|", station_col, ")"))
    
    mod <- lmerTest::lmer(f, data = df)
  }
  
  summary(mod)
}


#4. Wrapping the functions together 
 run_fd_fr_workflow <- function(abundance_data, trait_matrix_std, trait_groups,
                               metadata_data, 
                               group_col = "period",
                               taxon_col_abund = "accepted_name",
                               taxon_col_trait = "bta_name",
                               sample_col = "GrabSite",
                               station_col    = NULL,
                               station_filter = NULL,
                               phi_col   = "Mean_phi",
                               depth_col = "Depth",
                               ref_level = NULL,
                               fr_random_station = FALSE) {
  #' @description 
  #' @param 
  #' @param 
  #' 
  #' @return A list of outputs. Scatter plot, statistics, and linear models
  
  # Optional station filtering
  if (!is.null(station_filter)) {
    if (is.null(station_col)) {
      stop("station_col must be supplied (e.g. 'GrabSite_station') when using station_filter.")
    }
    message("Filtering to stations starting with: '", station_filter, "'")
    abundance_data <- abundance_data %>%
      dplyr::filter(stringr::str_starts(.data[[station_col]], station_filter))
    metadata_data <- metadata_data %>%
      dplyr::filter(stringr::str_starts(.data[[station_col]], station_filter))
    n_sites <- abundance_data %>% dplyr::distinct(.data[[sample_col]]) %>% nrow()
    message("  -> ", n_sites, " GrabSite(s) retained.")
    if (n_sites == 0) {
      stop("No GrabSites matched station_filter = '", station_filter,
           "' in column '", station_col, "'. Check spelling/column name.")
    }
  }
  
   # Building input for FD calculation
  trait_cols <- unlist(trait_groups, use.names = FALSE)
  sp_traits <- build_species_trait_matrix(trait_matrix_std, taxon_col_trait, trait_cols)
  tree      <- build_functional_dendrogram(sp_traits)
  comm      <- build_sample_species_matrix(abundance_data, taxon_col_abund, sample_col)
  
  # computing FD/FR
  fd_fr <- calc_fd_fr(comm, tree) %>%
    dplyr::left_join(metadata_data, by = c("sample" = sample_col))
  
  # Setting a reference level 
  if (!is.null(ref_level)) {
    fd_fr[[group_col]] <- stats::relevel(factor(fd_fr[[group_col]]), ref = ref_level)
  } else {
    fd_fr[[group_col]] <- factor(fd_fr[[group_col]])
  }
  
  # FD ~ H' scatter with R^2, faceted by group
  p_fd_h_scatter <- plot_fd_h_scatter(fd_fr, group_col)
  
  # Regression tables
  lm_by_group <- fd_fr %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_modify(~ broom::tidy(lm(FD ~ H, data = .x))) %>%
    dplyr::ungroup()
  
  lm_by_group_fit <- fd_fr %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_modify(~ broom::glance(lm(FD ~ H, data = .x))) %>%
    dplyr::ungroup()
  
  lm_by_group_std <- fd_fr %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(FD_z = as.numeric(scale(FD)), H_z = as.numeric(scale(H))) %>%
    dplyr::group_modify(~ broom::tidy(lm(FD_z ~ H_z, data = .x))) %>%
    dplyr::ungroup()
  
  # Interaction test - use the top-level function, pass ref_level through
  slope_diff_test <- test_slope_difference(fd_fr, group_col,
                                           station_col = if (fr_random_station) station_col else NULL,
                                           ref_level = ref_level)
  
  # FR model - pass ref_level through
  fr_model <- run_fr_lm(fd_fr,
                        period_col = group_col,
                        phi_col    = phi_col,
                        depth_col  = depth_col,
                        station_col = if (fr_random_station) station_col else NULL,
                        ref_level  = ref_level)
  
  list(
    data              = fd_fr,
    plots             = list(FD_H_scatter = p_fd_h_scatter),
    lm_by_group       = lm_by_group,
    lm_by_group_fit   = lm_by_group_fit,
    lm_by_group_std   = lm_by_group_std,
    slope_diff_test   = slope_diff_test,
    fr_model          = fr_model
  )
}