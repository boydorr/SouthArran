# Species accumulation curves 
library(iNEXT)
library(dplyr)
library(tidyverse)

# Cleaning the data- using df_meta
# Data preparation specifically for presence-absence analysis (e.g. SAC)
prepare_PA_data <- function(df_survey,
                            df_env,
                            species_col = "accepted_name",
                            abundance_col = "value",
                            status_col = "worms_status",
                            survey_site_col = "GrabSite",
                            survey_year_col = "year",
                            env_site_col = "GrabSite",
                            env_year_col = "Year",
                            exclude_samples = NULL) {
  
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
  
  # Cleaning the data - "P" (present, not counted) is treated as a presence, not dropped
  df_clean <- df_survey %>%
    filter(.data[[status_col]] == "accepted",
           str_count(.data[[species_col]], "\\S+") > 2) %>%  # keep species-level records only
    mutate(!!abundance_col := case_when(
      .data[[abundance_col]] == "P" ~ 1,
      TRUE ~ suppressWarnings(as.numeric(.data[[abundance_col]]))
    )) %>%
    group_by(sample_id, .data[[species_col]]) %>%
    summarise(abundance = sum(.data[[abundance_col]], na.rm = TRUE),
              .groups = "drop")                        # combine duplicated records
  
  # pivot to wide form
  community_wide <- df_clean %>%
    pivot_wider(names_from = all_of(species_col),
                values_from = abundance,
                values_fill = 0) %>%
    arrange(sample_id) %>%
    mutate(across(-sample_id, ~ as.integer(. > 0)))     # force strict 0/1
  
  # Joining grab data and adding period information (creating a meta date for presence)
  meta <- community_wide %>%
    select(sample_id) %>%
    left_join(df_env, by = "sample_id") %>%
    left_join(df_survey %>% distinct(sample_id, .data[[survey_year_col]]), by = "sample_id") %>%
    mutate(period = if_else(.data[[survey_year_col]] < cutoff_year, "Before", "After"),
           period = factor(period, levels = c("Before", "After")))
  stopifnot(nrow(meta) == nrow(community_wide))
  
  # Combine metadata/environmental data with the presence-absence matrix
  combined <- meta %>%
    left_join(community_wide, by = "sample_id")
  
  list(
    matrix     = community_wide %>% select(-sample_id) %>% as.data.frame(),
    meta       = meta,
    sample_ids = community_wide$sample_id,
    combined   = combined
  )
}

# 2ND FUNCTION
# Running the SAC 
run_station_sac <- function(PA_combined_table,
                            year_min = 2015, # cut off from this study
                            exclude_stations = "FS", # Excluding FS survey station as not relevant to this study
                            q = 0,
                            nboot = 50,
                            endpoint_multiplier = 2) {
  
  meta_cols <- c("sample_id", "GrabSite", "Date", "Latitude", "Longitude",
                 "Depth", "Sediment_type", "Median_phi", "Mean_phi", "Year",
                 "GrabSite_station", "GrabSite_base", "GrabNumber", "year", "period")
  
  species_cols <- setdiff(names(PA_combined_table), meta_cols)
  
  #1. Filter to years with consistent station design, drop excluded/missing stations
  dat <- PA_combined_table%>%
    filter(year >= year_min,
           !is.na(GrabSite_station),
           !GrabSite_station %in% exclude_stations)
  
  #2. Build one incidence matrix (species x grabs, 0/1) per station
  stations <- unique(dat$GrabSite_station)
  
  incidence_data <- lapply(stations, function(st) {
    sub <- dat %>% filter(GrabSite_station == st)
    mat <- as.matrix(sub[, species_cols])
    mat[is.na(mat)] <- 0   #treat missing as absent
    mat <- (mat > 0) * 1  #ensure binary presence-absence
    t(mat) 
  })
  names(incidence_data) <- stations
  
  #3. Dropping any station with fewer than 2 grabs -> for iNEXT to work
  n_grabs <- sapply(incidence_data, ncol)
  if (any(n_grabs < 2)) {
    warning("Dropping stations with < 2 grabs: ",
            paste(names(n_grabs)[n_grabs < 2], collapse = ", "))
    incidence_data <- incidence_data[n_grabs >= 2]
  }
 
  #4. Set extrapolation endpoint
  endpoint <- max(sapply(incidence_data, ncol)) * endpoint_multiplier
  
  #5. Run iNEXT
  inext_result <- iNEXT(incidence_data,
                        q = q, # Hill number - can do to give c(0, 1, 2) richness, shannon and simpson diversity
                        datatype = "incidence_raw", # 0/1 incidence counts
                        endpoint = endpoint,
                        nboot = nboot)
  
  #6. Plot
  sac_plot <- ggiNEXT(inext_result, type = 1) +
    theme_bw() +
    labs(x = "Number of grabs", y = "Species richness")
  
  list(incidence_data = incidence_data,
       inext_result = inext_result,
       plot = sac_plot)
}