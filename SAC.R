# Species accumulation curves 
library(iNEXT)
library(dplyr)
library(tidyverse)

# Running the SAC 
run_station_sac <- function(PA_combined_table,
                            year_min = 2015, # cut off from this study
                            exclude_stations = "FS", 
                            station_prefix = NULL,
                            split_by_year = TRUE,
                            q = 0, # hill number (species richness)
                            nboot = 50,
                            endpoint_multiplier = 2) {
  #'Cleaning presence and absence data for iNEXT to create a species accumulation curve (SAC) based on survey stations
  #'@description A function that uses previously made presence and absence dataframe to create SAC to assess survey effort (number of grabs)
  #'
  #'@param PA_combined_table The presence and absence data
  #'@param year_min The cut off year (in this study it is because of when survey stations were set)
  #'@param exclude_stations Which stations to not include in SAC (Excluding FS survey station as not relevant to this study)
  #'@param station_prefix Prefix of survey station in what to restrict by 
  #'@param split_by_year (examine survey station by year or across the whole study period (=FALSE))
  #'@param q The order of hill number (0 = species richness)
  #'@param nboot 
  #'@param endpoint_multiplier 
  #'
  #'@return A species accumulation curve based on indicidence data by survey stations
  
  meta_cols <- c("sample_id", "GrabSite", "Date", "Latitude", "Longitude",
                 "Depth", "Sediment_type", "Median_phi", "Mean_phi", "Year",
                 "GrabSite_station", "GrabSite_base", "GrabNumber", "year", "period",
                 "GrabSite_number", "protection_level", "depth_category", "sediment_type_simple")
  
  species_cols <- setdiff(names(PA_combined_table), meta_cols)
  
  #1. Filter to years with consistent station design, drop excluded/missing stations
  dat <- PA_combined_table%>%
    filter(year >= year_min,
           !is.na(GrabSite_station),
           !GrabSite_station %in% exclude_stations)
  
  # Restricting by survey station
  if (!is.null(station_prefix)) {
    dat <- dat %>% filter(str_starts(GrabSite_station, station_prefix))
  }
  
  if (nrow(dat) == 0) {
    stop("No rows remain after filtering - check station_prefix matches your GrabSite_station values.")
  }
  
  #2. Build one incidence matrix
      #(species x grabs, 0/1) per station
  stations <- unique(dat$GrabSite_station)
  
  dat <- dat %>%
    mutate(group_id = if (split_by_year) paste(GrabSite_station, year, sep = "_") else GrabSite_station)
  
  groups <- unique(dat$group_id)
  
  incidence_data <- lapply(groups, function(g) {
    sub <- dat %>% filter(group_id == g)
    mat <- as.matrix(sub[, species_cols])
    mat[is.na(mat)] <- 0
    mat <- (mat > 0) * 1
    t(mat)
  })
  names(incidence_data) <- groups
  
  #3. Dropping any station with fewer than 2 grabs or zero observed species -> for iNEXT to work
  keep <- sapply(incidence_data, function(m) {
    ncol(m) >= 2 && sum(rowSums(m) > 0) > 0
  })
  
  if (any(!keep)) { # telling what groups are being dropped
    warning(
      "Dropping groups: ",
      paste(names(incidence_data)[!keep], collapse = ", "))
  }
  
  incidence_data <- incidence_data[keep]
 
  #4. Set extrapolation endpoint
  endpoint <- max(sapply(incidence_data, ncol)) * endpoint_multiplier # making a comparable x axis limit
  
  #5. Run iNEXT
  inext_result <- iNEXT(incidence_data,
                        q = q, # Hill number, set to q to change if q value changes (eg Shannon D(q=1))
                        datatype = "incidence_raw", # 0/1 incidence counts
                        endpoint = endpoint,
                        nboot = nboot)
  
  #6. Plot
  sac_plot <- ggiNEXT(inext_result, type = 1) +
    theme_bw() +
    labs(x = "Number of grabs", y = "Species richness")
  
  station_year_colours <- c(
      # T sites 
      "T1_2015" = "#A6CEE3",
      "T1_2022" = "#1F78B4",
      "T2_2015" = "#8DD3C7",
      "T2_2022" = "#008080",
      "T3_2015" = "#B2DF8A",
      "T3_2022" = "#33A02C",
      "T4_2015" = "#9ECAE1",
      "T4_2022" = "#08519C",
      # D sites 
      "D1_2015" = "#CAB2D6",
      "D1_2022" = "#6A3D9A",
      "D3_2015" = "#FCCDE5",
      "D3_2022" = "#C51B8A",
      "D4_2015" = "#DADAEB",
      "D4_2022" = "#54278F",
      "D5_2015" = "#FBB4AE",
      "D5_2022" = "#E31A1C",
      "D6_2015" = "#E0C2FF",
      "D6_2022" = "#7B3294")
  
  sac_plot <- sac_plot +
    scale_colour_manual(values = station_year_colours) +
    scale_fill_manual(values = station_year_colours)
  
  sac_plot$data <- sac_plot$data %>%
    mutate(station_type = if_else(str_starts(Assemblage, "D"), "D sites", "T sites"))
  
  sac_plot <- sac_plot + facet_wrap(~ station_type, scales = "free_x")
  
  list(incidence_data = incidence_data,
       inext_result = inext_result,
       plot = sac_plot)
}