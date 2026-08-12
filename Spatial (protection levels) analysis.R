# Adding protection levels - Spatial analysis 
library(sf)
library(dplyr)
library(stringr)

# Getting one row per unique grab site with its coordinates
# Assign spatial protection zones (High/Medium/Low) to grab sites based on an MPA management-zone shapefile, and attach the result to grab dataframe by site
add_protection_zones <- function(df_env,
                                 mpa_shapefile,
                                 df_sites,
                                 site_col = "GrabSite",
                                 lon_col = "Longitude",
                                 lat_col = "Latitude",
                                 area_col = "area_name",
                                 area_filter = "South Arran",
                                 crs = 4326,
                                 high_pattern = "No Take Zone",
                                 low_pattern = "excepted area",
                                 outside_label = "Outside MPA") {
  
  # getting one row per unique grab site with its coordinates
  site_locations <- df_sites %>%
    distinct(.data[[site_col]], .data[[lon_col]], .data[[lat_col]]) %>%
    st_as_sf(coords = c(lon_col, lat_col), crs = crs)
  
  #Read in shapefile and classify zones
  mpa <- if (inherits(mpa_shapefile, "sf")) {
    mpa_shapefile
  } else {
    st_read(mpa_shapefile, quiet = TRUE)
  }
  
  mpa_zones <- mpa %>%
    filter(str_detect(.data[[area_col]], area_filter)) %>%
    mutate(protection_level = case_when(
      str_detect(.data[[area_col]], high_pattern) ~ "No Take Zone",
      str_detect(.data[[area_col]], low_pattern)  ~ "Dredge Prohibited, Trawl permitted subject to conditions",
      TRUE                                         ~ "Demersal Trawl and Dredge Prohibited"
    ))
  
  # Assign each grabs its protection level
  site_zones <- st_join(site_locations, mpa_zones["protection_level"]) %>%
    st_drop_geometry() %>%
    distinct(.data[[site_col]], .keep_all = TRUE) %>%   # guard against boundary overlaps
    mutate(protection_level = as.character(protection_level),
           protection_level = if_else(is.na(protection_level), outside_label, protection_level),
           protection_level = factor(protection_level,
                                     levels = c("Dredge Prohibited", "Demersal Trawl and Dredge Prohibited", "No Take Zone", outside_label))) %>%
    select(all_of(site_col), protection_level)
  
  #Attach to df_env
  df_env <- df_env %>%
    left_join(site_zones, by = site_col)
  
  df_env
}


# Changing Depth and phi size to categorical levels
add_zones_protect_depth_sedi <- function(df_env,
                                 mpa_shapefile,
                                 df_sites,
                                 site_col = "GrabSite",
                                 lon_col = "Longitude",
                                 lat_col = "Latitude",
                                 area_col = "area_name",
                                 area_filter = "South Arran",
                                 crs = 4326,
                                 high_pattern = "No Take Zone",
                                 low_pattern = "excepted area",
                                 outside_label = "Outside MPA",
                                 depth_col = "Depth",
                                 sediment_col = "Sediment_type") {
  
  # getting one row per unique grab site with its coordinates
  site_locations <- df_sites %>%
    distinct(.data[[site_col]], .data[[lon_col]], .data[[lat_col]]) %>%
    st_as_sf(coords = c(lon_col, lat_col), crs = crs)
  
  # Read in shapefile and classify zones
  mpa <- if (inherits(mpa_shapefile, "sf")) {
    mpa_shapefile
  } else {
    st_read(mpa_shapefile, quiet = TRUE)
  }
  
  mpa_zones <- mpa %>%
    filter(str_detect(.data[[area_col]], area_filter)) %>%
    mutate(protection_level = case_when(
      str_detect(.data[[area_col]], high_pattern) ~ "No Take Zone",
      str_detect(.data[[area_col]], low_pattern)  ~ "Dredge Prohibited, Trawl permitted subject to conditions",
      TRUE                                         ~ "Demersal Trawl and Dredge Prohibited"
    ))
  
  # Assign each grab its protection level
  site_zones <- st_join(site_locations, mpa_zones["protection_level"]) %>%
    st_drop_geometry() %>%
    distinct(.data[[site_col]], .keep_all = TRUE) %>%   # guard against boundary overlaps
    mutate(protection_level = as.character(protection_level),
           protection_level = if_else(is.na(protection_level), outside_label, protection_level),
           protection_level = factor(protection_level,
                                     levels = c("Dredge Prohibited, Trawl permitted subject to conditions", "Demersal Trawl and Dredge Prohibited", "No Take Zone", outside_label))) %>%
    select(all_of(site_col), protection_level)
  
  # Attach to df_env
  df_env <- df_env %>%
    left_join(site_zones, by = site_col)
  
  # Derive depth category from Depth
  df_env <- df_env %>%
    mutate(
      depth_category = case_when(
        is.na(.data[[depth_col]])                          ~ NA_character_,
        .data[[depth_col]] <= 27                            ~ "Shallow",
        .data[[depth_col]] > 27  & .data[[depth_col]] <= 90 ~ "Intermediate",
        .data[[depth_col]] > 90                              ~ "Deep",
        TRUE                                                  ~ NA_character_
      ),
      depth_category = factor(depth_category,
                              levels = c("Shallow", "Intermediate", "Deep"))
    )
  
  # Derive simplified sediment type from Sediment_type
  df_env <- df_env %>%
    mutate(
      sediment_type_simple = case_when(
        is.na(.data[[sediment_col]])                          ~ NA_character_,
        str_detect(str_to_lower(.data[[sediment_col]]), "silt")   ~ "Silt",
        str_detect(str_to_lower(.data[[sediment_col]]), "sand")   ~ "Sand",
        str_detect(str_to_lower(.data[[sediment_col]]), "gravel") ~ "Gravel",
        TRUE                                                       ~ NA_character_
      ),
      sediment_type_simple = factor(sediment_type_simple,
                                    levels = c("Silt", "Sand", "Gravel"))
    )
  
  df_env
}
