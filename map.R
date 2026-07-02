# Creating a map of grabsites over the years and the new boundaries of the MPA fishing restrictions 

library(sf)        # For handling vector data
library(terra)     # For handling raster data
library(ggplot2)   # For data visualization
library(mapview)   # For interactive mapping
library(tmap)

df_map <- read_excel("data/grab_data/grabsites.xlsx") %>%
  mutate(
    Year = case_when(
      # Already a 4-digit year (e.g. "2022")
      str_detect(as.character(Date), "^\\d{4}$") ~ as.integer(Date),
      # Excel serial number (e.g. 40980)
      is.numeric(Date) ~ year(as.Date(as.numeric(Date), origin = "1899-12-30")),
      # Full date string (e.g. "12/03/2012")
      TRUE ~ year(dmy(as.character(Date), quiet = TRUE))
    )
  ) %>%
  select(-Date)  # drop original Date column, keep Year

# Transform the dataframe into a spatial object.
# Use EPSG:4326 as the Coordinate Reference System (CRS)
df_map_plot <- st_as_sf(df_map, coords = c("Longitude", "Latitude"), crs = 4326)

# Import the boundary of the study area
mpa_restrictions <- st_read("data/Fishing measures relevant to South Arran/Fishing measures relevant to South Arran.shp")

# setting year as a factor
df_map_plot$Year <- as.factor(df_map_plot$Year)

# filtering restrictions for the south arran area
south_arran <- mpa_restrictions %>%
  dplyr::filter(
    str_detect(area_name, "South Arran"))

# how many grab sites fall inside each management zone
grab_zone <- st_join(
  df_map_plot,
  south_arran["area_name"])

table(grab_zone$area_name, useNA = "ifany")

south_arran <- south_arran %>%
  mutate(zone = case_when(
      str_detect(area_name, "No Take Zone") ~ "No-Take Zone",
      str_detect(area_name, "excepted area") ~ "Trawl permitted subject to conditions",
      TRUE ~ "Demersal trawl and dredge prohibited" ))

# Visualize the spatial object
tm_shape(south_arran) +
  tm_polygons(
    fill = "zone",
    fill.scale = tm_scale_categorical(
      values = c(
        "No-Take Zone" = "red",
        "Demersal trawl and dredge prohibited" = "orange",
        "Trawl permitted subject to conditions" = "seagreen" )),
    fill_alpha = 0.5,
    col = "black"
  ) +
  tm_shape(df_map_plot) +
  tm_dots(
    fill = "Year",
    fill.scale = tm_scale_categorical(
      values = c(
        "2012" = "cyan",
        "2013" = "deepskyblue",
        "2015" = "blue3",
        "2022"= "maroon1")),
    size = 0.2
  ) +
  tm_layout(legend.outside = TRUE)
