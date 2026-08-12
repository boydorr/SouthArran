# Creating a map of grabsites over the years and the new boundaries of the MPA fishing restrictions 

library(sf)        # For handling vector data
library(terra)     # For handling raster data
library(ggplot2)   # For data visualization
library(mapview)   # For interactive mapping
library(tmap)
library(leaflet)

# Creating a map of grabsites over the years and the new boundaries of the MPA fishing restrictions 
# Getting Grab site from grab survey data 
df_map <- read_excel("data/grab_data/grabsites.xlsx") %>%
  filter(!is.na(GrabSite)) %>%                      # drop blank row
  mutate(
    Longitude = as.numeric(Longitude),               # convert to numeric FIRST
    Year = case_when(
      str_detect(as.character(Date), "^\\d{4}$") ~ as.integer(Date),
      is.numeric(Date) ~ year(as.Date(as.numeric(Date), origin = "1899-12-30")),
      TRUE ~ year(dmy(as.character(Date), quiet = TRUE))
    )
  ) %>%
  filter(Longitude < -3, Longitude > -8) %>%
  select(-Date)

df_map_plot <- st_as_sf(df_map, coords = c("Longitude", "Latitude"), crs = 4326)
df_map_plot$Year <- as.factor(df_map_plot$Year)

# Import the boundary of the study area
mpa_restrictions <- st_read("data/Fishing measures relevant to South Arran/Fishing measures relevant to South Arran.shp")

# Import survey stations in study area
mpa_surveystations16 <- st_read("data/Fishing measures relevant to South Arran/Survey boxes - 2021 and 2016/Arran_survey_boxes_2016.shp")
mpa_surveystations22 <- st_read("data/Fishing measures relevant to South Arran/Survey boxes - 2021 and 2016/Arran_survey_boxes_2021_defined.shp")

# filtering restrictions for the south arran area
south_arran <- mpa_restrictions %>%
  dplyr::filter(str_detect(area_name, "South Arran"))

# Transform to match the boundary CRS
points_sf <- st_transform(df_map_plot, st_crs(south_arran))

# Remove FS before transforming/using stations at all
stations_sf <- mpa_surveystations22 %>%
  filter(Name != "FS") %>%
  st_transform(st_crs(south_arran))

# How many grab sites fall inside each management zone?
grab_zone <- st_join(df_map_plot, south_arran["area_name"])
table(grab_zone$area_name, useNA = "ifany")

# Defining the zones of the South Arran MPA from shapefiles
south_arran <- south_arran %>%
  mutate(zone = case_when(
    str_detect(area_name, "No Take Zone") ~ "No-Take Zone",
    str_detect(area_name, "excepted area") ~ "Trawl permitted subject to conditions",
    TRUE ~ "Demersal trawl and dredge prohibited"
  ))

# Static map
ggplot() +
  geom_sf(data = south_arran, aes(fill = zone), colour = "black", alpha = 0.5) +
  scale_fill_manual(
    values = c(
      "No-Take Zone" = "red",
      "Demersal trawl and dredge prohibited" = "orange",
      "Trawl permitted subject to conditions" = "seagreen")
  ) +
  geom_sf(data = points_sf, aes(colour = factor(Year)), size = 0.5) +
  scale_colour_manual(
    values = c(
      "2012" = "cyan",
      "2013" = "deepskyblue",
      "2015" = "blue3",
      "2022" = "maroon1")
  ) +
  geom_sf(data = stations_sf, colour = "black", fill = NA, linewidth = 0.3) +
  labs(title = "Sampling Locations", colour = "Year") +
  theme_minimal()

# Interactive map
map_zones <- mapview(
  south_arran,
  zcol = "zone",
  col.regions = c(
    "No-Take Zone" = "red",
    "Demersal trawl and dredge prohibited" = "orange",
    "Trawl permitted subject to conditions" = "seagreen"
  ),
  alpha.regions = 0.5,
  layer.name = "MPA Zone"
)

map_points <- mapview(
  points_sf,
  zcol = "Year",
  col.regions = c(
    "2012" = "cyan",
    "2013" = "deepskyblue",
    "2015" = "blue3",
    "2022" = "maroon1"
  ),
  cex = 2,
  layer.name = "Grab Site Year"
)

map_stations <- mapview(
  stations_sf,
  alpha.regions = 0,
  color = "black",
  lwd = 1,
  layer.name = "Survey Stations"
)

# Centroids for labels (FS already excluded upstream)
stations_centroids <- stations_sf %>%
  st_transform(4326) %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  setNames(c("lon", "lat")) %>%
  bind_cols(Name = stations_sf$Name)

combined_map <- map_zones + map_points + map_stations

combined_map@map <- combined_map@map %>%
  addLabelOnlyMarkers(
    data = stations_centroids,
    lng = ~lon, lat = ~lat,
    label = ~Name,
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "centre",
      offset = c(-10, 0),
      textOnly = TRUE,
      style = list("font-size" = "10px", "font-weight" = "bold")
    )
  )

combined_map