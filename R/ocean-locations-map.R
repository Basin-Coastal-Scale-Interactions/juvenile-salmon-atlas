### Plots of data summary

library(tidyverse)
library(ggrepel)
library(sf)
library(here)
library(rnaturalearth)

# needed to add depth raster
# grid_list <- read_csv(here::here("data-raw", "spatial", "grid_extended20241211.csv")) %>% 
#   mutate(utm_x_1000 = X / 1000,
#          utm_y_1000 = Y / 1000,
#          depth = -elevation_meter,
#          dist_to_coast_km = coast_distance_meter / 1000)  %>%
#   select(X, Y, lon, lat, depth, utm_x_1000, utm_y_1000)
# ggplot(grid_list) + geom_raster(aes(X, Y, fill = depth))


min_lat <- 46.75
max_lat <-  57.248
min_lon <- -136.25
max_lon <- -122.08

coast <- rbind(ne_states( "United States of America",
                                         returnclass = "sf"),
               ne_states( "Canada", returnclass = "sf")) %>%
  st_crop(., xmin = min_lon, ymin = min_lat-2, xmax = max_lon+10, ymax = max_lat) %>%
  #  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))


waterlabels <- read_csv("data-raw/map-locations.csv")  |>
  bind_rows(read_csv("data-raw/islands.csv"))

pwater <- waterlabels |>
 # filter(lon > -136) |>
ggplot() +
 # geom_tile(aes(x = X, y = Y, fill = depth)) +
  geom_sf(data = coast, color = "gray70", fill = "gray90",  size = 1.25) +
  #scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank(),
        panel.background = element_rect( fill = "lightblue")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_viridis_c(direction = -1, option = "G") +
  # theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8),
  #       legend.title = element_text(size=15),
  #       legend.text = element_text(size=13)) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15),
           xlim = c(min_lon + 1.7, max_lon + 0.15)) +
  geom_point(data = waterlabels, aes(x, y),size = waterlabels$pointsize) +
  #geom_text(data = waterlabels, aes(x, y, label = name),size = 4, fontface = "italic") +
 annotate("text", x = -132.5, y = 48, label = "Pacific\nOcean", fontface = "italic", size = 7) +
  annotate("text", x = -124, y = 53, label = "British\nColumbia", family = "sans", size = 6) +
  annotate("text", x = -131.5, y = 55.9, label = "Alaska", family = "sans", size = 6) +
  annotate("text", x = -125.2, y = 47.5, label = "Washington", family = "sans", size = 6) +
  # annotate("text", x = -132.2, y = 53.7, label = "Haida\nGwaii", family = "sans", size = 5) +
  # annotate("text", x = -125.9, y = 49.7, label = "Vancouver\nIsland", family = "sans", size = 5) +
  geom_text_repel(aes(x, y, label = name,  fontface = fontface), size = waterlabels$size,
                  min.segment.length = 0.2, force_pull = 1, force = 4, 
                  box.padding = 0.25,
                  verbose = TRUE, seed = 2) +
  NULL


#pwater

ggsave(pwater,
       filename = here::here("figs", "ocean-locations-map.png"), 
       width = 5.1, height = 6)


waterlabels <- data.frame(
    x = c(-130, -133.51, -132, -130.7, -129.4, -125.5, -123.7, -122.8),
    y = c(49, 56.393333, 54.45, 53, 51.5, 48.2, 47.8, 49.9),
    name = c(
      "Pacific\nOcean",
      "Sumner\nStrait",
      "Dixon Entrance",
      "Hecate\nStrait",
      "Queen\nCharlotte\nStrait",
      "Juan de\nFuca\nStrait",
      "Puget\nSound",
      "Strait of\nGeorgia"
    ),
    fontface = rep("italic", 8),
    size = c(6, 3, 3, 3, 3, 3, 3, 3)
  )  
  
write_csv(waterlabels, "data-raw/map-locations.csv")  
