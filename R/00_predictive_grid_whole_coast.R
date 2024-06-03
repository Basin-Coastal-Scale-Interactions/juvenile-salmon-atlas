library(tidyverse)
library(ggplot2)
library(ggspatial)
library(ggsidekick)
library(sf)
library(here)

##################################
### Whole coast predictive grid from Kelsey
##################################

# loading whole grid
grid_list <- read_csv(here::here("data-raw", "spatial", "coastwide_grid_20231130.csv")) %>% 
  mutate(utm_x_1000 = X / 1000,
         utm_y_1000 = Y / 1000,
         depth = -elevation_meter,
         dist_to_coast_km = coast_distance_meter / 1000) %>%
  

# plotting grid
ggplot(grid_list) + geom_raster(aes(X, Y, fill = depth))

range(grid_list$lat)
range(grid_list$lon)

options(scipen = 999)
# loading polygon with outline for subsetting grid and adding coordinate labels
subpoly <- read_csv(here::here("data-raw", "subset-polygon-nosog-offshore-synoptic.csv")) %>%
  mutate(label = paste(X,Y, sep =","))
options(scipen = 0)

# Plotting  polygon for subsetting predictive grid
ggplot(subpoly) +
  geom_point(aes(X, Y)) + geom_label(aes(X,Y, label = label), nudge_y = 40000)


# converting to spatial object and spatial polygon, respectively
gridsp <- st_as_sf(grid_list, coords = c("X", "Y"), 
                   crs =32609, #26909, # UTM Zone 9M, sp::CRS("+proj=utm +zone=9 +units=m"),
                   remove = FALSE)
polysp  <- st_polygon(list(as.matrix(subpoly[, -3]))) # removing third column with labels

# subsetting grid using st_intersects() with polygon
gridsp_sub <- gridsp[apply(st_intersects(gridsp, polysp), 1, any), ] %>%
  mutate(label = paste(X,Y, sep =",")) 

dim(gridsp)
dim(gridsp_sub)

# Loading chinook data and sourcing coast files for plotting
source("R/plot_map.R")
chinook_dat <- readRDS(here::here("data", "cleaned_index_dat_apr_2024.rds")) %>%
  filter(species == "chinook")


# plotting subsetted grid and chinook points
gridsp_sub %>%
  ggplot() +
  theme_sleek() +
  coord_fixed() +
  geom_sf(data = all_coast, color = "gray80", fill = "gray90")+
  geom_raster(aes(X, Y, fill = depth), alpha = 0.6) +
  geom_point(
    data = filter(chinook_dat, n_juv == 0),
    aes(utm_x, utm_y, pch = as.factor(ex)),
    pch = 4,
    col = "grey50",
    size = 0.5, alpha = 0.6
  ) +
  geom_point(
    data = filter(chinook_dat, n_juv > 0),
    aes(utm_x, utm_y,
        size = n_juv,
    ), # fill = pos_pt_fill,
    col = "red", pch = 21, alpha = 0.6
  ) +
  # xlim(c(650000, 710000)) +
  # ylim(5470000, 5530000)+
  annotation_scale(location = "bl", height = unit(0.15, "cm")) +
  scale_fill_continuous(trans = 'reverse') +
  #   scale_y_continuous(breaks = seq(5200000, 6100000, by = 50000)) +
  # scale_x_continuous(breaks = seq(0, 1000000, by = 50000)) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #   geom_polygon(data = subpoly, aes(X, Y), alpha = 0.3) +
  #xlim(c(600000,750000)) + ylim(c(5600000, 5800000)) +
  #xlim(c(250000,350000)) + ylim(c(5950000, 6050000)) +
  geom_label(data = subpoly, aes(X,Y, label = label), size = 3) +
  NULL 
ggsave(here("figs", "grid-subset-polygon.png"), width = 9, height =8)

saveRDS(st_drop_geometry(gridsp_sub), 
        file = here("data", "pred_grid_all_coast.rds"))
