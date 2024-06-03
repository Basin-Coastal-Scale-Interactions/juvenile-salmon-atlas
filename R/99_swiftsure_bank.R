### Maps of juvenile salmon abundance for Swiftsure Bank, BC

library(sdmTMB)
library(sdmTMBextra)
library(tidyverse)
library(ggsidekick)
library(here)
library(sf)

source(here("R", "plot_map.R"))

dat_in <- readRDS(here::here("data", "catch_survey_sbc_withnorth.rds")) 

dat_tbl_dmesh <- readRDS(here("data", "fits", "all_dmesh_allcoast_withpred.rds"))

# Swiftsure Bank area GIS layer
swsu <- read_sf("../gis-layers/swiftsure/JdFEddy_EBSA_Sept6.shp") %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>%
  st_zm()

# US-Canada border layer
usca <- read_sf("../gis-layers/US_CAN_border/US_CAN_border.shp") %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>%
  st_zm()

# US-Canada EEZ layer
eez <- read_sf("../gis-layers/BC_EEZ/BC_EEZ.shp") %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>%
  st_zm()

ggplot() + geom_sf(data = eez)

library(rnaturalearth)
library(rnaturalearthhires)

# rough WA coastline
wa_coast <- rnaturalearth::ne_states( "United States of America", 
                                      returnclass = "sf") %>% 
  sf::st_crop(., xmin = -135, ymin = 46.25, xmax = -122.25, ymax = 55.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# Canada-US polygon
caus <- read_sf("../gis-layers/Canada_US/Canada_US.shp") %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>%
  st_zm()

ggplot() + geom_sf(data = caus)

caus_c <- st_crop(caus, ymin = 5300000, ymax = 5500000, xmin = 700000, xmax = 920000)
ggplot() + geom_sf(data = caus_c)


# high water line layer
hwl <- read_sf("../gis-layers/CHS_high_water_line/CHS_high_water_line.shp") %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>%
  st_zm()

hwl_poly <- hwl %>%
  st_cast('POLYGON')

raster::extent(hwl)
ggplot() + geom_sf(data = hwl)

# cropping layer so everything runs faster
hwl_c <- st_crop(hwl, ymin = 5346995, ymax = 5500000, xmin = 700000, xmax = 961697.2)
ggplot() + geom_sf(data = hwl_c)




# zooming into Swiftsure Bank for chinook
species <- dat_tbl_dmesh$species[2]

x_range <- pm(794406, 50000)
y_range <- pm(5396664, 50000)
 
pred_grid_filtered <- dat_tbl_dmesh$pred_grid[[2]] %>%
  filter(X >= x_range[1] & X <= x_range[2]  &
           Y >= y_range[1] & Y <= y_range[2]) %>%
  filter(month_f %in% c("Jul", "Oct"))

plot_map(pred_grid_filtered, exp(est), 
         area = "BC", show_coast = FALSE, show_raw_data = FALSE) + 
  annotation_scale(location = "tr", height = unit(0.1, "cm")) +
  ggtitle(paste0("Predicted distribution of juvenile ", species," around Swiftsure Bank")) +
  # added new plot ranges based on UTM
  scale_x_continuous(name = NULL, limits = pm(794406, 50000), expand = c(0, 0), 
                     breaks = seq(-124.6, -125.6, by = -0.2)) +
  scale_y_continuous(name = NULL, limits = pm(5396664, 50000), expand = c(0, 0)) +
  theme(legend.key.height = unit(0.06, 'npc')) +
  geom_sf(data = caus_c, color = "gray80", lty = 1, lwd = 0.5, fill = "gray90") +
  #geom_sf(data = wa_coast, color = "gray80", lty = 1, lwd = 0.5, fill = "white") +
  geom_sf(data = eez, color = "gray80", lty = 2, fill = "NA") +
  geom_sf(data = swsu, color = "black", lwd = 1.3, fill = NA) + 
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuveniles\nper tow",
                       labels = comma,
                       trans = fourth_root_power_trans(),
                       limits = c(0, 85)) +
  facet_wrap(~month_f, labeller = as_labeller(c("Jul" = "July", "Oct" = "October"))) +
  NULL

# All species loop
for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
  x_range <- pm(794406, 50000)
  y_range <- pm(5396664, 50000)
  
  pred_grid_filtered_l <- dat_tbl_dmesh$pred_grid[[i]] %>%
    filter(X >= x_range[1] & X <= x_range[2]  &
             Y >= y_range[1] & Y <= y_range[2]) %>%
    filter(month_f %in% c("Jul", "Oct"))
  
  p <- plot_map(pred_grid_filtered_l, exp(est), 
           area = "BC", show_coast = FALSE, show_raw_data = FALSE) + 
    annotation_scale(location = "tr", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," around Swiftsure Bank")) +
    # added new plot ranges based on UTM

    scale_y_continuous(name = NULL, limits = y_range, expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc')) +
    geom_sf(data = caus_c, color = "gray70", lty = 1, fill = "gray90") +
    #geom_sf(data = hwl_c, color = "gray80", lty = 1, lwd = 0.3, fill = "NA") +
    #geom_sf(data = wa_coast, color = "gray80", lty = 1, lwd = 0.3, fill = "white") +
    geom_sf(data = eez, color = "gray80", lty = 2, fill = "NA") +
    geom_sf(data = swsu, color = "black", lwd = 1.3, fill = NA) + 
    scale_fill_viridis_c(name = "Predicted\nnumber of\njuveniles\nper tow",
                         labels = comma,
                         trans = fourth_root_power_trans(),
                         limits = c(0, 85)) +
    facet_wrap(~month_f, labeller = as_labeller(c("Jul" = "July", "Oct" = "October"))) +
    scale_x_continuous(name = NULL, limits = x_range, expand = c(0, 0),
                       breaks = rev(c(-124.5, -125.0, -125.5))) +
    NULL
  
  png(here::here("figs", "swiftsure", paste0("Swiftsure_", species, "_njuv.png")),
      height = 5, width = 8, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off()
  # 
  # p <- plot_map(pred_grid_filtered %>% filter(month_f == "Oct"), 
  #               exp(est), area = "BC", show_coast = FALSE,
  #               show_raw_data = FALSE, dat_tbl_dmesh$data[[i]] %>% filter(month_f == "Oct")) + 
  #   annotation_scale(location = "tr", height = unit(0.1, "cm")) +
  #   ggtitle(paste0("Predicted distribution of juvenile ", species," in October around Swiftsure Bank")) +
  #   # added new plot ranges based on UTM
  #   scale_x_continuous(name = NULL, limits = x_range, expand = c(0, 0)) +
  #   scale_y_continuous(name = NULL, limits = y_range, expand = c(0, 0)) +
  #   theme(legend.key.height = unit(0.06, 'npc')) +
  #   geom_sf(data = hwl_c, color = "gray80", lty = 1, fill = "NA") +
  #   geom_sf(data = wa_coast, color = "gray80", lty = 1, fill = "gray90") +
  #   geom_sf(data = eez, color = "gray80", lty = 2, fill = "NA") +
  #   geom_sf(data = swsu, color = "black", lwd = 1.3, fill = NA) + 
  #   scale_fill_viridis_c(name = "Predicted\nnumber of\njuveniles\nper tow",
  #                        labels = comma,
  #                        trans = fourth_root_power_trans(),
  #                        limits = c(0, 85)) +
  #   NULL
  # 
  # png(here::here("figs", paste0("Swiftsure_", species, "_njuv_Oct.png")),
  #     height = 7, width = 8, units = "in", res = 300)
  # print(p) # need to print() plot inside for loop!
  # dev.off() 
  
  print("Done!")
}

g + scale_x_continuous(name = NULL, limits = pm(794406, 60000)/1000, expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = pm(5396664, 60000)/1000, expand = c(0, 0)) 


ggplot() + geom_sf(data = eez) + facet_wrap(~territory1)
ggplot() + geom_sf(aes(data = eez, color = mrgid))
ggplot() + geom_sf(data = eez) 
#### Swiftsure bank maps
ggplot() + geom_sf(data = usca) 

ggplot() + geom_sf(data = hwl) + facet_wrap(~REGION)

ggplot() + geom_sf(data = hwl) + facet_wrap(~REGION)


library(leaflet)
library(leaflet.esri)

leaflet() %>%
  addEsriBasemapLayer(esriBasemapLayers$Gray, autoLabels = F)
