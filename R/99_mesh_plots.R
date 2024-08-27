# Diagnostic plots

library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggsidekick)
library(sdmTMB)
library(here)
library(sf)

source(here("R", "plot_map.R"))
mdmesh <- readRDS(here("data", "fits", "dmesh_allcoast.rds")) 
chinook_dat <- readRDS(here("data", "chinook_dat_allcoast.rds"))
gsimesh <- readRDS("data/gsi-prop-mesh.rds")

index_grid <- readRDS(here("data", "pred_grid_all_coast.rds")) %>%
  mutate(year = 2019L, 
         year_f = as.factor(year),
         survey_f = as.factor("hss"),
         target_depth = 0,
         day_night = "DAY",
         scale_depth = -0.5427145)

all_coast_mesh <- rbind(rnaturalearth::ne_states( "United States of America", 
                                                returnclass = "sf"), 
                      rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -139, ymin = 43, xmax = -115, ymax = 59) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=km")) %>%
  sf::st_crop(., xmin=  53.19298, ymin= 4927.892, xmax = 1200, ymax= 6550)


png(here::here("figs", "mesh-comp.png"), width = 9, height = 4, units = "in", res = 300)

par(mfrow =c(1,2))
par(mar = c(0,0,1,0))

plot(sf::st_geometry(all_coast_mesh), col = "#DDDDDD80", lwd = 1, asp = 1, border = "#BBBBBB", main = "A) Species-wide model")
plot(mdmesh$spde$mesh, edge.color = "blue", col = NA, lwd = 1,add = TRUE)
mtext(paste0(mdmesh$spde$mesh$n, " knots"), line =-2, adj = 0.6)

plot(sf::st_geometry(all_coast_mesh), col = "#DDDDDD80", lwd = 1, asp = 1, border = "#BBBBBB", main = "B) Stock-specific model")
plot(gsimesh$mesh, edge.color = "blue", col = NA, lwd = 1,add = TRUE)
mtext(paste0(gsimesh$mesh$n, " knots"), line =-2, adj = 0.6)

dev.off()
