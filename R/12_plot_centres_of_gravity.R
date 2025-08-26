library(tidyverse)
library(patchwork)
library(here)

pred_cmb <- readRDS(here("data", "pred_cmb.rds"))

rotated_coast <- readRDS(here("data", "rotated_coast_outline.rds"))
source("R/plot_map.R")

grid_df <- pred_cmb %>%
  mutate(month_f = month(month, label = TRUE),
         #Xr = x * 1000, Yr = y * 1000,
         region = fct_recode(region,
                             "Northern BC/AK" = "NBC/SEAK"           
         )) %>%
  mutate(region = fct_relevel(region,
                              c("Northern BC/AK", "Central BC", "WCVI", "Salish Sea", 
                                "Fraser Fall 4.1",
                                "Fraser Summer 4.1",
                                "Fraser Yearling",
                                "WA/OR Coastal", 
                                "Lower Col.", 
                                "Upper Col. Subyearling",
                                "Upper Col. Yearling"
                              )
  )) %>%
  dplyr::select(lon, lat, utm_x_1000, utm_y_1000, Xr, Yr, pred_cmb, salmon_region = region, month_f, month)

filter(grid_df, is.na(Xr) | is.na(Yr))



grid_cog <- grid_df |>
  group_by(salmon_region, month_f) |>
  mutate(temp = (utm_x_1000  * pred_cmb ) / sum(pred_cmb )) |>
  mutate(temp2 = (utm_y_1000  * pred_cmb ) / sum(pred_cmb )) |>
  summarise(mean_lon = sum(temp), mean_lat = sum(temp2), month = unique(month)) 


grid_rel_cog <- grid_df %>%
  group_by(salmon_region) %>%
  mutate(rel_pred = pred_cmb/max(pred_cmb)) %>%
  ungroup() |>
  group_by(salmon_region, month_f) |>
  mutate(temp = (utm_x_1000  * pred_cmb ) / sum(pred_cmb )) |>
  mutate(temp2 = (utm_y_1000  * pred_cmb ) / sum(pred_cmb )) |>
  summarise(mean_lon = sum(temp), mean_lat = sum(temp2), month = unique(month)) 

grid_cog |>
  ungroup() |>
  ggplot() +
  geom_sf(data = all_coast_km, color = "gray83") +
  geom_path(aes(mean_lon, mean_lat, color = factor(month_f), group = 1), linewidth = 1.4)  +
 # geom_point(aes(mean_lon, mean_lat, fill = month_f)) +
  scale_x_continuous(name = NULL, limits = range(grid_df$utm_x_1000)+c(20,-10), expand = c(0, 0),
                     breaks = seq(-134, -122, by = 4)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$utm_y_1000)+c(100,-60), expand = c(0, 0)) +
  ggsidekick::theme_sleek() +
  scale_colour_viridis_d(name = "Month",
    alpha = 1,
    begin = 0,
    end = 1,
    direction = -1,
  ) +
# scale_fill_viridis_d(month(5:12, label = T, abbr = T), direction = -1) +
  facet_wrap(~salmon_region) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.87, 0.15),
        legend.key.size = unit(0.5,"cm"),
        legend.text =  element_text(size = 10),
        legend.title = element_text(size = 14))

### rotated CoG

grid_cogr <- grid_df |>
  group_by(salmon_region, month_f) |>
  mutate(temp = (Xr  * pred_cmb ) / sum(pred_cmb )) |>
  mutate(temp2 = (Yr  * pred_cmb ) / sum(pred_cmb )) |>
  summarise(mean_xr = sum(temp), mean_yr = sum(temp2), month = unique(month)) 

grid_cogr |>
  ungroup() |>
  ggplot() +
  geom_sf(data = rotated_coast, color = "gray83") +
  geom_path(aes(mean_xr, mean_yr, color = factor(month_f), group = 1), linewidth = 1.4)  +
  # geom_point(aes(mean_lon, mean_lat, fill = month_f)) +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(20,-10), expand = c(0, 0),
                     breaks = seq(-134, -122, by = 4)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(100,-60), expand = c(0, 0)) +
  ggsidekick::theme_sleek() +
  scale_colour_viridis_d(name = "Month",
                         alpha = 1,
                         begin = 0,
                         end = 1,
                         direction = -1,
  ) +
  # scale_fill_viridis_d(month(5:12, label = T, abbr = T), direction = -1) +
  facet_wrap(~salmon_region, nrow =2) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.25),
        legend.key.size = unit(0.5,"cm"),
        legend.text =  element_text(size = 10),
        legend.title = element_text(size = 14))
