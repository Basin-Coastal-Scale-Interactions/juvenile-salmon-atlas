library(tidyverse)
library(patchwork)

pred_cmb <- readRDS(here("data", "pred_cmb.rds"))

grid_df <- pred_cmb %>%
  mutate(month_f = month(month, label = TRUE),
         #Xr = x * 1000, Yr = y * 1000,
         region = fct_recode(region,
                             "Northern BC/AK" = "NBC/SEAK"           
         )) %>%
  mutate(region = fct_relevel(region,
                              c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI",
                                "Col. Coastal", "Col. North", "Fraser Yearling", 
                                "Fraser Summer 4.1", "Upper Col. Yearling")
  )) %>%
  dplyr::select(Xr, Yr, pred_cmb, salmon_region = region, month_f)


prot_1 <- grid_df %>%
  filter(salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       labels = scales::comma,
                       limits = c(0, max(grid_df$pred_cmb)),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        legend.position="none",
        plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution") +
  NULL

prot_2 <- grid_df %>%
  filter(!salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nabundance\nof juvenile\nchinook\nper tow",
                       labels = scales::comma,
                       limits = c(0, max(grid_df$pred_cmb)),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        #legend.margin=margin(t = 0, b=0, unit='cm'),
        plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  # ggtitle("Juvenile chinoook stock distribution") +
  NULL

# This produces a nice output using the patchwork pkg
ggsave(prot_1 + prot_2,
       filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_rotated_no_cv_cutoff.png"),  
       width = 15, height = 11)


### Relative abundance plots

grid_rel <- grid_df %>%
  group_by(salmon_region) %>%
  mutate(rel_pred = pred_cmb/max(pred_cmb)) %>%
  dplyr::select(Xr, Yr, rel_pred, salmon_region, month_f)

prel_1 <- grid_rel %>%
  filter(salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
  #filter(salmon_region %in% c("SoG Coastal", "WA/OR Coastal"), month_f %in% c("Jun", "Jul")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Relative\nabundance of\njuvenile\nchinook\nper tow",
                       labels = scales::comma,
                       limits = c(0, 1),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        legend.position="none",
        plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution")+
  NULL
#prel_1

prel_2 <- grid_rel %>%
  filter(!salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
  #filter(salmon_region %in% c("Central BC", "WCVI"), month_f %in% c("Jun", "Jul")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nrelative\nabundance\nof juvenile\nchinook\nper tow",
                       labels = scales::comma,
                       limits = c(0, 1),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        plot.margin=margin(b = 0, t = 0, l = 0, r = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  # ggtitle("Juvenile chinoook stock distribution") +
  NULL

# This produces a nice output using the patchwork pkg
ggsave(prel_1 + prel_2,
       filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_rotated_no_cv_cutoff.png"), 
       width = 15, height = 11)
