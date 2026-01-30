library(tidyverse)
library(patchwork)
library(here)

pred_cmb <- readRDS(here("data", "pred_cmb.rds"))

rotated_coast <- readRDS(here("data", "rotated_coast_outline.rds"))

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
  dplyr::select(Xr, Yr, pred_cmb, salmon_region = region, month_f)

filter(grid_df, is.na(Xr) | is.na(Yr))

#### 

# Absolute abundance - resident groups (WCVI, CBC, NBC/SEAK)

pabs_residents <- grid_df %>%
  filter(salmon_region %in% c("Northern BC/AK", "Central BC", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       option = "rocket",
                       labels = scales::comma,
                       limits = c(0, max(grid_df$pred_cmb)),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        plot.margin=margin(b = 0, t = 0.2, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray83") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Distribution of juvenile chinook - \"Resident\" stock groupings") +
  NULL

ggsave(pabs_residents,
       filename = here::here("figs", "chinook_abs_abundance_resident.png"), 
       width = 7, height = 6.9)


# Absolute abundance - slow migration groups (Fraser, Salish Sea)

grid_slowmigs <- grid_df %>% filter(grepl( "Fraser|Salish", salmon_region))
  
pabs_slowmigs <- grid_slowmigs %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       option = "rocket",
                       labels = scales::comma,
                       limits = c(0, max(grid_slowmigs$pred_cmb)),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        plot.margin=margin(b = 0, t = 0.2, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray83") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Distribution of juvenile chinook - Salish Sea \"migrating\" stock groupings") +
  NULL

ggsave(pabs_slowmigs,
       filename = here::here("figs", "chinook_abs_abundance_slowmigs_freescale.png"), 
       width = 7, height = 9)


# Absolute abundance - fast migration groups (Columbias, WA/OR Coastal)

grid_fastmigs <- grid_df %>% filter(grepl( "Col.|OR", salmon_region))

pabs_fastmigs <- grid_fastmigs %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       option = "rocket",
                       labels = scales::comma,
                       limits = c(0, max(grid_fastmigs$pred_cmb)),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        plot.margin=margin(b = 0, t = 0.2, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray83") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Distribution of juvenile chinook - Pacific Ocean \"migrating\" stock groupings") +
  NULL

ggsave(pabs_fastmigs,
       filename = here::here("figs", "chinook_abs_abundance_fastmigs_freescale.png"), 
       width = 7, height = 9)




# 
# 
# prot_1 <- grid_df %>%
#   #filter(salmon_region %in% c("WCVI", "Salish Coastal"), month_f %in% c("May", "Jun")) %>%
#   filter(salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Sea", "WA/OR Coastal", "WCVI")) %>%
#   ggplot(data = .) +
#   #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render
#   geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
#   #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
#   scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
#                        labels = scales::comma,
#                        limits = c(0, max(grid_df$pred_cmb)),
#                        trans = ggsidekick::fourth_root_power_trans()) +
#   ggsidekick::theme_sleek() +
#   coord_fixed() +
#   theme(legend.key.height = unit(0.8, "cm"),
#         legend.position="none",
#         plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
#   #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
#   scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
#   scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
#   facet_grid(salmon_region ~ month_f)  +
#   #ggtitle("Juvenile chinoook stock distribution") +
#   NULL

# prot_2 <- grid_df %>%
#   filter(!salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
#   ggplot(data = .) +
#   #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
#   geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
#   #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
#   scale_fill_viridis_c(name = "Predicted\nabundance\nof juvenile\nchinook\nper tow",
#                        labels = scales::comma,
#                        limits = c(0, max(grid_df$pred_cmb)),
#                        trans = ggsidekick::fourth_root_power_trans()) +
#   ggsidekick::theme_sleek() +
#   coord_fixed() +
#   theme(legend.key.height = unit(0.8, "cm"),
#         #legend.position="none",
#         #legend.margin=margin(t = 0, b=0, unit='cm'),
#         plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
#   #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
#   scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
#   scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
#   facet_grid(salmon_region ~ month_f)  +
#   # ggtitle("Juvenile chinoook stock distribution") +
#   NULL
# 
# # This produces a nice output using the patchwork pkg
# ggsave(prot_1 + prot_2,
#        filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_rotated_no_cv_cutoff_11.png"),  
#        width = 11, height = 11)


### Relative abundance plots

grid_rel <- grid_df %>%
  group_by(salmon_region) %>%
  mutate(rel_pred = pred_cmb/max(pred_cmb)) %>%
  dplyr::select(Xr, Yr, rel_pred, salmon_region, month_f)

#########
# col_opts <- c("plasma", "inferno", "magma", "rocket", "mako")

prel_residents <- grid_rel %>%
  filter(salmon_region %in% c("Northern BC/AK", "Central BC", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Relative\nabundance of\njuvenile\nchinook\nwithin\nregional group",
                       option = "viridis",
                       labels = scales::comma,
                       limits = c(0, 1),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        plot.margin=margin(b = 0, t = 0.2, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray83") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Relative distribution of juvenile chinook - \"Resident\" stock groupings") +
  NULL


ggsave(prel_residents,
       filename = here::here("figs", "chinook_relative_abundance_resident.png"), 
       width = 7, height = 6.9)

  
prel_slowmigs <- grid_rel %>%
    filter(grepl( "Fraser|Salish", salmon_region)) %>%
    ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Relative\nabundance of\njuvenile\nchinook\nper tow",
                       option = "viridis",
                       labels = scales::comma,
                       limits = c(0, 1),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
       # legend.position="none",
        plot.margin=margin(b = 0, t = 0.2, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray83") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Relative distribution of juvenile chinook - Salish Sea \"migrating\" stock groupings") +
  NULL
#prel_1

ggsave(prel_slowmigs,
       filename = here::here("figs", "chinook_relative_abundance_slowmigs.png"), 
       width = 7, height = 9)


prel_fastmigs <- grid_rel %>%
  filter(grepl( "Col.|OR", salmon_region)) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Relative\nabundance of\njuvenile\nchinook\nper tow",
                       option = "viridis",
                       labels = scales::comma,
                       limits = c(0, 1),
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
      #  legend.position="none",
        plot.margin=margin(b = 0, t = 0.2, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray83") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Relative distribution of juvenile chinook - Pacific Ocean \"migrating\" stock groupings") +
  NULL
#prel_1

ggsave(prel_fastmigs,
       filename = here::here("figs", "chinook_relative_abundance_fastmigs.png"), 
       width = 7, height = 9)




# prel_1 <- grid_rel %>%
#   filter(salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
#   #filter(salmon_region %in% c("SoG Coastal", "WA/OR Coastal"), month_f %in% c("Jun", "Jul")) %>%
#   ggplot(data = .) +
#   #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render
#   geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
#   #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
#   scale_fill_viridis_c(name = "Relative\nabundance of\njuvenile\nchinook\nper tow",
#                        labels = scales::comma,
#                        limits = c(0, 1),
#                        trans = ggsidekick::fourth_root_power_trans()) +
#   ggsidekick::theme_sleek() +
#   coord_fixed() +
#   theme(legend.key.height = unit(0.8, "cm"),
#         legend.position="none",
#         plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
#   #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
#   scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
#   scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
#   facet_grid(salmon_region ~ month_f)  +
#   #ggtitle("Juvenile chinoook stock distribution")+
#   NULL
# #prel_1
# 
# prel_2 <- grid_rel %>%
#   filter(!salmon_region %in% c("Northern BC/AK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
#   #filter(salmon_region %in% c("Central BC", "WCVI"), month_f %in% c("Jun", "Jul")) %>%
#   ggplot(data = .) +
#   #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render
#   geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
#   #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
#   scale_fill_viridis_c(name = "Predicted\nrelative\nabundance\nof juvenile\nchinook",
#                        labels = scales::comma,
#                        limits = c(0, 1),
#                        trans = ggsidekick::fourth_root_power_trans()) +
#   ggsidekick::theme_sleek() +
#   coord_fixed() +
#   theme(legend.key.height = unit(0.8, "cm"),
#         #legend.position="none",
#         plot.margin=margin(b = 0, t = 0, l = 0, r = 0, unit='cm'),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
#   #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
#   scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
#   scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
#   facet_grid(salmon_region ~ month_f)  +
#   # ggtitle("Juvenile chinoook stock distribution") +
#   NULL
# 
# # This produces a nice output using the patchwork pkg
# ggsave(prel_1 + prel_2,
#        filename = here::here("figs", "chinook_relative_abundance_by_stock_sdmTMB_rotated_no_cv_cutoff_lowerfraser.png"),
#        width = 11, height = 11)


### Fraser stocks

fraser_rel <- grid_rel %>%
  filter(grepl( "Fraser", salmon_region)) %>%
  #filter(salmon_region %in% c("Central BC", "WCVI"), month_f %in% c("Jun", "Jul")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nrelative\nabundance\nof juvenile\nchinook",
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
ggsave(fraser_rel,
       filename = here::here("figs", "chinook_relative_abundance_fraser.png"),
       width = 11, height = 8)

fraser_abs <- grid_df %>%
  filter(grepl( "Fraser", salmon_region)) %>%
  #filter(salmon_region %in% c("Central BC", "WCVI"), month_f %in% c("Jun", "Jul")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nabundance\nof juvenile\nchinook\nper tow",
                       labels = scales::comma,
                       option = "plasma",
                       limits = c(0, max(filter(grid_df, grepl( "Fraser", salmon_region))$pred_cmb)),
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
