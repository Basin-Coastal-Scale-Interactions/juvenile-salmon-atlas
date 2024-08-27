# generate joint prediction of stock abundance from both models
library(tidyverse)
library(sdmTMB)
library(here)
library(sf)

library(egg)
library(patchwork)

source("R/99_rotation_functions.R")

# loading stock component model predictions and grid
pred_stock <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB.rds")

prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))


ggplot(prop_grid) +
  geom_tile(aes(X,Y))

dat_tbl_dmesh <- readRDS(here::here("data", "fits", "all_dmesh_allcoast.rds"))
chinook_dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))        

# functions for quickly applying and unapplying scale()
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

# grid for predicting 
grid_months <- prop_grid %>%
  replicate_df(., "month", c(3:11)) %>%
  mutate(month = as.numeric(month), 
         #yday = 160,
         utm_x_1000 = X/1000,
         utm_y_1000 = Y/1000,
         year_adj_f = as.factor(2019),
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         scale_month_adj = scale_est(month_adj, chinook_dat$month_adj),
         NULL) %>%
  filter(!month_adj %in% c(11)) #

glimpse(grid_months)

pred_sp <- predict(dat_tbl_dmesh$fit$chinook, newdata = grid_months,
                   se_fit = FALSE, re_form_iid = NA,
                offset = rep.int(median(chinook_dat$effort), nrow(grid_months)))

pred_sp$pred_i <- pred_sp$est %>% dat_tbl_dmesh$fit$chinook$family$linkinv() #pred_linkinv
glimpse(pred_sp)


glimpse(pred_stock)

pred_stock_sub <- pred_stock %>%
  select(utm_x_1000, utm_y_1000, lon, lat, region, month, month_adj, prob_i)

pred_sp_sub <- pred_sp %>%
  select(utm_x_1000, utm_y_1000, month, pred_i)

pred_cmb <- left_join(pred_stock_sub, pred_sp_sub, by = c("utm_x_1000", "utm_y_1000", "month")) %>%
  arrange(region, month, utm_x_1000, utm_y_1000) %>%
  mutate(pred_cmb = pred_i * prob_i)


sum(range(grid_months$utm_x_1000))/2

rotation_angle <- 45
rotation_centre <- c(sum(range(grid_months$utm_x_1000))/2, sum(range(grid_months$utm_y_1000))/2)

pred_cmb_rot <- bind_cols(pred_cmb, 
                          gfplot:::rotate_coords(pred_cmb$utm_x_1000, pred_cmb$utm_y_1000, 45, 
                                                 c(sum(range(grid_months$utm_x_1000))/2, 
                                                   sum(range(grid_months$utm_y_1000))/2)))


pred_cmb_sub <- filter(pred_cmb, region == "Columbia", month == 6)

pred_cmb_sub_rot <- bind_cols(pred_cmb_sub, 
                                 rotate_coords(pred_cmb_sub$utm_x_1000, pred_cmb_sub$utm_y_1000, 
                                               rotation_angle, rotation_centre))

# cell_size <- 1
#
# pred_cmb_sub_sq <- lapply(seq_len(nrow(pred_cmb_sub)), function(i) {
#   #browser()
#   row_dat <- pred_cmb_sub[i, , drop = FALSE]
#   X <- row_dat$utm_x_1000 
#   Y <- row_dat$utm_y_1000 
#   data.frame(
#     X = c(
#       X - cell_size / 2, X + cell_size / 2,
#       X + cell_size / 2, X - cell_size / 2
#     ),
#     Y = c(
#       Y - cell_size / 2, Y - cell_size / 2,
#       Y + cell_size / 2, Y + cell_size / 2
#     ),
#     month = row_dat$month,
#     region = row_dat$region,
#     pred_cmb = row_dat$pred_cmb
#     
#   )
# }) %>% bind_rows()

# pred_cmb_sub_sq_rot <- bind_cols(pred_cmb_sub_sq, 
#                                  rotate_coords(pred_cmb_sub_sq$X, pred_cmb_sub_sq$Y, 
#               rotation_angle, rotation_centre))


source("R/plot_map.R")

# all_coast_km_trans <- all_coast_km %>% st_transform("+proj=omerc +lonc=-129 +lat_0=51.5 +gamma=45 +alpha=0")
# crs <- "+proj=omerc +lonc=-129 +lat_0=51.5 +gamma=45 +alpha=0"

p <- pred_cmb %>%
  filter(region == "WCVI", month == 6) %>% # filtering for faster plotting
  mutate(month_f = month(month, label = TRUE)) %>%
  dplyr::select(utm_x_1000, utm_y_1000, pred_cmb, salmon_region = region, month_f) %>%
  ggplot(data = .) +
  geom_tile( aes(x=utm_x_1000, y=utm_y_1000, fill = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ndistribution\nof juvenile\nchinook by\nstock",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  #coord_fixed() +
  geom_sf(data = all_coast_km, color = "gray80", fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(pred_cmb$utm_x_1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(pred_cmb$utm_y_1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #coord_sf(crs = crs) +
  ggtitle("Juvenile chinoook stock distribution")

p

ggsave(p, filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_sdmTMB.png"), 
       width = 20, height = 17)

ggsave(p, filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_sdmTMB.png"), 
       width = 14, height = 12)

# rotating coastline 
sf::sf_use_s2(FALSE)
map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
coast <- sf::st_crop(
  map_data,
  c(xmin = -135, ymin = 47, xmax = -120, ymax = 57)
)
plot(st_geometry(coast))
coast_proj <- sf::st_transform(coast, crs = 32609)
# make them in the right format
coast_proj4 <-
  st_cast(
    coast_proj,
    "POLYGON"
  )

rotated_coast <- list()

for (i in seq_len(dim(coast_proj4)[1])) {
  rotated_coast[[i]] <- splitrotatepolygon(
    coast_proj4,
    rotation_angle,
    rotation_centre[1]*1000,
    rotation_centre[2]*1000
  )
}

rotated_coast <- do.call(rbind, rotated_coast)
plot(rotated_coast)

##############
# Wrangling data frame one last time for plotting

grid_df <- pred_cmb_rot %>%
  mutate(month_f = month(month, label = TRUE),
         Xr = x * 1000, Yr = y * 1000,
         region = fct_recode(region,
                             "Fraser Summer 4.1" = "Fraser_Summer_4.1",
                             "Northern BC/AK" = "NBC_SEAK"           
                             )) %>%
  mutate(region = fct_relevel(region,
                              c("Northern BC/AK", "Central BC", "SoG Coastal", "WA/OR Coastal", "WCVI",
                                "Col. Coastal", "Col. North", "Fraser Yearling", 
                                "Fraser Summer 4.1", "Upper Col. Yearling")
                              )) %>%
  dplyr::select(Xr, Yr, pred_cmb, salmon_region = region, month_f)
  


saveRDS(rotated_coast, here("data", "rotated_coast.rds"))
saveRDS(grid_df, here("data", "grid_df.rds"))

#prot <- pred_cmb_sub_sq_rot %>%
prot <- grid_df %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  scale_color_viridis_c(name = "Predicted\ndistribution\nof juvenile\nchinook by\nstock",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinoook stock distribution")

prot

ggsave(prot, filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_sdmTMB_rotated_20240726.png"), 
       width = 14, height = 30)

# Dividing plot into 2 for better veiwing

prot_1 <- grid_df %>%
  filter(salmon_region %in% c("Northern BC/AK", "Central BC", "SoG Coastal", "WA/OR Coastal", "WCVI")) %>%
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
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution") +
  NULL

prot_1

prot_2 <- grid_df %>%
  filter(!salmon_region %in% c("Northern BC/AK", "Central BC", "SoG Coastal", "WA/OR Coastal", "WCVI")) %>%
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
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  # ggtitle("Juvenile chinoook stock distribution") +
  NULL

prot_2

library(cowplot)

plot_grid(prot_1, prot_2)

ggsave(plot_grid(prot_1, prot_2), filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_sdmTMB_rotated_20240729.png"), 
       width = 20, height = 16)



# png(here::here("figs", "chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_sdmTMB_rotated_png.png"), 
#           width = 9, height = 10, units = "in", res = 300)
# prot
# dev.off()

### Relative abundance plots

grid_rel <- grid_df %>%
  group_by(salmon_region) %>%
  mutate(rel_pred = pred_cmb/max(pred_cmb)) %>%
  dplyr::select(Xr, Yr, rel_pred, salmon_region, month_f)

prel_1 <- grid_rel %>%
  filter(salmon_region %in% c("Northern BC/AK", "Central BC", "SoG Coastal", "WA/OR Coastal", "WCVI")) %>%
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
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution")+
  NULL
prel_1

prel_2 <- grid_rel %>%
  filter(!salmon_region %in% c("Northern BC/AK", "Central BC", "SoG Coastal", "WA/OR Coastal", "WCVI")) %>%
  #filter(salmon_region %in% c("Central BC", "WCVI"), month_f %in% c("Jun", "Jul")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = rel_pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Relative\nabundance\nof juvenile\nchinook\nper tow",
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
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
 # ggtitle("Juvenile chinoook stock distribution") +
  NULL

prel_2
# 
# ggsave(plot_grid(prel_1, prel_2), filename = here::here("figs", "chinook_relative_abundance_by_stock_rotated_20240729.png"), 
#        width = 20, height = 16)


ggarrange(prel_1, prel_2, ncol = 2, heights = c(60), widths = c(20,20))



plot_grid(prel_1, prel_2, rel_widths = c(1,1.4))


# this produces unevenly sized plots
ggsave(plot_grid(prel_1, prel_2, rel_heights = c(1,1)),
       filename = here::here("figs", "chinook_relative_abundance_by_stock_rotated_20240729-h11.png"), 
      width = 20, height = 16)


prel_1 + prel_2

# This produces a nice output using the patchwork pkg
ggsave(prel_1 + prel_2,
       filename = here::here("figs", "chinook_relative_abundance_by_stock_rotated_20240730-pw.png"), 
       width = 15, height = 11)

ggsave(prot_1 + prot_2,
       filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_rotated_20240730-pw.png"), 
       width = 15, height = 11)

#prot <- pred_cmb_sub_sq_rot %>%
#prot <- 
grid_df %>% 
  filter(salmon_region == "WCVI") %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  scale_color_viridis_c(name = "Predicted\ndistribution\nof juvenile\nchinook by\nstock",
                        labels = scales::comma,#) +#,
                        trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinook stock distribution")

prot

ggsave(prot, filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_sdmTMB_rotated.png"), 
       width = 14, height = 30)







pred_cmb_sub_rot %>%
  #filter(region == "Columbia", month == 6) %>%
  mutate(month_f = month(month, label = TRUE)) %>%
  dplyr::select(x, y, pred_cmb, salmon_region = region, month_f) %>%
  mutate(xr_1000= x *1000, yr_1000 = y *1000) %>% 
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) +
  geom_tile(aes(x = xr_1000, y = yr_1000, fill = pred_cmb), height = 1000, width = 1500) +
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ndistribution\nof juvenile\nchinook by\nstock",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  scale_color_viridis_c(name = "Predicted\ndistribution\nof juvenile\nchinook by\nstock",
                        labels = scales::comma,#) +#,
                        trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm")) +
  geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(pred_cmb_rot$x*1000)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(pred_cmb_sub_rot$y*1000), expand = c(0, 0)) +
  annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  
 # facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinook stock distribution")

prot
