###
###  Generate joint prediction of stock abundance from both models
### 

library(tidyverse)
library(sdmTMB)
library(here)
library(sf)
library(cowplot)
library(egg)
library(patchwork)

source("R/plot_map.R")
source("R/99_rotation_functions.R")

# loading stock component model predictions and grid
# pred_stock <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB.rds")
pred_stock <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB-no05prop.rds")
table(pred_stock$region)
prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))

# map rotation for better plotting
rotated_coast <- readRDS(here::here("data", "rotated_coast_outline.rds"))

# cv cutoffs for both models
pred_chinook_cutoff <- readRDS(here("data", "pred_chinook_cutoff.rds"))
pred_stock_cutoff <- readRDS(here("data", "pred_stock_cutoff.rds"))

table(pred_stock_cutoff$region)

# estimated predicted abundance values across grid
pred_chinook_cutoff$pred %>% range() %>% format(scientific = FALSE)
filter(pred_chinook_cutoff, pred == max(pred_chinook_cutoff$pred))

table(pred_stock_cutoff$month_f)

# ggplot(prop_grid) +# gmonth_fgplot(prop_grid) +
#   geom_tile(aes(X,Y))

dat_tbl_dmesh <- readRDS(here::here("data", "fits", "all_dmesh_allcoast.rds"))
# chinook_dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))        
chinook_dat <- dat_tbl_dmesh[dat_tbl_dmesh$species == "chinook", "data"][[1]][[1]]

# highest chinook catch
chinook_dat$n_juv %>% max()


# functions for quickly applying and unapplying scale()
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

# # grid for predicting 
# grid_months <- prop_grid %>%
#   select(-FID, -elevation_meter, -dist_to_coast_km, -slope_degree, -coast_distance_meter, -depth, -label) %>%
#   replicate_df(., "month", c(5:12)) %>%
#   mutate(month = as.numeric(month), 
#          #yday = 160,
#          utm_x_1000 = X/1000,
#          utm_y_1000 = Y/1000,
#          survey_f = as.factor("hss"),
#          day_night = as.factor("DAY"),
#          scale_depth = -0.5427145, # Change to 0m based on data
#          year_adj_f = as.factor(2019),
#          month_f = month(month, label = TRUE),
#          month_adj = ifelse(month > 3, month - 3, month + 9), # adjusting to start in April
#          scale_month_adj = scale_est(month_adj, chinook_dat$month_adj),
#          NULL) #%>%
#   #filter(!month_adj %in% c(11)) #
# 
# glimpse(grid_months)
# 
# pred_sp <- predict(dat_tbl_dmesh$fit$chinook, newdata = grid_months,
#                    se_fit = FALSE, re_form_iid = NA,
#                 offset = rep.int(median(chinook_dat$effort), nrow(grid_months)))
# 
# pred_sp$pred_i <- pred_sp$est %>% dat_tbl_dmesh$fit$chinook$family$linkinv() #pred_linkinv
# 
# glimpse(pred_sp)
# glimpse(pred_stock)

# adding cv cutoffs to species-wide model
pred_chinook <- rename(pred_chinook_cutoff, sp_cutoff_keep = cutoff_keep)

# pred_chinook <- left_join(pred_sp, select(pred_chinook_cutoff, 
#                                           lon, lat, month_adj, Xr, Xr, pred, 
#                                           cvs, sp_cutoff_keep = cutoff_keep),
#                           by = c("lon", "lat", "month_adj"))

#mutate(pred_chinook, diff = pred_i - pred) %>% pull(diff) %>% range()

# plot(pred_i ~ pred, pred_chinook)
# ggplot(pred_chinook) +
#   geom_point(aes(pred, pred_i))
#lm(pred_i ~ pred, pred_chinook) %>% summary()

# this one requires pred_stock_cutoff, which is computationally intensive as it
# calculates the cvs for the stock composition model
# pred_stock_sub <- pred_stock %>%
#   select(utm_x_1000, utm_y_1000, lon, lat, region, month, month_adj, prob_i) %>%
#   left_join(select(pred_stock_cutoff, lon, lat, region, month_adj, Xr, Yr, cvs, 
#                    stock_cutoff_keep = cutoff_keep),
#             by = c("lon", "lat", "month_adj", "region"))

XYrot <- filter(pred_stock_cutoff, month_adj == 4, region == "WCVI") %>%
  select(lon, lat,  Xr, Yr)
pred_stock_sub <-  pred_stock %>%
    select(utm_x_1000, utm_y_1000, lon, lat, region, month, month_adj, prob_i) %>%
     left_join(XYrot, by = c("lon", "lat"))

filter(pred_stock_sub, is.na(Xr) | is.na(Yr))

pred_sp_sub <- pred_chinook %>%
  select(utm_x_1000, utm_y_1000, month, pred, sp_cutoff_keep)

pred_cmb <- left_join(pred_stock_sub, pred_sp_sub, by = c("utm_x_1000", "utm_y_1000", "month")) %>%
  arrange(region, month, utm_x_1000, utm_y_1000) %>%
  mutate(pred_cmb = pred * prob_i)

filter(pred_cmb, is.na(Xr) | is.na(Yr))

options(pillar.sigfig = 5)
pred_cmb %>%
  group_by(region) %>%
  filter(stock_cutoff_keep) %>%
  filter(pred_cmb == max(pred_cmb)) %>%
  select(lat, lon, region, month, sp_cutoff_keep, stock_cutoff_keep, 
         cvs, pred, prob_i, pred_cmb)


# pred_cmb_rot <- bind_cols(pred_cmb, 
#                           gfplot:::rotate_coords(pred_cmb$utm_x_1000, pred_cmb$utm_y_1000, 45, 
#                                                  c(sum(range(grid_months$utm_x_1000))/2, 
#                                                    sum(range(grid_months$utm_y_1000))/2)))
# 
# 
# pred_cmb_sub <- filter(pred_cmb, region == "Columbia", month == 6)
# 
# pred_cmb_sub_rot <- bind_cols(pred_cmb_sub, 
#                                  rotate_coords(pred_cmb_sub$utm_x_1000, pred_cmb_sub$utm_y_1000, 
#                                                rotation_angle, rotation_centre))

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
ggsave(p, filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_WCVI_Jun.png"), 
       width = 14, height = 12)


##############
# Wrangling data frame one last time for plotting

saveRDS(pred_cmb, file = here("data", "pred_cmb_fraserfall.rds"))

pred_cmb_cut <- filter(pred_cmb, stock_cutoff_keep == TRUE & sp_cutoff_keep == TRUE)

grid_df <- pred_cmb_cut %>%
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

saveRDS(grid_df, here("data", "grid_df.rds"))

#### 
grid_df <- readRDS(here("data", "grid_df.rds"))

# rm(pred_cmb)
# rm(pred_cmb_cut)
# rm(pred_stock_cutoff)
# rm(pred_stock)

#prot <- pred_cmb_sub_sq_rot %>%
# prot <- grid_df %>%
#   ggplot(data = .) +
#   #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
#   geom_tile(aes(x = Xr, y = Yr, fill = pred_cmb), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
#   #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
#   scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
#                        labels = scales::comma,#) +#,
#                        trans = ggsidekick::fourth_root_power_trans()) +
#   scale_color_viridis_c(name = "Predicted\ndistribution\nof juvenile\nchinook by\nstock",
#                        labels = scales::comma,#) +#,
#                        trans = ggsidekick::fourth_root_power_trans()) +
#   ggsidekick::theme_sleek() +
#   coord_fixed() +
#   theme(legend.key.height = unit(0.8, "cm"),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   geom_sf(data = rotated_coast, color = "gray80", fill = "gray90") +
#   annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
#   scale_x_continuous(name = NULL, limits = range(grid_df$Xr)+c(-1000,1000), expand = c(0, 0)) +
#   scale_y_continuous(name = NULL, limits = range(grid_df$Yr)+c(-1000,1000), expand = c(0, 0)) +
#   facet_grid(salmon_region ~ month_f)  +
#   ggtitle("Juvenile chinoook stock distribution")
# 
# prot
# 
# ggsave(prot, 
#        filename = here::here("figs", 
#                              paste0("chinook_abundance_by_stock_sdmTMB_cmb_svc_mapped_tauZ_rotated_no05_",
#                                     ymd(Sys.Date()),".png")), 
#        width = 14, height = 30)

# Dividing plot into 2 for better viewing

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

gc()

# This produces a nice output using the patchwork pkg
ggsave(prot_1 + prot_2,
       filename = here::here("figs", "chinook_abundance_by_stock_sdmTMB_rotated_pw_no05.png"),  
       width = 11, height = 11)


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

#prel_2
# 
# ggsave(plot_grid(prel_1, prel_2), filename = here::here("figs", "chinook_relative_abundance_by_stock_rotated_20240729.png"), 
#        width = 20, height = 16)


# ggarrange(prel_1, prel_2, ncol = 2, heights = c(60), widths = c(20,20))
# 
# 
# 
# plot_grid(prel_1, prel_2, rel_widths = c(1,1.4))
# 
# 
# # this produces unevenly sized plots
# ggsave(plot_grid(prel_1, prel_2, rel_heights = c(1,1)),
#        filename = here::here("figs", "chinook_relative_abundance_by_stock_rotated_20240729-h11.png"), 
#       width = 20, height = 16)
# 
# prel_1 + prel_2

# This produces a nice output using the patchwork pkg
ggsave(prel_1 + prel_2,
       filename = here("figs", "chinook_relative_abundance_by_stock_rotated_pw_no05.png"), 
       width = 11, height = 11)

