###
### Estimating uncertainty of both species-wide and stock-specific predictions 
### across whole grid
###

library(tidyverse)
library(sdmTMB)
#library(sdmTMBextra)
library(here)
library(patchwork)
library(sf)

cutoff <- 5

source("R/cleave_by.R")
source("R/99_rotation_functions.R")

prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))
pred_grid <- readRDS(here("data", "pred_grid_all_coast.rds"))
m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))
m_svc$version
dat_trim <- m_svc$data
table(dat_trim$region)
#dat_trim <- readRDS(here::here("data", "chinook_gsi_counts_fitted.rds"))
m_chinook <- readRDS(here("data", "fits", "fits_list_mdmesh.rds"))$chinook

###  rotating coastline 
sf::sf_use_s2(FALSE)
map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
coast <- sf::st_crop(
  map_data,
  c(xmin = -135, ymin = 46, xmax = -120, ymax = 58.5)
)
plot(st_geometry(coast))
coast_proj <- sf::st_transform(coast, crs = 32609)
# make them in the right format
coast_proj4 <- st_cast(coast_proj, "POLYGON")

rotated_coast <- list()
rotation_angle <- 45
rotation_centre <- c(sum(range(pred_grid$utm_x_1000))/2, sum(range(pred_grid$utm_y_1000))/2)
saveRDS(rotation_centre, here("data", "rot_centre.rds"))

for (i in seq_len(dim(coast_proj4)[1])) {
  rotated_coast[[i]] <- splitrotatepolygon(
    coast_proj4,
    rotation_angle,
    rotation_centre[1] * 1000,
    rotation_centre[2] * 1000
  )
}

rotated_coast <- do.call(rbind, rotated_coast)

#saveRDS(rotated_coast, here::here("data", "rotated_coast_outline.rds"))

# functions for quickly applying and unapplying scale()
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

grid_months_chinook <- prop_grid %>%
  select(-FID, -elevation_meter, -slope_degree, -label, -coast_distance_meter, -depth, - dist_to_coast_km) %>% # removing unused columns to reduce size
  # replicate_df(., "region", levels(dat_trim$region)) %>% # c("WCVI", "North/Central BC", "Columbia")) %>%
  replicate_df(., "month", c(5:12)) %>%
  mutate(month = as.numeric(month), 
         #region = as.factor(region),
         #yday = 160,
         survey_f = as.factor("hss"),
         day_night = as.factor("DAY"),
         scale_depth = scale_est(0, m_chinook$data$target_depth), #-0.5427145, # Change to 0m based on data
         utm_x_1000 = X/1000,
         utm_y_1000 = Y/1000,
         year_adj_f = as.factor(2019),
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 3, month - 3, month + 9), # adjusting to start in April
         scale_month_adj = scale_est(month_adj, m_chinook$data$month_adj),
         NULL) %>%
  select(-X, -Y, -year_adj) 

glimpse(grid_months_chinook)

# 200 simulations to estimate CV for species-wide model
grid_months_list <- grid_months_chinook %>% cleave_by(month_adj)
length(grid_months_list)

cv_list_chinook200 <- map(grid_months_list, function (newdata) {
  gc()
  pred <- predict(m_chinook, newdata = newdata, nsim = 200,
  offset = rep.int(median( m_chinook$data$effort), nrow(newdata)))
  #offset = rep.int(median(dat_trim$effort), nrow(newdata)))
  pred_i <- predict(m_chinook, newdata = newdata, se_fit = FALSE, re_form_iid = NA,
                    # offset = rep.int(median(dat_trim$effort), nrow(newdata)),
                    offset = rep.int(median( m_chinook$data$effort), nrow(newdata)))  %>% 
    pull(est) %>%
    m_chinook$family$linkinv() %>% 
    as.numeric()
  
  month_id <- as.character(unique(newdata$month_adj))
  print(month_id)

  sds <- apply(pred, 1, sd) %>% as.numeric()
  means <- apply(pred, 1, mean) %>% as.numeric()
  
  cvs <- apply(pred, 1, \(x) sd(exp(x)) / mean(exp(x))) %>% as.numeric()
  tibble(pred = pred_i, means = means %>% m_chinook$family$linkinv(), sds = sds, cvs = cvs)
})


all_list_chinook <- map2(grid_months_list, cv_list_chinook200, function (grid, cvs) {
  bind_cols(grid, cvs)
})

cvall_chinook <- bind_rows(all_list_chinook)

cvall_chinook_rot <- bind_cols(cvall_chinook, 
                       gfplot:::rotate_coords(cvall_chinook$utm_x_1000, cvall_chinook$utm_y_1000, 45, 
                                              rotation_centre)) %>%
                                              # c(sum(range(cvall_chinook$utm_x_1000))/2, 
                                              #   sum(range(cvall_chinook$utm_y_1000))/2))) %>%
  mutate(Xr = x * 1000, Yr = y * 1000)

cvall_chinook_rot$cvs %>% range()
# CV cutoff value to remove grid points with cvs higher than cutoff



pred_chinook_cutoff <- cvall_chinook_rot %>%
  mutate(cutoff_keep = if_else(cvs < cutoff, TRUE, FALSE))

# saving cutoffs to join with normalized predictions
saveRDS(pred_chinook_cutoff, here("data", "cvs", "pred_chinook_cutoff.rds"))
#pred_chinook_cutoff <- readRDS(here("data", "cvs", "pred_chinook_cutoff.rds"))

# saving cvs data frame for faster replotting
saveRDS(cvall_chinook_rot, here("data", "cvs", "cvs_chinook.rds"))
#cvall_chinook_rot <- readRDS(here("data", "cvs", "cvs_chinook.rds"))


cvhists <- cvall_chinook_rot %>%
  dplyr::select(Xr, Yr, cvs, month_f) %>%
  ggplot(data = .) +
  geom_histogram(aes(cvs), binwidth = 0.1) +
  facet_wrap(~month_f)  +
  ylab("Count of grid units") +
  xlab("Coefficient of variation") +
  ggtitle("Species-wide model") +
  theme_bw() +
  geom_vline(xintercept = 5, color = "red") +
  NULL

cvhists 

ggsave(cvhists, filename = here("figs", "cv_hist_species.png"), width = 6, height = 5)

chinook_ggdata <- cvall_chinook_rot %>%
  dplyr::select(Xr, Yr, cvs, pred, month_f) %>%
  #filter(cvs < cutoff)
  mutate(cvs2 = if_else(cvs <= cutoff, cvs, NA)) 

cvplot_chinook_rot <- ggplot(data = chinook_ggdata) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = cvs2), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ncoefficient\nof variation",
                       labels = scales::comma,
                       na.value = "orange",
                       # trans = ggsidekick::fourth_root_power_trans(),
                       limits = c(0, cutoff)#max(cvall_chinook_rot$cvs))
  ) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        # for when n_months = 9
        # legend.position = "inside",
        # legend.position.inside = c(0.91, 0.25),
        plot.margin = margin(b = 0, t = 0, r = 0, l = 0, unit = 'cm'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, 
                     limits = range(cvall_chinook_rot$Xr, na.rm = TRUE) + c(-1000,1000), 
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, 
                     limits = range(cvall_chinook_rot$Yr, na.rm = TRUE) + c(-1000,1000), 
                     expand = c(0, 0)) +
  facet_wrap(~month_f, nrow = 2)  
  #ggtitle("Juvenile chinoook stock distribution") +

cvplot_chinook_rot

ggsave(cvplot_chinook_rot, filename = here("figs", "cv_month_species.png"), 
       width = 5, height = 6)


pred_plot_chinook_rot <- ggplot(data = chinook_ggdata) +
  geom_tile(aes(x = Xr, y = Yr, fill = pred), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  scale_fill_viridis_c(name = "Predicted\nnumber of\njuvenile\nchinook\nper tow",
                       labels = scales::comma,
                        trans = ggsidekick::fourth_root_power_trans(),
                       limits = c(0, max(chinook_ggdata$pred))
  ) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position = "inside",
    #    legend.position.inside = c(0.91, 0.25),
        plot.margin = margin(b = 0, t = 0, r = 0, l = 0, unit = 'cm'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(cvall_chinook_rot$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(cvall_chinook_rot$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_wrap(~month_f, nrow = 2)  

ggsave(pred_plot_chinook_rot, filename = here("figs", "pred_month_species.png"), width = 5, height = 6)


### Stock-specific model CVs

grid_months <- prop_grid %>%
  select(-FID, -elevation_meter, -slope_degree, -label, -coast_distance_meter, -depth, - dist_to_coast_km) %>% # removing unused columns to reduce size
  replicate_df(., "region", levels(dat_trim$region)) %>% # c("WCVI", "North/Central BC", "Columbia")) %>%
  replicate_df(., "month", c(5:12)) %>%
  mutate(month = as.numeric(month), 
         region = as.factor(region),
         #yday = 160,
         utm_x_1000 = X/1000,
         utm_y_1000 = Y/1000,
         year_adj_f = as.factor(2019),
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 3, month - 3, month + 9), # adjusting to start in April
         scale_month_adj = scale_est(month_adj, dat_trim$month_adj),
         NULL) %>%
  select(-X, -Y, -year_adj) 

glimpse(grid_months)
table(grid_months$month_f)
table(grid_months$month)
table(grid_months$region)

#cleave_by allows to chop by more than one variable
grid_month_region <- grid_months %>% cleave_by(region, month_adj)
length(grid_month_region)
gc()

# predicting on grid chunks for each region and month
# with 100 simulations
cv_list <- map(grid_month_region, function (newdata) {
  
  gc()
  #browser()
  pred <- predict(m_svc, newdata = newdata, nsim = 100)
                #  offset = rep.int(median(dat_trim$effort), nrow(newdata)))
  pred_i <- predict(m_svc, newdata = newdata, se_fit = FALSE, re_form_iid = NA) %>%
                 # offset = rep.int(median(dat_trim$effort), nrow(newdata))) %>%
    pull(est) %>%
    m_svc$family$linkinv() %>% 
    as.numeric()

  group_id <- as.character(unique(newdata$region))
  month_id <- as.character(unique(newdata$month_adj))
  print(group_id)
  print(month_id)
  
  sds <- apply(pred, 1, sd) %>% as.numeric()
  means <- apply(pred, 1, mean) %>% as.numeric()
  
  
  cvs <- apply(pred, 1, \(x) sd(exp(x)) / mean(exp(x))) %>% as.numeric()
  tibble(pred = pred_i, means = means %>% m_svc$family$linkinv(), sds = sds, cvs = cvs)
})

saveRDS(cv_list, here("data", "cvs", "cv_stock_grid_list.rds"))
#cv_list <- readRDS(here("data", "cvs", "cv_stock_grid_list.rds"))

all_list <- map2(grid_month_region, cv_list, function (grid, cvs) {
  bind_cols(grid, cvs)
})

cvall <- bind_rows(all_list)

cvall_stock_rot <- bind_cols(cvall, 
                          rotate_coords(cvall$utm_x_1000, cvall$utm_y_1000, 45, 
                                                 rotation_centre)) %>%
                                                 # c(sum(range(cvall$utm_x_1000))/2, 
                                                 #   sum(range(cvall$utm_y_1000))/2))) %>%
  mutate(Xr = x * 1000, Yr = y * 1000)


cv_stock_hists <- cvall_stock_rot %>%
  mutate(cvs = Mod(cvs)) %>%
 # filter(cvs < 25) %>%
  dplyr::select(Xr, Yr, cvs, salmon_region = region, month_f) %>%
  ggplot(data = .) +
  geom_histogram(aes(cvs), binwidth = 0.2) +
  facet_grid(salmon_region ~ month_f)  +
  ylab("Count of grid units") +
  xlab("Coefficient of variation") +
  ggtitle("Stock-specific model") +
  theme_bw() +
  geom_vline(xintercept = 5, color = "red") +
  NULL

ggsave(cv_stock_hists, filename = here::here("figs", "stock_cvs_hists.png"), 
       width = 12, height = 14)

pred_stock_cutoff <- cvall_stock_rot %>%
  mutate(cutoff_keep = if_else(cvs < cutoff, TRUE, FALSE))

# saving cutoffs to join with normalized predictions
saveRDS(pred_stock_cutoff, here("data", "cvs","pred_stock_cutoff.rds"))

pred_stock_cutoff <- readRDS(here("data", "cvs","pred_stock_cutoff.rds"))
rotated_coast <- readRDS(here("data", "rotated_coast_outline.rds"))

stock_ggdata <- pred_stock_cutoff %>%
  dplyr::select(Xr, Yr, stock_cvs = cvs, pred, salmon_region = region, month_f,cutoff_keep) %>%
  #filter(cvs < cutoff)
  mutate(stock_cvs2 = if_else(cutoff_keep, stock_cvs, NA)) %>%
  select(-cutoff_keep)



cvplot_stock_rot_1 <- stock_ggdata %>%
  filter(salmon_region %in% c("NBC/SEAK", "Central BC", "Salish Sea", "WA/OR Coastal", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = stock_cvs2), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ncoefficient\nof variation",
                       labels = scales::comma,
                       na.value = "orange",
                       # trans = ggsidekick::fourth_root_power_trans(),
                       limits = c(0, cutoff)
                       ) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        legend.position="none",
        plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(stock_ggdata$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(stock_ggdata$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution") +
  NULL

cvplot_stock_rot_2 <- stock_ggdata %>%
  filter(!salmon_region %in% c("NBC/SEAK", "Central BC", "Salish Sea", "WA/OR Coastal", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = stock_cvs2), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ncoefficient\nof variation",
                       labels = scales::comma,
                       na.value = "orange",
                       # trans = ggsidekick::fourth_root_power_trans(),
                       limits = c(0, cutoff)
  ) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        #legend.position="none",
        plot.margin=margin(b = 0, t = 0, r = 0, l = 0, unit='cm'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(stock_ggdata$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(stock_ggdata$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution") +
  NULL

# rm(all_list)


layout <- "
AB
AB
AB
AB
AB
#B
"


ggsave(cvplot_stock_rot_1 + cvplot_stock_rot_2 + plot_layout(design = layout) +
         plot_annotation(title = 'Stock-specific model',
                         theme = theme(plot.title = element_text(size = 18))),
       filename = here::here("figs", "cv_stock_region_month.png"), 
       width = 11, height = 11)


gc()
