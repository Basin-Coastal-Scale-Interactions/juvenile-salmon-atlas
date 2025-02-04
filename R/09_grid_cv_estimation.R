###
### Estimating uncertainty of both species-wide and stock-specific predictions 
### across whole grid
###

library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)
library(here)
library(patchwork)

rotated_coast <- readRDS(here::here("data", "rotated_coast_outline.rds"))


#-----------------------------------------------------------------------------
#' tidy version of the base `split()` function
#'
#' @param df data.frame to split
#' @param ... unquoted column names from the 'df' data.frame to specify 
#'            the splitting groups
#'
#' @return list of data.frames, one for each group
#-----------------------------------------------------------------------------
cleave_by <- function(df, ...) {
  stopifnot(inherits(df, "data.frame"))
  
  # use tidyeval to get the names ready for dplyr
  grouping <- quos(...)
  
  # Calculate a single number to represent each group
  group_index <- df %>%
    group_by(!!!grouping) %>%
    group_indices()
  
  # do the split by this single group_index variable and return it
  split(df, group_index)
}

prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))
m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB_2025-01-27.rds"))
dat_trim <- readRDS(here::here("data", "chinook_gsi_counts_fitted_2025-01-27.rds"))
m_chinook <- readRDS(here("data", "fits", "fits_list_mdmesh.rds"))$chinook

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
  replicate_df(., "month", c(4:12)) %>%
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

grid_months_list <- grid_months_chinook %>% cleave_by(month_adj)
length(grid_months_list)

cv_list_chinook200 <- map(grid_months_list, function (newdata) {
  gc()
  pred <- predict(m_chinook, newdata = newdata, nsim = 200,
                  offset = rep.int(median(dat_trim$effort), nrow(newdata)))
  pred_i <- predict(m_chinook, newdata = newdata, se_fit = FALSE, re_form_iid = NA,
                    offset = rep.int(median(dat_trim$effort), nrow(newdata))) %>%
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
                                              c(sum(range(cvall_chinook$utm_x_1000))/2, 
                                                sum(range(cvall_chinook$utm_y_1000))/2))) %>%
  mutate(Xr = x * 1000, Yr = y * 1000)

cvall_chinook_rot$cvs %>% range()

cvall_chinook_rot %>%
  ggplot() +
  geom_histogram(aes(cvs)) +
  theme_bw() +
  facet_wrap(~month)

# CV cutoff value to remove grid points with cvs higher than cutoff
cutoff <- 5

pred_chinook_cutoff <- cvall_chinook_rot %>%
  mutate(cutoff_keep = if_else(cvs < cutoff, TRUE, FALSE))

# saving cutoffs to join with normalized predictions
saveRDS(pred_chinook_cutoff, here("data", "pred_chinook_cutoff.rds"))


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
  filter(cvs < cutoff)

cvplot_chinook_rot <- ggplot(data = chinook_ggdata) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = cvs), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ncoefficient\nof variation",
                       labels = scales::comma,
                       # trans = ggsidekick::fourth_root_power_trans(),
                       limits = c(0, cutoff)#max(cvall_chinook_rot$cvs))
  ) +
  ggsidekick::theme_sleek() +
  coord_fixed() +
  theme(legend.key.height = unit(0.8, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.91, 0.25),
        plot.margin = margin(b = 0, t = 0, r = 0, l = 0, unit = 'cm'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  #annotate("text", x = 616000, y = 5290000, color = "red", size = 11, label = "?") +
  scale_x_continuous(name = NULL, limits = range(cvall_chinook_rot$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(cvall_chinook_rot$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_wrap(~month_f, nrow = 2)  
  #ggtitle("Juvenile chinoook stock distribution") +

cvplot_chinook_rot

ggsave(cvplot_chinook_rot, filename = here("figs", "cv_month_species.png"), width = 8, height = 6)

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
        legend.position = "inside",
        legend.position.inside = c(0.91, 0.25),
        plot.margin = margin(b = 0, t = 0, r = 0, l = 0, unit = 'cm'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  geom_sf(data = rotated_coast, color = NA, fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(cvall_chinook_rot$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(cvall_chinook_rot$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_wrap(~month_f, nrow = 2)  

ggsave(pred_plot_chinook_rot, filename = here("figs", "pred_month_species.png"), width = 8, height = 6)


### Stock-specific model CVs

grid_months <- prop_grid %>%
  select(-FID, -elevation_meter, -slope_degree, -label, -coast_distance_meter, -depth, - dist_to_coast_km) %>% # removing unused columns to reduce size
  replicate_df(., "region", levels(dat_trim$region)) %>% # c("WCVI", "North/Central BC", "Columbia")) %>%
  replicate_df(., "month", c(4:12)) %>%
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

grid_month_region <- grid_months %>% cleave_by(region, month_adj)

gc()

# predicting on grid chunks for each region and month
cv_list <- map(grid_month_region, function (newdata) {
  
  gc()
  pred <- predict(m_svc, newdata = newdata, nsim = 100,
                  offset = rep.int(median(dat_trim$effort), nrow(newdata)))
  pred_i <- predict(m_svc, newdata = newdata, se_fit = FALSE, re_form_iid = NA,
                  offset = rep.int(median(dat_trim$effort), nrow(newdata))) %>%
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

saveRDS(cv_list, here("data", "fits", "cv_stock_grid_list.rds"))
cv_list <- readRDS(here("data", "fits", "cv_stock_grid_list.rds"))

all_list <- map2(grid_month_region, cv_list, function (grid, cvs) {
  bind_cols(grid, cvs)
})

cvall <- bind_rows(all_list)

cvall_stock_rot <- bind_cols(cvall, 
                          gfplot:::rotate_coords(cvall$utm_x_1000, cvall$utm_y_1000, 45, 
                                                 c(sum(range(cvall$utm_x_1000))/2, 
                                                   sum(range(cvall$utm_y_1000))/2))) %>%
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
       width = 12, height = 12)

stock_ggdata <- cvall_stock_rot %>%
  dplyr::select(Xr, Yr, stock_cvs = cvs, pred, salmon_region = region, month_f) %>%
  filter(stock_cvs < cutoff)

pred_stock_cutoff <- cvall_stock_rot %>%
  mutate(cutoff_keep = if_else(cvs < cutoff, TRUE, FALSE))


# saving cutoffs to join with normalized predictions
saveRDS(pred_stock_cutoff, here("data", "pred_stock_cutoff.rds"))

cvplot_stock_rot_1 <- stock_ggdata %>%
  filter(salmon_region %in% c("NBC/SEAK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = stock_cvs), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ncoefficient\nof variation",
                       labels = scales::comma,
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
  scale_x_continuous(name = NULL, limits = range(cvall_stock_rot$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(cvall_stock_rot$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution") +
  NULL

cvplot_stock_rot_2 <- stock_ggdata %>%
  filter(!salmon_region %in% c("NBC/SEAK", "Central BC", "Salish Coastal", "WA/OR Coastal", "WCVI")) %>%
  ggplot(data = .) +
  #geom_point(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb), shape = 15) + # this hack works but is slower to render 
  geom_tile(aes(x = Xr, y = Yr, fill = stock_cvs), height = 1000, width = 1500) + # this height/width combo fixes plotting issues
  #geom_polygon(aes(x = xr_1000, y = yr_1000, fill = pred_cmb, color = pred_cmb)) +
  scale_fill_viridis_c(name = "Predicted\ncoefficient\nof variation",
                       labels = scales::comma,
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
  scale_x_continuous(name = NULL, limits = range(cvall_stock_rot$Xr)+c(-1000,1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(cvall_stock_rot$Yr)+c(-1000,1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  #ggtitle("Juvenile chinoook stock distribution") +
  NULL

ggsave(cvplot_stock_rot_1 + cvplot_stock_rot_2 +
         plot_annotation(title = 'Stock-specific model',
                         theme = theme(plot.title = element_text(size = 18))),
       filename = here::here("figs", paste0("cv_stock_region_month_",
                                            ymd(Sys.Date()),".png")), 
       width = 15, height = 11)


