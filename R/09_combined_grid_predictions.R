###
###  Generate joint prediction of stock groupings abundances from both models
### 

library(tidyverse)
library(sdmTMB)
library(here)
library(sf)
library(cowplot)
library(egg)
library(patchwork)
library(sp)

get_ll <- function(x, y, zone, loc){
  points_utm <- SpatialPoints(
    cbind(x, y), 
    proj4string = CRS(paste0("+proj=utm +zone=",zone[1]," +units=m"))
  )
  points_ll <- spTransform(
    points_utm, 
    CRS("+proj=longlat +datum=WGS84")
  )
  
  if (loc == "x") {
    return(coordinates(points_ll)[,1])
  } else if (loc == "y") {
    return(coordinates(points_ll)[,2])
  }
}

source("R/cleave_by.R")
source("R/plot_map.R")
source("R/99_rotation_functions.R")

# loading stock component model predictions and grid
# pred_stock <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB.rds")
pred_stock <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB-no05prop.rds")
table(pred_stock$region)
prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))
pred_grid <- readRDS(here("data", "pred_grid_all_coast.rds"))

# map rotation for better plotting
rotated_coast <- readRDS(here::here("data", "rotated_coast_outline.rds"))

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


centres_m <- rbind(c(sum(range(pred_grid$utm_x_1000))/2, sum(range(pred_grid$utm_y_1000))/2),
c(sum(range(prop_grid$utm_x_1000))/2, sum(range(prop_grid$utm_y_1000))/2))

colnames(centres_m) <- c("x", "y")
centres <- as_tibble(centres_m) %>%
  mutate(source = c("pred_grid", "prop_grid"),
         lon = get_ll(x, y, "9", "x"),
         lat = get_ll(x, y, "9", "y"))       

plot(st_geometry(coast))

ggplot(coast) +
  geom_sf() +
  geom_point(data = centres, aes(lon, lat, color = source))

coast_proj <- sf::st_transform(coast, crs = 32609)
# make them in the right format
coast_proj4 <- st_cast(coast_proj, "POLYGON")

rotated_coast <- list()
rotation_angle <- 45
rotation_centre <- c(sum(range(pred_grid$utm_x_1000))/2, sum(range(pred_grid$utm_y_1000))/2)

saveRDS(rotation_centre, here("data", "rot_centre.rds"))

# function for quickly applying scale(), needed for predictive grid using scaled 
# model parameters
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

# # grid for predicting 
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
  # pred <- predict(m_chinook, newdata = newdata, nsim = 200,
  #                 offset = rep.int(median( m_chinook$data$effort), nrow(newdata)))
  #offset = rep.int(median(dat_trim$effort), nrow(newdata)))
  pred_i <- predict(m_chinook, newdata = newdata, se_fit = FALSE, re_form_iid = NA,
                    # offset = rep.int(median(dat_trim$effort), nrow(newdata)),
                    offset = rep.int(median(m_chinook$data$effort), nrow(newdata)))  %>% 
    pull(est) %>%
    m_chinook$family$linkinv() %>% 
    as.numeric()
  
  month_id <- as.character(unique(newdata$month_adj))
  print(month_id)
  
  tibble(pred = pred_i)
})

# pred_sp <- predict(dat_tbl_dmesh$fit$chinook, newdata = grid_months,
#                    se_fit = FALSE, re_form_iid = NA,
#                 offset = rep.int(median(chinook_dat$effort), nrow(grid_months)))
# 
# pred_sp$pred_i <- pred_sp$est %>% dat_tbl_dmesh$fit$chinook$family$linkinv() #pred_linkinv


all_list_chinook <- map2(grid_months_list, cv_list_chinook200, function (grid, cvs) {
  bind_cols(grid, cvs)
})

cvall_chinook <- bind_rows(all_list_chinook)

cvall_chinook_rot <- bind_cols(cvall_chinook, 
                               gfplot:::rotate_coords(cvall_chinook$utm_x_1000, 
                                                      cvall_chinook$utm_y_1000, 
                                                      45, 
                                                      rotation_centre)) %>%
  # c(sum(range(cvall_chinook$utm_x_1000))/2, 
  #   sum(range(cvall_chinook$utm_y_1000))/2))) %>%
#  mutate(Xr = x * 1000, Yr = y * 1000)

pred_stock_sub <-  pred_stock %>%
  select(utm_x_1000, utm_y_1000, lon, lat, region, month, month_adj, prob_i)

filter(pred_stock_sub, is.na(Xr) | is.na(Yr))

pred_sp_sub <- cvall_chinook_rot %>%
  select(utm_x_1000, utm_y_1000, month, pred, Xr, Yr)

filter(pred_sp_sub, is.na(Xr) | is.na(Yr))

pred_stock_sub

pred_cmb <- left_join(pred_stock_sub, pred_sp_sub, 
                      by = c("utm_x_1000", "utm_y_1000", "month")) %>%
  arrange(region, month, utm_x_1000, utm_y_1000) %>%
  mutate(pred_cmb = pred * prob_i)

saveRDS(pred_cmb, file = here("data", "pred_cmb.rds"))
