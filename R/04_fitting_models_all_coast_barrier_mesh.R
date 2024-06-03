library(tidyverse)
library(ggplot2)
library(sdmTMB)
library(sdmTMBextra)
library(here)
library(purrr)
library(ggsidekick)
library(inlabru)
library(sf)

source(here("R", "plot_map.R"))

# Need coast polygon in km for some graphs
all_coast_km <- rbind(rnaturalearth::ne_states( "United States of America", 
                                                returnclass = "sf"), 
                      rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -135, ymin = 46.25, xmax = -122.25, ymax = 55.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=km"))

dat_in <- readRDS(here::here("data", "catch_survey_sbc_withnorth.rds")) 
dat_tbl_mout3 <- readRDS(here("data", "fits", "all_mout3_allcoast_withpred.rds"))
dat_tbl_dmesh <- readRDS(here("data", "fits", "all_dmesh_allcoast_withpred.rds"))

spde <- readRDS(here::here("data","spde_mesh_allcoast.rds"))        
spde2 <- readRDS(here::here("data","spde2_mesh_allcoast.rds"))        

dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    month_f = month(date, label = TRUE),
    month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night)) %>% 
  filter(!species == "steelhead") %>% 
  droplevels()

dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest()

dat_tbl

# Adding barriers to mesh using coastline
bmesh01 <- sdmTMBextra::add_barrier_mesh(
  spde2,
  all_coast, # from plot_mat.R
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot=TRUE
)
bmesh01$mesh$n
spde2$mesh$n # same number of knots

# Second and third meshes with higher range fractions
bmesh02 <- sdmTMBextra::add_barrier_mesh(
  spde2,
  all_coast, # from plot_mat.R
  range_fraction = 0.2,
  proj_scaling = 1000,
  plot=TRUE
)
bmesh02$mesh$n

bmesh03 <- sdmTMBextra::add_barrier_mesh(
  spde2,
  all_coast, # from plot_mat.R
  range_fraction = 0.4,
  proj_scaling = 1000,
  plot=TRUE
)
bmesh03$mesh$n

mesh_df_water <- bmesh01$mesh_sf[bmesh01$normal_triangles, ]
mesh_df_land <- bmesh01$mesh_sf[bmesh01$barrier_triangles, ]


g <- ggplot() +
  inlabru::gg(bmesh01$mesh) +
  coord_fixed(clip = "off") +
  geom_point(aes(utm_x_1000, utm_y_1000),
             shape = "x",
             size = 0.75,
             data = dat,
             inherit.aes = FALSE) +
  geom_sf(data = all_coast_km, color = "gray80", fill = "darkgreen", alpha = 0.5) +
  geom_sf(data = mesh_df_land, size = 1, colour = "darkgreen")+
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  NULL
g

g$layers <- g$layers[c(4,5,6,3,1,2)]
g


png(here("figs", "barrier_mesh.png"), 
    height = 10, width = 10, units = "in", res = 300)
g
dev.off() 



#####################

# pararellizing species models 
fits_list_barrier01 <- map( #future_map(
  dat_tbl$data,
  function(x) {
    sdmTMB(n_juv ~ 0 + s(month_adj, bs = "tp", k = 6) + # tp smoother is default 
             day_night + survey_f + 
             # scale_dist + 
             scale_depth + (1 | year_f),
           offset = "effort",
           # knots = list(month = c(0, 12)),
           data = x,
           mesh = bmesh01, # barrier mesh
           family = sdmTMB::nbinom2(),
           spatial = "on", # No need to estimate "average"/spatial only RF, 
           #tends to collapse with RW/also when first time step is flat
           share_range = FALSE,
           time = "month_adj",
           extra_time = c(2,11), 
           spatiotemporal = "RW", # RW = AR1 with cor of 1, Philina said only 
           # viable option with temporally patchy data (Greenland Halibut analysis)
           anisotropy = FALSE, # No anisotropy when barrier mesh is used
           silent = FALSE,
           priors = sdmTMBpriors(
             matern_ = pc_matern(
               range_gt = 10, # distance in km 
               sigma_lt = 3
             ), matern_st = pc_matern(
               range_gt = 10, # distance in km 
               sigma_lt = 3
             )
           )
    )
  }
)

# pararellizing species models 
fits_list_barrier04 <- map( #future_map(
  dat_tbl$data,
  function(x) {
    sdmTMB(n_juv ~ 0 + s(month_adj, bs = "tp", k = 6) + # tp smoother is default 
             day_night + survey_f + 
             # scale_dist + 
             scale_depth + (1 | year_f),
           offset = "effort",
           # knots = list(month = c(0, 12)),
           data = x,
           mesh = spde2, # barrier mesh
           family = sdmTMB::nbinom2(),
           spatial = "off", # No need to estimate "average"/spatial only RF, 
           #tends to collapse with RW/also when first time step is flat
           #share_range = FALSE,
           time = "month_adj",
           extra_time = c(2,11), 
           spatiotemporal = "RW", # RW = AR1 with cor of 1, Philina said only 
           # viable option with temporally patchy data (Greenland Halibut analysis)
           anisotropy = TRUE, # No anisotropy when barrier mesh is used
           silent = FALSE #,
           # priors = sdmTMBpriors(
           #   matern_st = pc_matern(
           #     range_gt = 10, # distance in km 
           #     sigma_lt = 3
           #   ),
           #)
    )
  }
)

sanity_check_barrier01 <- map(fits_list_barrier01, 
                    function (x) unlist(sanity(x)) %>% 
                                  all()) %>% 
                                  unlist

sanity_check_barrier04 <- map(fits_list_barrier04, 
                            function (x) unlist(sanity(x)) %>% 
                              all()) %>% 
                                unlist

fits_list_barrier04[[1]] %>% sanity

dat_tbl_barrier <- dat_tbl %>%
  mutate(fit = fits_list_barrier,
         pass_sanity = sanity_check_barrier)

### Fitting models with second with range_fraction = 0.2
# pararellizing species models 
fits_list_barrier2 <- map( #future_map(
  dat_tbl$data,
  function(x) {
    sdmTMB(n_juv ~ 0 + s(month_adj, bs = "tp", k = 6) + # tp smoother is default 
             day_night + survey_f + 
             # scale_dist + 
             scale_depth + (1 | year_f),
           offset = "effort",
           # knots = list(month = c(0, 12)),
           data = x,
           mesh = spde2b2, # barrier mesh range_fraction = 0.2
           family = sdmTMB::nbinom2(),
           spatial = "off", # No need to estimate "average"/spatial only RF, 
           #tends to collapse with RW/also when first time step is flat
           time = "month_adj",
           extra_time = c(2,11), 
           spatiotemporal = "RW", # RW = AR1 with cor of 1, Philina said only 
           # viable option with temporally patchy data (Greenland Halibut analysis)
           anisotropy = FALSE, # No anisotropy when barrier mesh is used
           silent = FALSE
    )
  }
)

sanity_check_barrier2 <- map(fits_list_barrier2, 
                            function (x) unlist(sanity(x)) %>% 
                              all()) %>% unlist

dat_tbl_barrier <- dat_tbl_barrier %>%
  mutate(fit2 = fits_list_barrier2,
         pass_sanity2 = sanity_check_barrier2)

### Plotting distribution maps 


index_grid <- readRDS(here("data", "pred_grid_all_coast.rds")) %>%
  mutate(year = 2019L, 
         year_f = as.factor(year),
         survey_f = as.factor("hss"),
         target_depth = 0,
         day_night = "DAY",
         scale_depth = -0.5427145)

grid_months <- index_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)

dim(grid_months)

# Monthly predictions for all species
all_pred_months_barrier <- purrr::map(
  dat_tbl_barrier$fit,
  ~ predict(.x, newdata = grid_months, 
            offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 
)

all_pred_months_barrier2 <- purrr::map(
  dat_tbl_barrier$fit2,
  ~ predict(.x, newdata = grid_months, 
            offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 
)

# saving to df ############# and adding adjusted month factor levels so maps start from May (except chum)
dat_tbl_barrier$pred_grid <- all_pred_months_barrier
dat_tbl_barrier$pred_grid2 <- all_pred_months_barrier2


saveRDS(dat_tbl_barrier, here::here("data", "fits", "dat_tbl_barrier.rds"))

### Looking at smoother fits

dat_tbl_barrier <- readRDS(here::here("data", "fits", "dat_tbl_barrier.rds"))



plot_smoother <- function (fit) {

  # create data frame with time steps being smoothed to predict with 
  nd <- data.frame(
    month = 1:12L,
    survey_f = "hss",
    day_night = "DAY",
    scale_depth = 0,
    scale_dist = 0,
    year_f = "2010"
  ) %>%
    mutate(month_adj = ifelse(month > 2, # adjusting to start in March
                              month - 2, 
                              month + 10))
  
  # predicting using NA for both re_form and re_form_iid
  p <- predict(fit, 
               newdata = nd, se_fit = TRUE, offset = rep(-6.287114, nrow(nd)),
               re_form = NA, re_form_iid = NA)
  
  # plotting...
  ggplot(p, aes(month_adj, exp(est),
                ymin = exp(est - 1.96 * est_se),
                ymax = exp(est + 1.96 * est_se))) +
    geom_line() +
    geom_ribbon(alpha = 0.4) +
    scale_x_continuous() +
    coord_cartesian(expand = F) +
    labs(x = "month_adj", y = "Count")
}



smb01 <- plot_smoother(dat_tbl_barrier$fit[[1]]) + 
  ggtitle("Smoother of month_adj in chinook model with barrier = 0.1")

smb02 <- plot_smoother(dat_tbl_barrier$fit2[[1]]) + 
  ggtitle("Smoother of month_adj in chinook model with barrier = 0.2")

sman <- plot_smoother(dat_tbl_mout3$fit[[1]]) + 
  ggtitle("Smoother of month_adj in chinook model without barrier (anisotropy)")

png(here::here("figs", paste0("smoother_comp_chinook_anisotropy_barrier.png")),
    height = 10, width = 7, units = "in", res = 300)
cowplot::plot_grid(sman, smb02, smb01, nrow = 3)
dev.off() 

nd <- data.frame(
  #week = seq(4:48),
  #week = seq(1:52),
  month = 1:12L, #c(1,1,rep(1:12, each = 4), 12,12),
  #year_c = 0, # a chosen year
  survey_f = "hss",
  day_night = "DAY",
  scale_depth = 0,
  scale_dist = 0,
  year_f = "2010"
) %>%
  mutate(month_adj = ifelse(month > 2,  # adjusting to start in March
                            month - 2, 
                            month + 10))

p <- predict(dat_tbl_barrier$fit[[4]], 
             newdata = nd, se_fit = TRUE, offset = rep(-6.287114, nrow(nd)),
             re_form = NA, re_form_iid = NA)

ggplot(p, aes(month_adj, exp(est),
              ymin = exp(est - 1.96 * est_se),
              ymax = exp(est + 1.96 * est_se))) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_x_continuous() +
  coord_cartesian(expand = F) +
  labs(x = "month_adj", y = "Count")


### looping through to obtain distribution maps for all species

for (i in 1:nrow(dat_tbl_barrier)) {
  species <- dat_tbl_barrier$species[i]
  print(species)
  
  p <- plot_map(dat_tbl_barrier$pred_grid[[i]], exp(est), area = "BC",
                show_raw_data = FALSE, dat_tbl_barrier$data[[i]]) + 
    annotation_scale(location = "bl", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month (barrier mesh (0.1)")) +
    #scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
    #scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  png(here::here("figs", paste0("Barrier_mesh_atlas_noraw_", species, "_juv_months_barrier01.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}

### looping through to obtain distribution maps for all species

for (i in 1:nrow(dat_tbl_barrier)) {
  species <- dat_tbl_barrier$species[i]
  print(species)
  
  p <- plot_map(dat_tbl_barrier$pred_grid2[[i]], exp(est), area = "BC",
                show_raw_data = FALSE, dat_tbl_barrier$data[[i]]) + 
    annotation_scale(location = "bl", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month (barrier mesh 0.2)")) +
    #scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
    #scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  png(here::here("figs", paste0("Barrier_mesh_atlas_noraw_", species, "_juv_months_barrier02.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}


### DETAIL OF INLETS AND SOUNDS

# Barkley Sound
barkley_grid <- filter(index_grid,
                       lon <= -124.5 & lon >= -126 & 
                         lat >= 48.6 & lat <= 49.5)


# northern VI inlets
northvi_grid <- filter(index_grid,
                       lon <= -126 & lon >= -129 & 
                         lat >= 49 & lat <= 50.7)


# northern VI inlets
qcss_grid <- filter(index_grid, 
                    lon <= -125.5 & lon >= -132 & 
                      lat >= 50 & lat <= 53)



barkley_grid_months <- barkley_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)    

northvi_grid_months <- northvi_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)    


qcss_grid_months <- qcss_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)    



barkley_pred_months <- purrr::map(
  dat_tbl_barrier$fit,
  ~ predict(.x, newdata = barkley_grid_months, 
            offset = rep.int(-6.287114, nrow(barkley_grid_months)), re_form_iid = NA) 
)

northvi_pred_months <- purrr::map(
  dat_tbl_barrier$fit,
  ~ predict(.x, newdata = northvi_grid_months, 
            offset = rep.int(-6.287114, nrow(northvi_grid_months)), re_form_iid = NA) 
)

qcss_pred_months <- purrr::map(
  dat_tbl_barrier$fit,
  ~ predict(.x, newdata = qcss_grid_months, 
            offset = rep.int(-6.287114, nrow(qcss_grid_months)), re_form_iid = NA) 
)


### Creating inlet plots

# looking at min and max lat/lon to see equivalent UTM coordinates
barkley_grid %>%
  slice(c(which.min(lon), which.max(lon)))
barkley_grid %>%
  slice(c(which.min(lat), which.max(lat)))


### Barkley
for (i in 1:nrow(dat_tbl_barrier)) {
  species <- dat_tbl_barrier$species[i]
  print(species)
  
  p <- plot_map(barkley_pred_months[[i]], exp(est), area = "Barkley",
                show_raw_data = TRUE, dat_tbl_barrier$data[[i]]) + 
    annotation_scale(location = "tr", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month around Barkley Sound")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = c(718406, 812406), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5388664, 5477664), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
  png(here::here("figs", paste0("Barrier_mesh_Barkley_atlas_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}


northvi_grid %>%
  slice(c(which.min(lon), which.max(lon)))
northvi_grid %>%
  slice(c(which.min(lat), which.max(lat))) 

### North VI
for (i in 1:nrow(dat_tbl_barrier)) {
  species <- dat_tbl_barrier$species[i]
  print(species)
  
  p <- plot_map(northvi_pred_months[[i]], exp(est), area = "BC",
                show_raw_data = TRUE, dat_tbl_barrier$data[[i]]) + 
    annotation_scale(location = "tr", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month around northern WCVI")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = c(500406, 718406), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5429664, 5618664), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
  png(here::here("figs", paste0("Barrier_mesh_northvi_atlas_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}

# Queen Charlotte Strait & Sound

qcss_grid %>%
  slice(c(which.min(lon), which.max(lon)))
qcss_grid %>%
  slice(c(which.min(lat), which.max(lat))) 

for (i in 1:nrow(dat_tbl_barrier)) {
  species <- dat_tbl_barrier$species[i]
  print(species)
  
  p <- plot_map(qcss_pred_months[[i]], exp(est), area = "BC",
                show_raw_data = TRUE, dat_tbl_barrier$data[[i]]) + 
    annotation_scale(location = "bl", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month around QC Strait and Sound")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = c(292406, 724406), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5538664, 5872664), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
  png(here::here("figs", paste0("Barrier_mesh_qcss_atlas_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}


# Finding and plotting grid cell with max chinook predicted abundance
maxchinook <- dat_tbl_barrier$pred_grid[[1]] %>% 
  mutate(expest = exp(est)) %>% 
  filter(expest == max(expest)) %>% 
  select(month_f, X, Y, utm_x_1000, utm_y_1000, expest)


maxchinook2 <- dat_tbl_barrier$pred_grid2[[1]] %>% 
  mutate(expest = exp(est)) %>% 
  filter(expest == max(expest)) %>% 
  select(month_f, X, Y, utm_x_1000, utm_y_1000, expest)

g + geom_point(aes(utm_x_1000, utm_y_1000), 
               data = maxchinook, size =2, color = "red")

# function to quickly provide plus and minus values around given value
# useful for plotting around specific UTM coordinates
# pm <- function(x, pm) {
#   return(c(x - pm, x + pm))
# }

pm(697406, 10000)

# zooming into chinook max predicted abundance
pcm <- plot_map(dat_tbl_barrier$pred_grid[[1]], exp(est), area = "BC",
                show_raw_data = TRUE, dat_tbl_barrier$data[[1]],
               show_coast = "coastline") + 
    annotation_scale(location = "bl", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile chinook around area of max pred abundance (barrier 0.1)")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = pm(maxchinook$X, 30000), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = pm(maxchinook$Y, 30000), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
png(here::here("figs", paste0("Barrier_mesh_zoom_maxpred_atlas_chinook_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(pcm) # need to print() plot inside for loop!
dev.off() 

# zooming into chinook max predicted abundance
pcm2 <- plot_map(dat_tbl_barrier$pred_grid2[[1]], exp(est), area = "BC",
                show_raw_data = TRUE, dat_tbl_barrier$data[[1]],
                show_coast = "coastline") + 
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  ggtitle(paste0("Predicted distribution of juvenile chinook around area of max pred abundance (barrier 0.2)")) +
  # added new plot ranges based on UTM
  scale_x_continuous(name = NULL, limits = pm(maxchinook$X, 30000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = pm(maxchinook$Y, 30000), expand = c(0, 0)) +
  theme(legend.key.height = unit(0.06, 'npc'))

png(here::here("figs", paste0("Barrier_mesh2_zoom_maxpred_atlas_chinook_juv_months.png")),
    height = 7, width = 10, units = "in", res = 300)
print(pcm2) # need to print() plot inside for loop!
dev.off() 

# zooming into chinook max predicted abundance
pcan <- plot_map(dat_tbl_mout3$pred_grid[[1]], exp(est), area = "BC",
                 show_raw_data = TRUE, dat_tbl_barrier$data[[1]],
                 show_coast = "coastline") + 
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  ggtitle(paste0("Predicted distribution of juvenile chinook around area of max pred abundance (no barrier)")) +
  # added new plot ranges based on UTM
  scale_x_continuous(name = NULL, limits = pm(maxchinook$X, 30000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = pm(maxchinook$Y, 30000), expand = c(0, 0)) +
  theme(legend.key.height = unit(0.06, 'npc'))

png(here::here("figs", paste0("No_barrier_zoom_maxpred_atlas_chinook_juv_months.png")),
    height = 7, width = 10, units = "in", res = 300)
pcan # need to print() plot inside for loop!
dev.off() 



png(here::here("figs", paste0("barrier_mesh_around_max_chinook_pred.png")),
    height = 7, width = 10, units = "in", res = 300)
g + scale_x_continuous(name = NULL, limits = pm(maxchinook$X, 90000)/1000, expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = pm(maxchinook$Y, 90000)/1000, expand = c(0, 0)) +
  geom_point(aes(utm_x_1000, utm_y_1000),data = maxchinook, size =3, color = "red") +
  geom_sf(data = mesh_df_land, size = 2, colour = "darkgreen") +
  geom_sf(data = mesh_df_water, size = 2, colour = "blue")
dev.off()


filter(index_grid, lat == 49.0)


index_grid[which.min(abs(index_grid$lat - 48.7)), "Y"]
index_grid[which.min(abs(index_grid$lon - -125)), "X"]

table(index_grid$lat)


# zooming into Swiftsure Bank
pswi <- plot_map(dat_tbl_mout3$pred_grid[[1]], exp(est), area = "BC",
                 show_raw_data = TRUE, dat_tbl_barrier$data[[1]],
                 show_coast = TRUE) + 
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  ggtitle(paste0("Predicted distribution of juvenile chinook around area of max pred abundance (no barrier)")) +
  scale_x_continuous(name = NULL, limits = pm(794406, 60000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = pm(5396664, 60000), expand = c(0, 0)) +
  theme(legend.key.height = unit(0.06, 'npc'))

png(here::here("figs", paste0("No_barrier_zoom_maxpred_atlas_chinook_juv_months.png")),
    height = 7, width = 10, units = "in", res = 300)
pswi
dev.off() 




dat_tbl_mout3$fit[[1]]
dat_tbl_barrier$fit[[1]]
dat_tbl_barrier$fit2[[1]]

#####

# 
# TODO FROM PHILINAS CHAT:
#   
#   - [ ] different ranges for spatial and st fields
#   - [ ] keep prior for sigma
#   - [ ] soften mesh barrier (0.3-0.4?)
