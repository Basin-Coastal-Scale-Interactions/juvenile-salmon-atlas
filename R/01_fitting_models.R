# Playing with sdmTMB for juvenile salmon index

library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)
library(here)
library(ggsidekick)
library(tidyverse)
library(lubridate)
library(scales)
library(rnaturalearth)
library(future)
library(furrr)

dat_in <- readRDS(here::here("data-raw", "catch_survey_sbc.rds")) 

# downscale data and predictive grid
dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    month_f = month(date, label = TRUE),
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


## plotting color palette
# col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
# names(col_pal) <- c('chinook','pink','chum','coho','sockeye')


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = 5) # one per species


## mesh shared among species
dat_coords <- dat %>% 
  filter(species == "chinook") %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()

inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(2, 10) * 500,
  cutoff = 30,
  offset = c(10, 50)
)  
# plot(inla_mesh_raw)

# Cam's original mesh
spde <- make_mesh(
  dat %>% 
    filter(species == "chinook"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 
# plot(spde)
spde$mesh$n

# New mesh with more knots and denser around data rich sounds
spde2 <- make_mesh(
  dat %>% 
    filter(species == "chinook"),
  c("utm_x_1000", "utm_y_1000"),
  type = "kmeans",
  n_knots = 100)
spde2$mesh$n
plot(spde2)

png(here("figs", "99-mesh-comparison.png"), width = 8, height = 5, units = "in",
    res = 300)
par(mfrow=c(1,2))
plot(spde)
mtext(paste0("knots = ",spde$mesh$n), 1, cex =2)
plot(spde2)
mtext(paste0("knots = ",spde2$mesh$n), 1, cex =2)
par(mfrow=c(1,1))
dev.off()

dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest()

dat_tbl

# Extracting chinook data first to be able to play with code more easily
chinook_dat <- dat_tbl[[1,"data"]][[1]]

# chinook_pred <- dat_tbl[[1,"sims"]][[1]]
#saveRDS(chinook_pred, here::here("data", "chinook_pred.rds"))
#chinook_pred <- readRDS(here::here("data", chinook_pred.rds"))        


# dir.create("data/", recursive = TRUE, showWarnings = FALSE)
saveRDS(chinook_dat, here::here("data", "chinook_dat.rds"))
saveRDS(spde, here::here("data", "spde_mesh.rds"))
saveRDS(spde2, here::here("data", "spde2_mesh.rds"))

chinook_dat <- readRDS(here::here("data", "chinook_dat.rds"))        
spde <- readRDS(here::here("data","spde_mesh.rds"))        
spde2 <- readRDS(here::here("data","spde2_mesh.rds"))        


# sdmTMB model run
mout3 <- sdmTMB(
  n_juv ~ 0 + s(month, bs = "cc", k = 12) +
    day_night + survey_f +
    scale_dist + scale_depth + (1 | year_f),
  offset = "effort",
  knots = list(month = c(0, 12)),
  # first and last knot are tied together
  # Dec is also month 0
  data = chinook_dat,
  mesh = spde2,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  time = "month",
  extra_time = c(1, 4),
  spatiotemporal = "AR1",
  anisotropy = TRUE,
  silent = FALSE
)

sanity(mout3)
mout3

saveRDS(mout3, here::here("data", "fits", "mout3.rds"))
# mout3 <- readRDS(here::here("data", "fits", "mout3.rds")) 

# pararellizing species models 
fits_list_mout3 <- future_map(
  dat_tbl$data,
  function(x) {
    sdmTMB(n_juv ~ 0 + s(month, bs = "cc", k = 12) + 
             day_night + survey_f + scale_dist + scale_depth + (1 | year_f),
           offset = "effort",
           knots = list(month = c(0, 12)),
           data = x,
           mesh = spde2,
           family = sdmTMB::nbinom2(),
           spatial = "on",
           time = "month",
           extra_time = c(1,4), 
           spatiotemporal = "AR1",
           anisotropy = TRUE,
           silent = FALSE
    )
  }
)

sanity_check <- map(fits_list_mout3, function (x) unlist(sanity(x)) %>% all()) %>% unlist

dat_tbl_mout3 <- dat_tbl %>%
  mutate(fit = fits_list_mout3,
         pass_sanity = sanity_check)
dat_tbl_mout3


sims_mout3 <- predict(mout3, nsim = 50) #

dat_tbl_mout3$sims <- purrr::map2(
  dat_tbl_mout3$fit, dat_tbl_mout3$data, 
  ~ simulate(.x, .y, nsim = 50, , seed = 42, params = c("mle"))
)

# simulate(dat_tbl_mout3$fit[[1]], dat_tbl_mout3$data[[1]], nsim = 1, seed = 42, params = c("mle"))

saveRDS(dat_tbl_mout3, here::here("data", "fits", "all_mout3.rds"))
#dat_tbl_mout3 <- readRDS(here::here("data", "fits", "all_mout3.rds"))


## SPATIAL GRID  ---------------------------------------------------------------

# predictive grid
grid_list <- readRDS(here::here("data-raw", "spatial", "pred_ipes_grid.RDS")) %>% 
  purrr::map(
    .,
    ~ {.x %>% 
        mutate(utm_x_1000 = X / 1000,
               utm_y_1000 = Y / 1000,
               dist_to_coast_km = shore_dist / 1000)}
  )

summer_grid <- grid_list$ipes_grid
fall_grid <- grid_list$wcvi_grid

mutate(summer_grid, season_grid = "summer") %>%
  rbind(., mutate(fall_grid, season_grid = "fall")) %>%
ggplot(aes(X, Y, fill = season_grid)) +
         geom_tile(alpha = 0.6) +
  scale_fill_manual(values = c("blue", "red")) + theme_sleek()

# make key so that missing year-season combinations can be indexed
year_season_key <- expand.grid(
  year = unique(dat_in$year) %>% sort(),
  season_f = unique(dat_in$season_f) %>% sort()
) %>%
  arrange(year, season_f)


# remove season_f and only have a single year
pred_grid_list <- tibble(
    survey_f = "hss",
    year = 2019,
    year_f = as.factor(year)) %>% 
  mutate(id = row_number(),
         survey_f = as.factor(survey_f)) %>%
  split(., .$id)


## INDICES ---------------------------------------------------------------------

#add unique years and seasons
index_grid <- pred_grid_list %>%
  purrr::map(., function (x) {
    dum_grid <- fall_grid
    
    dum_grid %>%
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        target_depth = 0,
        day_night = "DAY",
        scale_depth = (target_depth - mean(dat$target_depth)) /
          sd(dat$target_depth),
        scale_dist = (dist_to_coast_km - mean(dat$dist_to_coast_km)) /
          sd(dat$dist_to_coast_km)
      )
  }) %>%
  bind_rows() %>%
  select(-c(depth, slope, shore_dist))

saveRDS(index_grid, here::here("data", "index_grid.rds"))
#index_grid <- readRDS(here::here("data", "index_grid.rds"))

##########


