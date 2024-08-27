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
library(fmesher)
library(st)
library(sf)

#dat_in <- readRDS(here::here("data", "catch_survey_sbc_withnorth.rds")) # models fitted Feb-Mar 2024 
#dat_trim <- readRDS(here::here("data", "cleaned_index_dat_feb_2024.rds")) # not used recently
#dat_in_apr <- readRDS(here::here("data", "catch_survey_sbc_withnorth_Mar_28_2024.rds")) # latest wrangling Mar 28 2024, 
                                                                          # not fitted yet, doesn't have 2023
dat_trim <- readRDS(here::here("data", "cleaned_index_dat_apr_2024.rds")) # latest wrangling Apr 2, 2024, has 2023 

# downscale data and predictive grid
dat <- dat_trim %>% 
  mutate(
    month_f = month(date, label = TRUE),
    month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
    scale_month_adj = scale(month_adj)[, 1],
    year_adj = ifelse(month > 2, year, year - 1), # adjusting year to represent fish cohorts
    year_f = as.factor(year),
    year_adj_f = as.factor(year_adj),
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

# To reverse month_adj scaling:
# dat$scale_month_adj * attr(scale(dat$month_adj),"scaled:scale") + attr(scale(dat$month_adj),"scaled:center")

## plotting color palette
col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')



# prep multisession
ncores <- parallel::detectCores() 
#plan(multisession, workers = 5) # one per species


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
plot(inla_mesh_raw)

# Cam's original mesh
spde <- make_mesh(
  dat %>% 
    filter(species == "chinook"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 
plot(spde)
spde$mesh$n

# New mesh with more knots and denser around data rich sounds
# adding ~90 more knots for a total of 341 to see if models run alright
spde2 <- make_mesh(
  dat %>% 
    filter(species == "chinook"),
  c("utm_x_1000", "utm_y_1000"),
  type = "kmeans",
  n_knots = 210)
spde2$mesh$n
plot(spde2)

st_dc <- st_multipoint(dat_coords)

max_edge <- diff(range(st_dc[,1]))/(3*5)
bound_outer <- diff(range(st_dc[,1]))/3

mesh4 <-  fm_mesh_2d_inla(
  loc = dat_coords,
  cutoff = 21,
  max.edge = c(1,2) * max_edge#,
  #offset = c(max_edge, bound_outer)
)
mesh4
plot(mesh4)



# larger offset
mesh5 <-  fm_mesh_2d_inla(
  loc = dat_coords,
  cutoff = 19,
  max.edge = c(1, 2) * max_edge,
  offset = c(max_edge, bound_outer)
)


ggplot() +
  inlabru::gg(mesh5) 
ggplot() +
  inlabru::gg(mesh4) 
mesh4$n
mesh5$n

dmesh <- make_mesh(dat %>% 
                     filter(species == "chinook"),
                   c("utm_x_1000", "utm_y_1000"), mesh = mesh4)
dmesh$mesh$n

plot(dmesh)
plot(all_coast_km)

plot_mesh <- function (x, ...) 
{
  plot(x$mesh, main = NA, edge.color = "grey60", asp = 1, ...)
  points(x$loc_xy, pch = 20, cex = 0.8, col = "red") #col = "#00000080")
  #points(x$loc_centers, pch = 20, col = "red")
}

plot_mesh(dmesh)



png(here("figs", "99-mesh-comparison-allcoast.png"), width = 8, height = 8, units = "in",
    res = 300)
par(mfrow=c(2,2))
plot(spde)
mtext(paste0("knots = ", spde$mesh$n), 1, cex = 2)
plot(spde2)
mtext(paste0("knots = ", spde2$mesh$n), 1, cex = 2)
plot(dmesh)
mtext(paste0("knots = ", dmesh$mesh$n), 1, cex = 2)
par(mfrow=c(1,1))
dev.off()

dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest()

dat_tbl



# Extracting chinook data first to be able to play with code more easily
chinook_dat <- dat_tbl[[1,"data"]][[1]]
glimpse(chinook_dat)

# all months except Jan
table(chinook_dat$month)


# dir.create("data/", recursive = TRUE, showWarnings = FALSE)
saveRDS(chinook_dat, here::here("data", "chinook_dat_allcoast.rds"))
saveRDS(spde, here::here("data", "spde_mesh_allcoast.rds"))
saveRDS(spde2, here::here("data", "spde2_mesh_allcoast.rds"))
saveRDS(dmesh, here::here("data", "dense_mesh_allcoast.rds"))

saveRDS(dat_tbl, here::here("data", "dat_tbl_allcoast.rds"))


chinook_dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))        
spde <- readRDS(here::here("data","spde_mesh_allcoast.rds"))        
spde2 <- readRDS(here::here("data","spde2_mesh_allcoast.rds"))        
dmesh <- readRDS(here::here("data", "dense_mesh_allcoast.rds"))

# sdmTMB model run
# mout3 <- sdmTMB(
#   n_juv ~ 0 + s(month_adj, bs = "tp", k = 6) +
#     day_night + survey_f +
#     #scale_dist + 
#     scale_depth + (1 | year_f),
#   offset = "effort",
#   knots = list(month = c(0, 12)),
#   # first and last knot are tied together
#   # Dec is also month 0
#   data = chinook_dat,
#   mesh = spde2,
#   family = sdmTMB::nbinom2(),
#   spatial = "off",
#   time = "month_adj",
#   extra_time = c(2, 11),
#   spatiotemporal = "RW",
#   anisotropy = TRUE,
#   silent = FALSE
# )
# sanity(mout3)
# mout3


# sdmTMB model run
mdmesh <- sdmTMB(
  n_juv ~ 0 + s(month_adj, bs = "tp", k = 6) +
    day_night + survey_f +
    #scale_dist +
    scale_depth + (1 | year_f),
  offset = "effort",
  knots = list(month = c(0, 12)),
  # first and last knot are tied together
  # Dec is also month 0
  data = chinook_dat,
  mesh = dmesh,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  time = "month_adj",
  extra_time = c(2, 11),
  spatiotemporal = "RW",
  anisotropy = TRUE,
  silent = FALSE
)
# 
# sanity(mdmesh)
# mdmesh

# saveRDS(mdmesh, here::here("data", "fits", "dmesh_allcoast.rds"))
# mdmesh <- readRDS(here::here("data", "fits", "dmesh_allcoast.rds"))

# checking if Hessian function is available in sdmTMB fit object
# obj <- mdmesh$tmb_obj
# obj$fn() # NEED TO RERUN THE FUNCTION TO REATTACH THE ENVIRONMENT, OTHERWISE SPHESS FUNCTION CRASHES!
# obj$env$spHess() %>% dim()
# 
# saveRDS(mout3, here::here("data", "fits", "mout3_allcoast.rds"))
# mout3 <- readRDS(here::here("data", "fits", "mout3_allcoast.rds")) 


RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
get_num_cores()
get_num_procs()
blas_get_num_procs()
omp_get_num_procs()
omp_get_max_threads()

# pararellizing species models 
fits_list_mout3 <- map( #future_map(
  dat_tbl$data,
  function(x) {
    sdmTMB(n_juv ~ 0 + s(scale_month_adj, bs = "tp", k = 6) + # tp smoother is default 
             day_night + survey_f + 
             # scale_dist + 
             scale_depth + (1 | year_adj_f),
           offset = "effort",
           # knots = list(month = c(0, 12)),
           data = x,
           mesh = spde2,
           family = sdmTMB::nbinom2(),
           spatial = "off", # No need to estimate "average"/spatial only RF, 
           #tends to collapse with RW/also when first time step is flat
           time = "month_adj",
           extra_time = 11, #c(2,11),  # month_adj 2 is April and 11 is January
           spatiotemporal = "RW", # RW = AR1 with cor of 1, Philina said only 
           #viable option with temporally patchy data (Greenland Halibut analysis)
           anisotropy = TRUE,
           silent = FALSE
    )
  }
)

names(sanity_check_mout3) <-dat_tbl$species

saveRDS(fits_list_mout3, here::here("data", "fits", "fits_list_mout3.rds"))


# pararellizing species models 
fits_list_mdmesh <- map( #future_map(
  dat_tbl$data,
  function(x) {
    sdmTMB(n_juv ~ 0 + s(scale_month_adj, bs = "tp", k = 6) + # tp smoother is default 
             day_night + survey_f + 
             # scale_dist + 
             scale_depth + (1 | year_adj_f),
           offset = "effort",
           # knots = list(month = c(0, 12)),
           data = x,
           mesh = dmesh,
           family = sdmTMB::nbinom2(),
           spatial = "off", # No need to estimate "average"/spatial only RF, 
           #tends to collapse with RW/also when first time step is flat
           time = "month_adj",
           extra_time = 11, #c(2,11),  # month_adj = 2 is April and 11 is January, 1 is March (start of cycle)
           spatiotemporal = "RW", # RW = AR1 with cor of 1, Philina said only 
           #viable option with temporally patchy data (Greenland Halibut analysis)
           anisotropy = TRUE,
           silent = FALSE
    )
  }
)
sanity_check_dmesh <- map(fits_list_mdmesh, function (x) unlist(sanity(x)) %>% all()) %>% unlist
names(fits_list_mdmesh) <-dat_tbl$species

saveRDS(fits_list_mdmesh, here::here("data", "fits", "fits_list_mdmesh.rds"))

# fits_extra_pink <- run_extra_optimization(fits_list_mout3[[4]], nlminb_loops = 3, newton_loops = 3)
# sanity(fits_extra_pink)

dat_tbl_dmesh <- dat_tbl %>%
  mutate(fit = fits_list_mdmesh,
         pass_sanity = sanity_check_dmesh)
dat_tbl_dmesh

dat_tbl_mout3 <- dat_tbl %>%
  mutate(fit = fits_list_mout3,
         pass_sanity = sanity_check_mout3)




saveRDS(dat_tbl_dmesh, here::here("data", "fits", "all_dmesh_allcoast.rds"))
saveRDS(dat_tbl_mout3, here::here("data", "fits", "all_mout3_allcoast.rds"))
#dat_tbl_mout3 <- readRDS(here::here("data", "fits", "all_mout3_allcoast.rds"))
