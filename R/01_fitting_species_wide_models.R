# Playing with sdmTMB for juvenile salmon index

library(sdmTMB)
library(ggplot2)
#library(sdmTMBextra)
library(here)
library(ggsidekick)
library(tidyverse)
library(scales)
library(rnaturalearth)
library(fmesher)
library(st)
library(sf)

source("R/plot_smoothers.R")

# latest wrangling has all of 2024
dat_in <- readRDS(here::here("data", "cleaned_atlas_catch_dat.rds")) 

# downscale data and predictive grid
dat <- dat_in %>% 
  mutate(
    month_f = month(date, label = TRUE),
    month_adj = case_when(
      # adjusting to start in April except chum which starts in March      
     species == "chum" ~ ifelse(month > 2, month - 2, month + 10), 
     species %in% c("coho", "pink","sockeye","chinook") ~ ifelse(month > 3, month - 3, month + 9)
     ),
    scale_month_adj = scale(month_adj)[, 1],
    year_adj = case_when(
      # adjusting year to represent fish cohorts
      species == "chum" ~ ifelse(month > 2, year, year - 1),
      species %in% c("coho", "pink","sockeye","chinook") ~ ifelse(month > 3, year, year - 1)
    ),
    year_f = as.factor(year),
    year_adj_f = as.factor(year_adj),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3),
   # scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night)) %>% 
  filter(!species == "steelhead") %>% 
  droplevels() %>%
  as_tibble() %>%
  # temporarily fixing NAs in day_night here (should be fixed in db)
  mutate(day_night = case_when(
    is.na(day_night) & unique_event == "IPES2024-026-126" ~ "NIGHT",
    is.na(day_night) & month == 10 ~ "DAY",
    .default = day_night
  ))
  
filter(dat, is.na(day_night))


# checking that month_adj makes sense
dat %>%
  filter(species == "chinook") %>%
  group_by(month_adj) %>%
  summarise(month = unique(month_f)) %>%
  as_tibble()

# looking at unique events around Barkley Sound 
barkley_ue <- filter(dat, lat > 48, lat < 49, 
                     lon < -125, lon > -126, 
                     species == "chinook") %>% 
  select(unique_event, date) %>%
  mutate(date = as_date(date)) %>%
  as_tibble()

saveRDS(barkley_ue, file = here("data", "barkley_unique_events.rds"))

# To reverse month_adj scaling:
# dat$scale_month_adj * attr(scale(dat$month_adj),"scaled:scale") + attr(scale(dat$month_adj),"scaled:center")

## plotting color palette
col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')

# prep multisession
# ncores <- parallel::detectCores() 
# plan(multisession, workers = 5) # one per species

## mesh shared among species
dat_coords <- dat %>% 
  filter(species == "chinook") %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()

# custom mesh:
domain <- fmesher::fm_nonconvex_hull_inla(
  dat_coords,
  # concave = -0.07, convex = -0.05, resolution = c(200, 200)
  #concave = 1, convex = 1, resolution = c(200, 200)
 concave = -4, convex = -0.02, resolution = c(500, 500)
)
glimpse(domain)
plot(domain)

## started with these: fit with regular lognormal_mix, but not with poisson-link
# min_edge <- 40 #try this as model doesn't converge
# max_edge <- 60
## so tried a finer mesh: fit for everything!
min_edge <- 21
max_edge <- 60

# if (file.exists("output/mesh-trawl-overall.rds")) {
#   # ensure consistent knot locations across platforms:
#   mesh3 <- readRDS("output/mesh-trawl-overall.rds")
# } else {
  mesh3 <- fm_mesh_2d_inla(
    loc = dat_coords,
    # boundary = domain,
    # interior = domain,
    # max.edge = c(150, 2000),
    max.edge = c(max_edge, 2*max_edge),
    offset = c(10, 70),
    cutoff = min_edge
  )
  #saveRDS(mesh3, "output/mesh-trawl-overall.rds")
#}

chmesh <- make_mesh(dat %>% filter(species == "chinook"), 
                    c("utm_x_1000", "utm_y_1000"), mesh = mesh3)
chmesh$mesh$n
plot(chmesh)



inla_mesh_raw <- fmesher::fm_mesh_2d_inla(
  loc = dat_coords,
  max.edge = c(2, 10) * 500,
  cutoff = 30,
  offset = c(10, 50)
)  
plot(inla_mesh_raw)

# CF's original mesh
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

mesh4 <-  fmesher::fm_mesh_2d_inla(
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
plot(mesh5)

# ggplot() +
#   inlabru::gg(mesh5) 
# ggplot() +
#   inlabru::gg(mesh4) 
mesh4$n
mesh5$n

dmesh <- make_mesh(dat %>% 
                     filter(species == "chinook"),
                   c("utm_x_1000", "utm_y_1000"), mesh = mesh4)
dmesh$mesh$n
plot(dmesh)
chmesh$mesh$n
plot(chmesh)

plot_mesh <- function (x, ...) 
{
  plot(x$mesh, main = NA, edge.color = "grey60", asp = 1, ...)
  points(x$loc_xy, pch = 20, cex = 0.8, col = "red") #col = "#00000080")
  #points(x$loc_centers, pch = 20, col = "red")
}

plot_mesh(dmesh)
plot_mesh(chmesh)


png(here("figs", "99-mesh-comparison-allcoast.png"), width = 8, height = 8, units = "in",
    res = 300)
par(mfrow=c(2,2))
plot(spde)
mtext(paste0("knots = ", spde$mesh$n), 1, cex = 2)
plot(spde2)
mtext(paste0("knots = ", spde2$mesh$n), 1, cex = 2)
plot(dmesh)
mtext(paste0("knots = ", dmesh$mesh$n), 1, cex = 2)
plot(chmesh)
mtext(paste0("knots = ", chmesh$mesh$n), 1, cex = 2)
par(mfrow=c(1,1))
dev.off()

dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest()

dat_tbl

# Extracting chinook data first to be able to play with code more easily
chinook_dat <- dat_tbl[[1,"data"]][[1]]
glimpse(chinook_dat)

chinook_dat %>%
  group_by(month_adj) %>%
  summarise(month = unique(month_f))





# all months except Jan
table(chinook_dat$month)


# looking at NAs in the data
filter(chinook_dat, is.na(day_night)) %>% print(n = 65)
filter(chinook_dat, year == 2024 & month == 10) %>% print(n = 65)
filter(chinook_dat, is.na(day_night) & month == 7)$unique_event




# dir.create("data/", recursive = TRUE, showWarnings = FALSE)
saveRDS(chinook_dat, here::here("data", "chinook_dat_allcoast.rds"))
saveRDS(spde, here::here("data", "spde_mesh_allcoast.rds"))
saveRDS(spde2, here::here("data", "spde2_mesh_allcoast.rds"))
saveRDS(dmesh, here::here("data", "dense_mesh_allcoast.rds"))
saveRDS(chmesh, here::here("data", "chmesh_allcoast.rds"))
saveRDS(dat_tbl, here::here("data", "dat_tbl_allcoast.rds"))

chinook_dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))        
spde <- readRDS(here::here("data","spde_mesh_allcoast.rds"))        
spde2 <- readRDS(here::here("data","spde2_mesh_allcoast.rds"))        
dmesh <- readRDS(here::here("data", "dense_mesh_allcoast.rds"))
chmesh <- readRDS(here::here("data", "chmesh_allcoast.rds"))

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
mchmesh <- sdmTMB(
  n_juv ~ 0 + s(scale_month_adj, bs = "tp", k = 5) +
    day_night + survey_f +
    #scale_dist +
    scale_depth + 
    (1 | year_adj_f),
  offset = "effort",
 # knots = list(month = c(0, 12)),
  # first and last knot are tied together
  # Dec is also month 0
  data = chinook_dat,
  mesh = chmesh,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  time = "month_adj",
  extra_time = 1:12, #c(10),
  spatiotemporal = "RW",
  anisotropy = TRUE,
  #control = sdmTMBcontrol(map = list(#ln_tau_Z = factor(c(rep(1, 10), rep(2, 10)))
 #                                    #b_smooth = factor(b_smooth_map),
 #                                    #ln_smooth_sigma = factor(c(1,2,3,4,5,6,7,8,NA)
                                     #ln_smooth_sigma = factor(NA))),
 silent = FALSE
)



sanity(mchmesh)
mchmesh
tidy(mchmesh, "ran_pars", conf.int = TRUE)
mchmesh$sd_report
mchmesh$tmb_obj$env$spHess() %>% dim()

psmooth_speciesch <- plot_smoother(mchmesh, species = "chinook")
psmooth_speciesch + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + ggtitle("chmesh")


plot_smoother(mchmesh) + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

### other mesh
mdmesh <- sdmTMB(
  n_juv ~ 0 + s(scale_month_adj, bs = "tp", k = 4) +
    day_night + survey_f +
    #scale_dist +
    scale_depth + 
    (1 | year_adj_f),
  offset = "effort",
  # knots = list(month = c(0, 12)),
  # first and last knot are tied together
  # Dec is also month 0
  data = chinook_dat,
  mesh = chmesh,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  time = "month_adj",
  extra_time = 1:12, #c(10),
  spatiotemporal = "RW",
  anisotropy = TRUE,
  #control = sdmTMBcontrol(map = list(#ln_tau_Z = factor(c(rep(1, 10), rep(2, 10)))
  #                                    #b_smooth = factor(b_smooth_map),
  #                                    #ln_smooth_sigma = factor(c(1,2,3,4,5,6,7,8,NA)
  #ln_smooth_sigma = factor(NA))),
  silent = FALSE
)

sanity(mdmesh)
summary(mdmesh)
tidy(mdmesh, "ran_pars", conf.int = TRUE)
mdmesh$sd_report
mdmesh$tmb_obj$env$spHess() %>% dim()

psmooth_speciesd <- plot_smoother(mdmesh, species = "chinook")
psmooth_speciesd + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + ggtitle("dmesh")

# saveRDS(mdmesh, here::here("data", "fits", "dmesh_allcoast.rds"))
# mdmesh <- readRDS(here::here("data", "fits", "dmesh_allcoast.rds"))

# checking if Hessian function is available in sdmTMB fit object
# obj <- mdmesh$tmb_obj
# obj$fn() # NEED TO RERUN THE FUNCTION TO REATTACH THE ENVIRONMENT, OTHERWISE SPHESS FUNCTION CRASHES!
# obj$env$spHess() %>% dim()

# RhpcBLASctl::blas_set_num_threads(1)
# RhpcBLASctl::omp_set_num_threads(1)
# get_num_cores()
# get_num_procs()
# blas_get_num_procs()
# omp_get_num_procs()
# omp_get_max_threads()

# pararellizing species models 
fits_list_mdmesh <- map2( #future_map(
  dat_tbl$data, dat_tbl$species,
  function(x, y) {
    sdmTMB(n_juv ~ 0 + s(scale_month_adj, bs = "tp", k = 4) + # tp smoother is default 
             day_night + survey_f + 
             # scale_dist + 
             scale_depth + (1 | year_adj_f),
           offset = "effort",
           # knots = list(month = c(0, 12)),
           data = x,
           mesh = chmesh,
           family = sdmTMB::nbinom2(),
           spatial = "off", # No need to estimate "average"/spatial only RF, 
           #tends to collapse with RW/also when first time step is flat
           time = "month_adj",
           extra_time = 1:12, #11, #c(2,11),  # month_adj = 2 is April and 11 is January, 1 is March (start of cycle)
           spatiotemporal = "RW", # RW = AR1 with cor of 1, Philina said only 
           #viable option with temporally patchy data (Greenland Halibut analysis)
           anisotropy = TRUE,
           # fixing ln_smooth_sigma for chinook
          # control = sdmTMBcontrol(map = list(ln_smooth_sigma = factor(if_else(y == "chinook", NA, 1)))),
           silent = FALSE
    )
  }
)

sanity_check_dmesh <- map(fits_list_mdmesh, function (x) unlist(sanity(x)) %>% all()) %>% unlist
names(fits_list_mdmesh) <- dat_tbl$species

saveRDS(fits_list_mdmesh, here::here("data", "fits", "fits_list_mdmesh.rds"))

dat_tbl_dmesh <- dat_tbl %>%
  mutate(fit = fits_list_mdmesh,
         pass_sanity = sanity_check_dmesh)

dat_tbl_dmesh

saveRDS(dat_tbl_dmesh, here::here("data", "fits", "all_dmesh_allcoast.rds"))

map2(fits_list_mdmesh, dat_tbl_dmesh$species,
     function(x,y) plot_smoother(x, species = y) + ggtitle(y) +
       theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))

