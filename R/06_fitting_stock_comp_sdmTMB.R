
library(tidyverse)
library(sdmTMB)
#library(sdmTMBextra)
#library(mgcv)
library(here)
library(sf)
source("R/plot_map.R")
source("R/cleave_by.R")

#  Load GSI data ------------------------------------------------------------

# doesn't include fish with max stock_prop <0.5
dat <- readRDS(here::here("data", "chinook_gsi_counts.rds")) %>% 
  select(-month_adj) %>%
  mutate(month_adj = ifelse(month > 3, month - 3, month + 9), # redo month_adj to start in April
        # season_n = as.numeric(season_f),
         scale_month_adj = scale(month_adj)[, 1],
         year_adj = ifelse(month > 3, year, year - 1), # adjusting year to represent fish cohorts
         # year_f = as.factor(year),
         year_adj_f = as.factor(year_adj))

# 19 tows above 57 N, all in nov and dec, mostly 2000-2002
filter(dat, lat > 57) %>%
  group_by(unique_event) %>%
  summarise(date = unique(date), totalfish = unique(n_fish), lat = unique(lat))

table(filter(dat, lat > 57)$month, filter(dat, lat > 57)$year)
range(dat$lat)

table(dat$region)
levels(dat$region)

# removing tows above 57 N latitude
dat_trim <- filter(dat, lat < 57)

# saving df used for model fititng
saveRDS(dat_trim, here::here("data", "chinook_gsi_counts_fitted.rds"))

dat_trim %>% group_by(region) %>% summarize(propsum = sum(stock_prop))

# only keeping four most common MUs
# dat <- filter(dat, region %in% c("WCVI", "North/Central BC", "Columbia")) %>%
#   droplevels()

glimpse(dat)
levels(dat$region)

# all observations
dat %>%
  filter(region == "WCVI") %>%
ggplot(aes(utm_x, utm_y, colour = stock_prop)) +
  geom_point() +
  coord_fixed() + 
  facet_wrap(~month)

ggplot() +
  geom_sf(data = all_coast, color = "gray80", fill = "gray90") +
  geom_point(data = filter(dat, region == "WCVI"), aes(utm_x, utm_y)) +
  xlim(range(dat$utm_x)) + ylim(range(dat$utm_y))

# trimmed observations < 57 N
ggplot() +
  geom_sf(data = all_coast, color = "gray80", fill = "gray90") +
  geom_point(data = filter(dat_trim, region == "WCVI"), aes(utm_x, utm_y)) +
  xlim(range(dat_trim$utm_x)) + ylim(range(dat_trim$utm_y)) +
  theme_bw()

ggplot(dat_trim, aes(utm_x, utm_y, colour = stock_prop)) +
  geom_point(alpha= 0.4, shape = 16) +
  facet_wrap(~region)

# Mesh and SPDE matrix construction -----------------------------------------

dat_coords <- dat_trim %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()

inla_mesh_raw <- fmesher::fm_mesh_2d_inla(
  loc = dat_coords,
  max.edge = c(60,100),
  cutoff = 20,
  offset = c(20, 50),
)  
plot(inla_mesh_raw)

sdmTMB_mesh <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), mesh = inla_mesh_raw)  
#readRDS(data/gsi-prop-mesh.rds")
spde <- sdmTMB_mesh$spde
mesh <- sdmTMB_mesh$mesh
inneridx <- mesh$segm$int$idx # index of points around inner boundary
innerbox <- mesh$loc[mesh$segm$int$idx[,1],-3] # points around inner edge of mesh (for overlaying to predictive grid)

plot(sdmTMB_mesh)
sdmTMB_mesh$mesh$n

saveRDS(sdmTMB_mesh, file = "data/gsi-prop-mesh.rds")

# New mesh with more knots and denser around data rich sounds
# spde2 <- make_mesh(
#   dat,
#   c("utm_x_1000", "utm_y_1000"),
#   type = "kmeans",
#   n_knots = 210)
# spde2$mesh$n
# plot(spde2)

st_dc <- sf::st_multipoint(dat_coords)

max_edge <- diff(range(st_dc[,1]))/(3*5)
bound_outer <- diff(range(st_dc[,1]))/3

mesh4 <-  fmesher::fm_mesh_2d_inla(
  loc = dat_coords,
  cutoff = 21,
  max.edge = c(1, 2) * max_edge#,
  #offset = c(max_edge, bound_outer)
)
mesh4
spde4 <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), mesh = mesh4) 
plot(spde4)

# compute bilinear interpolation matrix from mesh to data:
# interpolator_data <- fmesher::fm_basis(
#   mesh,
#   loc = as.matrix(dat[, c("utm_x_1000", "utm_y_1000")])
# )

dat2 <- dat_trim
dat2$region <- fct_rev(dat2$region)

#debugonce(sdmTMB)

table(dat_trim$month_adj)
table(dat_trim$month)
table(dat_trim$region)

m_svc <- sdmTMB(
  stock_prop ~ 0 + region + (1 | year_adj_f) + s(scale_month_adj, by = region, bs = "tp", k = 6),
  family = tweedie(),
  spatial_varying = ~ 0 + region * scale_month_adj,
 # offset = dat_trim$effort,
  data = dat_trim,
  time = "month_adj",
  extra_time = 2:12, 
  spatial = "off",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  mesh = sdmTMB_mesh,
  silent = FALSE,
  control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L),
  do_fit = FALSE
)

# b_smooth_map <-  1:40
# b_smooth_map[17:20] <- NA
# b_smooth_map

# Mapping same SD for both all coefficients and interactions (ln_tau_Z)
m_svc <- update(m_svc, do_fit = TRUE,
                control = sdmTMBcontrol(map = list(ln_tau_Z = factor(c(rep(1, 11), rep(2, 11))),
                                                  # b_smooth = factor(b_smooth_map),
                                                  ln_smooth_sigma = factor(rep(1, 11))
                                                  # ln_smooth_sigma = factor(c(1,2,3,4,NA,6,7,8,9,10))
                                                   )
                                        )
                )



sanity(m_svc)

m_svc
m_svc$version
m_svc$sd_report
tidy(m_svc, "ran_pars", conf.int = TRUE)

saveRDS(m_svc, file = here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))
#m_svc <- readRDS(here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))
m_svc$formula
m_svc$time
m_svc$extra_time
# this shows the 20 distinct parameters relevant to ln_tau_Z, which are the
# ten regions plus ten interactions of each region with scaled month, which is continuous
m_svc$spatial_varying

m_svc$data$region %>% table


### PREDICTING 

index_grid <- readRDS(here("data", "pred_grid_all_coast.rds")) %>%
  # filter(row_number() %% 5 == 1) %>% ### DOWNSAMPLING GRID %>%
  # arrange(X,Y) %>%
  # filter(row_number() %% 5 == 1) %>% ### DOWNSAMPLING GRID %>%
  mutate(year = 2019L, 
         year_f = as.factor(year)
        # survey_f = as.factor("hss"),
         #target_depth = 0,
         #scale_depth = -0.5427145,
        # day_night = "DAY"
         ) %>%
  rename(year_adj_f = year_f, year_adj = year)

glimpse(index_grid)

ggplot(index_grid) +
  geom_tile(aes(X,Y, fill = depth))

# Polygon of points from inner sdmTMB mesh
ib_poly <- st_polygon(list(inner = rbind(innerbox,innerbox[1,])))
st_crs(ib_poly)

grid_sf <- index_grid %>% 
  st_as_sf(., coords = c("utm_x_1000", "utm_y_1000"), remove = FALSE)

# extract grid points that are within inner boundary of mesh
prop_grid <- st_intersection(grid_sf, ib_poly) %>% sf::st_drop_geometry()
glimpse(prop_grid)


ggplot(prop_grid) +
  geom_tile(aes(X,Y, fill = depth))

dim(index_grid)
dim(prop_grid)

saveRDS(prop_grid, file = here("data", "gsi-prop-grid.rds"))
# check
# ggplot() +
#   geom_sf(data = st_intersection(bridge_sf, ipes_sf_poly))


# prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))
# m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))
# dat_trim <- readRDS(here::here("data", "chinook_gsi_counts_fitted_2025-02-03.rds"))

# functions for quickly applying and unapplying scale()
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

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


# Splitting grid data frame into manageable chunks for prediction
grid_region <- grid_months %>% cleave_by(region, month_adj)
#grid_region <- grid_months %>% group_split(group_id = region) 

length(grid_region)

# reduce memory usage
#rm(grid_months)
gc()

# predtest <- predict(m_svc, newdata = grid_region[[1]], nsim = 3,
#                 offset = rep.int(median(dat_trim$effort), nrow(grid_region[[1]])))
# glimpse(predtest)
# 
# simtest <- simulate(m_svc, nsim = 3, type = "mle-mvn",
#                     offset = rep.int(median(dat_trim$effort), nrow(grid_region[[1]])))
# glimpse(simtest)
# head(simtest)


# predicting on grid chunks for each region
predict_list <- map(grid_region, function (newdata) {
  
  # max_chunk_length <- 5000 
  # smaller_chunks <- group_split(newdata, group_id = row_number() %/% max_chunk_length)
  # pred_list <- map(smaller_chunks,
  #                  ~predict(m, newdata = .x, se_fit = FALSE, re_form_iid = NA, 
  #                           offset = rep.int(median(dat_trim$effort), nrow(.x))))
  
  gc()
  pred <- predict(m_svc, newdata = newdata, se_fit = FALSE, re_form_iid = NA)
  
  group_id <- as.character(unique(pred$region))
  print(group_id)

  pred$pred_i <- pred$est %>% m_svc$family$linkinv() #pred_linkinv
  pred
})


# creating objects for converting predictions to proportions
pred_ic <- array(NA, dim = c(nrow(grid_months)/nlevels(grid_months$region), 
                             nlevels(grid_months$region)),
                 dimnames = list(NULL, levels(grid_months$region)))
# se_pred_ic <- pred_ic
dim(pred_ic)

# joining prediction chunks and re-splitting predictions by region
pred <- bind_rows(predict_list)
predict_region <- group_split(pred, group_id = region) 
names(predict_region) <- levels(pred$region)

gc()

for (i in seq_along(levels(grid_months$region))) {
  cur_region <- levels(grid_months$region)[i]
  pred_ic[,cur_region] <- predict_region[[i]]$pred_i
}
head(pred_ic)

# Normalize probability for each observation and class
rowsum_pred_ic <- outer(rowSums(pred_ic), rep(1, ncol(pred_ic)))
head(rowsum_pred_ic)

prob_ic <- pred_ic / rowsum_pred_ic
head(prob_ic)

#prob_i <- prob_ic[cbind(1:nrow(pred_ic), match(grid_months$region, levels(dat$region)))]

# return prediction
#if(se.fit==TRUE) {
# Normalize SE-squared for each observation and class
# rowsum_se2_ic = outer( rowSums(se_pred_ic^2), rep(1,ncol(pred_ic)) )
# se2_prob_ic = prob_ic^2 * ( se_pred_ic^2/pred_ic^2 - 2*se_pred_ic^2/(pred_ic*rowsum_pred_ic) + rowsum_se2_ic/rowsum_pred_ic^2 )
# se_i = sqrt(se2_prob_ic[ cbind(1:nrow(se2_prob_ic), match(newdata$region,levels(dat$region))) ])
# out = list("fit"=prob_i, "se.fit"=se_i)
# } else {
#   out = prob_i


# reassigning probabilities back to prediction list
# for (i in seq_along(levels(grid_months$region))) {
#   cur_region <- levels(grid_months$region)[i]
#   predict_region[[cur_region]]$prob_i <-
# }

pred_region2 <- map(predict_region, function (x) {
  cur_region <- unique(x$region)
  print(cur_region)
  x$prob_i <- prob_ic[, cur_region]
  x2 <- select(x, -starts_with("zeta_s"))
  return(x2)
})


rm(predict_list)
gc()

pred2 <- bind_rows(pred_region2)


head(pred)
head(pred2)
head(pred_ic)
head(rowsum_pred_ic)
head(prob_ic)

head(pred2$prob_i)
length(pred2$prob_i)

pred2 %>%
  ggplot(data = .) +
  geom_histogram(aes(prob_i)) +
  facet_grid(region ~ month) 

glimpse(predict_list)
#glimpse(prob_i)

hist(pred$pred_i)
hist(pred$prob_i)

# new predict_list 
predict_list <- pred2 %>% cleave_by(region, month_adj)


saveRDS(pred2, file = "data/fits/gsi-prediction-normalized-svc-sdmTMB-no05prop.rds")
saveRDS(predict_list, file = "data/fits/gsi-prediction-normalized-svc-sdmTMB-list-no05prop.rds")

# pred2 <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB-no05prop.rds")

source("R/plot_map.R")

pred_ggdata <- pred2 %>%
  dplyr::select(utm_x_1000, utm_y_1000, prob_i, salmon_region = region, month_f)

p <- ggplot(data = pred_ggdata) + # (dat, aes(X, Y, color = {{ column }})) +
  geom_raster( aes(x=utm_x_1000, y=utm_y_1000, fill = prob_i)) +
  scale_fill_viridis_c(name = "Probability of\nGSI assignment",#"Predicted\nnumber of\njuveniles\nin GSI samples",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  geom_sf(data = all_coast_km, color = "gray80", fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(pred2$utm_x_1000),
                     labels = c("", "132°W", "", "128°W", "", "124°W", ""),
                                expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(pred2$utm_y_1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinoook GSI")


# rm(pred)
# rm(predict_region)

ggsave(p, filename = here("figs", 
                          "chinook_by_gsi_prob_i_svc_mapped_tauZ_sdmTMB_no05prop_fraserfall.png"), 
       width = 14, height = 16)


