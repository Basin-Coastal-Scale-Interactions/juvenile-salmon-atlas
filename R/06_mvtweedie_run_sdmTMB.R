
library(tidyverse)
#library(RTMB)
library(sdmTMB)
library(sdmTMBextra)
library(mgcv)
#library(mvtweedie)
library(here)
library(sf)
#  Load GSI data ------------------------------------------------------------

dat <- readRDS(here::here("data", "chinook_gsi_counts.rds")) %>%
  mutate(season_n = as.numeric(season_f))

table(dat$region)
  
levels(dat$region)

dat %>% group_by(region) %>% summarize(propsum = sum(stock_prop))

# only keeping four most common MUs
# dat <- filter(dat, region %in% c("WCVI", "North/Central BC", "Columbia")) %>%
#   droplevels()

glimpse(dat)
levels(dat$region)

# observations
ggplot(dat, aes(utm_x, utm_y, colour = stock_prop)) +
  geom_point() +
  facet_wrap(~month)

# Mesh and SPDE matrix construction -----------------------------------------

dat_coords <- dat %>% 
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
spde <- make_mesh(dat, c("utm_x_1000", "utm_y_1000"), mesh = inla_mesh_raw) 
plot(spde)
spde$mesh$n

saveRDS(spde, file = "data/gsi-prop-mesh.rds")

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

ggplot() +
  inlabru::gg(mesh4) 

mesh4$n

# dmesh <- make_mesh(dat %>% 
#                      filter(species == "chinook"),
#                    c("utm_x_1000", "utm_y_1000"), mesh = mesh4)

sdmTMB_mesh <- spde  #readRDS(data/gsi-prop-mesh.rds")
spde <- sdmTMB_mesh$spde
mesh <- sdmTMB_mesh$mesh
inneridx <- mesh$segm$int$idx # index of points around inner boundary
innerbox <- mesh$loc[mesh$segm$int$idx[,1],-3] # points around inner edge of mesh (for overlaying to predictive grid)

# plot(mesh)
# plot(sdmTMB_mesh)

# compute bilinear interpolation matrix from mesh to data:
# interpolator_data <- fmesher::fm_basis(
#   mesh,
#   loc = as.matrix(dat[, c("utm_x_1000", "utm_y_1000")])
# )


dat2 <- dat
dat2$region <- fct_rev(dat2$region)

#debugonce(sdmTMB)
m_svc <- sdmTMB(
  stock_prop ~ 0 + region + (1 | year_f) + s(month_adj, by = region, bs = "tp", k = 6), 
  family = tweedie(),
  spatial_varying = ~ 0 + region * month_adj,
  offset = dat$effort,
  data = dat2,
  time = "month_adj",
  extra_time = c(2,11),
  spatial = "off",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  mesh = sdmTMB_mesh,
  silent = FALSE,
  control = sdmTMBcontrol(newton_loops = 1L),
  do_fit = FALSE
)

sanity(m)

m
m$sd_report
tidy(m, "ran_pars", conf.int = TRUE)

saveRDS(m, file = here("data", "fits", "chinook_gsi_prop_svc_int_sdmTMB.rds"))

saveRDS(m, file = here("data", "fits", "chinook_gsi_prop_sdmTMB.rds"))
m <- readRDS(here::here("data", "fits", "chinook_gsi_prop_sdmTMB.rds"))
m$formula
m$time
m$extra_time



### PREDICTING 


ib_poly <- st_polygon(list(inner = rbind(innerbox,innerbox[1,])))
st_crs(ib_poly)

grid_sf <- index_grid %>% 
  st_as_sf(., coords = c("utm_x_1000", "utm_y_1000"))

# extract grid points that are within inner boundary of mesh
prop_grid <- st_intersection(grid_sf, ib_poly) %>% sf::st_drop_geometry()

saveRDS(prop_grid, file = here("data", "gsi-prop-grid.rds"))
# check
# ggplot() +
#   geom_sf(data = st_intersection(bridge_sf, ipes_sf_poly))

prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))

# extract events within ipes survey grid

grid_months <- prop_grid %>%
  replicate_df(., "region", levels(dat$region)) %>% # c("WCVI", "North/Central BC", "Columbia")) %>%
  replicate_df(., "month", c(3:11)) %>%
  mutate(month = as.numeric(month), 
         region = as.factor(region),
         yday = 160,
         utm_x_1000 = X/1000,
         utm_y_1000 = Y/1000,
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL) %>%
  filter(!month_adj %in% c(11)) #

grid_region <- grid_months %>% group_split(group_id = region) 

#newdata <- dat %>% filter(region == "WCVI")
#region <- "WCVI"

pred_ic <- array(NA, dim = c(nrow(grid_months)/nlevels(grid_months$region), nlevels(grid_months$region)),
                               dimnames = list(NULL, levels(grid_months$region)))
# se_pred_ic <- pred_ic
dim(pred_ic)

predict_list <- map(grid_region, function (newdata) {
  
  # max_chunk_length <- 5000 
  # smaller_chunks <- group_split(newdata, group_id = row_number() %/% max_chunk_length)
    
  # pred_list <- map(smaller_chunks,
  #                  ~predict(m, newdata = .x, se_fit = TRUE, re_form_iid = NA, 
  #                          offset = rep.int(-6.287114, nrow(.x))))
  
  # pred_list <- map(smaller_chunks,
  #                  ~predict(m, newdata = .x, se_fit = FALSE, re_form_iid = NA, 
  #                           offset = rep.int(-6.287114, nrow(.x))))
  
  pred <- predict(m, newdata = newdata, se_fit = FALSE, re_form_iid = NA,
          offset = rep.int(-6.287114, nrow(newdata)))
  
  #pred <- bind_rows(pred_list)
  group_id <- as.character(unique(pred$region))
  print(group_id)

  # print(head(pred))
  # pred_linkinv <- pred$est %>% m$family$linkinv()
  # print(head(pred_linkinv))

  #pred_ic[, group_id] 
  pred$pred_i <- pred$est %>% m$family$linkinv() #pred_linkinv
  pred
  #se_pred_ic[, group_id] <- pred$est_se
})

#pred <- bind_rows(predict_list)

for (i in seq_along(levels(grid_months$region))) {
  cur_region <- levels(grid_months$region)[i]
  pred_ic[,cur_region] <- predict_list[[i]]$pred_i
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


for (i in seq_along(levels(grid_months$region))) {
  #browser()
  cur_region <- levels(grid_months$region)[i]
   predict_list[[i]]$prob_i <- prob_ic[, cur_region]
}

pred <- bind_rows(predict_list)


head(pred)
head(pred)
head(pred_ic)
head(rowsum_pred_ic)
head(prob_ic)

head(pred$prob_i)
length(prob_i)

pred %>%
  ggplot(data = .) +
  geom_histogram(aes(prob_i)) +
  facet_grid(region ~ month) 



glimpse(predict_list)
glimpse(prob_i)

hist(pred$pred_i)
hist(pred$prob_i)

# pr1 <- predict(m, newdata = grid_months, 
#           offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 

saveRDS(pred, file = "data/fits/gsi-prediction-normalized-svc-int.rds")
saveRDS(predict_list, file = "data/fits/gsi-prediction-list.rds")

# pr1 <- readRDS("data/fits/gsi-prediction.rds")

pred <- readRDS("data/fits/gsi-prediction-normalized.rds")


source("R/plot_map.R")

p <- pred %>%
  dplyr::select(-X,-Y,-lat, -lon, -group_id, -est_non_rf, -est_rf, -epsilon_st) %>%
  rename(salmon_region = region) %>%
  ggplot(data = .) + # (dat, aes(X, Y, color = {{ column }})) +
  geom_raster( aes(x=utm_x_1000, y=utm_y_1000, fill = prob_i)) +
  scale_fill_viridis_c(name = "Probability of\nGSI assignment",#"Predicted\nnumber of\njuveniles\nin GSI samples",
                       labels = scales::comma,#) +#,
                       trans = fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  geom_sf(data = all_coast_km, color = "gray80", fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(pred$utm_x_1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(pred$utm_y_1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinoook GSI")


p

ggsave(p, filename = here::here("figs", "chinook_by_gsi_prob_i_svc_int.png"), width = 14, height = 12)

# +
#   ggtitle(paste0("Predicted distribution of juvenile ", species," by month (denser mesh)")) +
#   #scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
#   #scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
#   theme(legend.key.height = unit(0.06, 'npc'))

