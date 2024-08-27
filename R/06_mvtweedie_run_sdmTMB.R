
library(tidyverse)
#library(RTMB)
library(sdmTMB)
library(sdmTMBextra)
library(mgcv)
#library(mvtweedie)
library(here)
library(sf)
#  Load GSI data ------------------------------------------------------------

dat <- readRDS(here::here("data", "chinook_gsi_counts_20240725.rds")) %>%
  mutate(season_n = as.numeric(season_f),
         scale_month_adj = scale(month_adj)[, 1],
         year_adj = ifelse(month > 2, year, year - 1), # adjusting year to represent fish cohorts
         # year_f = as.factor(year),
         year_adj_f = as.factor(year_adj))

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

# observations
ggplot(dat, aes(utm_x, utm_y, colour = stock_prop)) +
  geom_point() +
  facet_wrap(~region)
# Mesh and SPDE matrix construction -----------------------------------------

dat_coords <- dat %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()

inla_mesh_raw <- fmesher::fm_mesh_2d_inla(
  loc = dat_coords,
  max.edge = 500,
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
  stock_prop ~ 0 + region + (1 | year_adj_f) + s(scale_month_adj, by = region, bs = "tp", k = 6),
  family = tweedie(),
  spatial_varying = ~ 0 + region * scale_month_adj,
  offset = dat$effort,
  data = dat,
  time = "month_adj",
  extra_time = c(2,11),
  spatial = "off",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  mesh = sdmTMB_mesh,
  silent = FALSE,
  control = sdmTMBcontrol(newton_loops = 0L, nlminb_loops = 1L),
  do_fit = FALSE
)

b_smooth_map <-  1:36
b_smooth_map[33:36] <- NA
b_smooth_map

m_svc <- update(m_svc, do_fit = TRUE,
                control = sdmTMBcontrol(map = list(ln_tau_Z = factor(c(rep(1, 10), rep(2, 10)))
                                                   #b_smooth = factor(b_smooth_map),
                                                   #ln_smooth_sigma = factor(c(1,2,3,4,5,6,7,8,NA)
                                                   ))
                                        )
                )


sanity(m_svc)

m_svc
m_svc$sd_report
tidy(m_svc, "ran_pars", conf.int = TRUE)

saveRDS(m_svc, file = here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))


m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))
m_svc$formula
m_svc$time
m_svc$extra_time

m_svc$data$region %>% table


### PREDICTING 

# Whole coast index grid that will get cut
index_grid_old <- readRDS(here("data", "index_grid.rds")) 
glimpse(index_grid_old)

index_grid <- readRDS(here("data", "pred_grid_all_coast.rds")) %>%
  # filter(row_number() %% 5 == 1) %>% ### DOWNSAMPLING GRID %>%
  # arrange(X,Y) %>%
  # filter(row_number() %% 5 == 1) %>% ### DOWNSAMPLING GRID %>%
  mutate(year = 2019L, 
         year_f = as.factor(year),
         survey_f = as.factor("hss"),
         target_depth = 0,
         day_night = "DAY",
         scale_depth = -0.5427145) %>%
  rename(year_adj_f = year_f, year_adj = year)
glimpse(index_grid)

# Polygon of points from inner sdmTMB mesh
ib_poly <- st_polygon(list(inner = rbind(innerbox,innerbox[1,])))
st_crs(ib_poly)

grid_sf <- index_grid %>% 
  st_as_sf(., coords = c("utm_x_1000", "utm_y_1000"), remove = FALSE)

# extract grid points that are within inner boundary of mesh
prop_grid <- st_intersection(grid_sf, ib_poly) %>% sf::st_drop_geometry()
glimpse(prop_grid)

saveRDS(prop_grid, file = here("data", "gsi-prop-grid.rds"))
# check
# ggplot() +
#   geom_sf(data = st_intersection(bridge_sf, ipes_sf_poly))

prop_grid <- readRDS(here("data", "gsi-prop-grid.rds"))


# functions for quickly applying and unapplying scale()
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}


grid_months <- prop_grid %>%
  replicate_df(., "region", levels(dat$region)) %>% # c("WCVI", "North/Central BC", "Columbia")) %>%
  replicate_df(., "month", c(3:11)) %>%
  mutate(month = as.numeric(month), 
         region = as.factor(region),
         #yday = 160,
         utm_x_1000 = X/1000,
         utm_y_1000 = Y/1000,
         year_adj_f = as.factor(2019),
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         scale_month_adj = scale_est(month_adj, dat$month_adj),
         NULL) %>%
  filter(!month_adj %in% c(11)) #

glimpse(grid_months)


grid_region <- grid_months %>% group_split(group_id = region) 
glimpse(grid_region)
#newdata <- dat %>% filter(region == "WCVI")
#region <- "WCVI"



predict_list <- map(grid_region, function (newdata) {
  
  # max_chunk_length <- 5000 
  # smaller_chunks <- group_split(newdata, group_id = row_number() %/% max_chunk_length)
    
  # pred_list <- map(smaller_chunks,
  #                  ~predict(m, newdata = .x, se_fit = TRUE, re_form_iid = NA, 
  #                          offset = rep.int(median(dat$effort), nrow(.x))))
  
  # pred_list <- map(smaller_chunks,
  #                  ~predict(m, newdata = .x, se_fit = FALSE, re_form_iid = NA, 
  #                           offset = rep.int(median(dat$effort), nrow(.x))))
  
  pred <- predict(m_svc, newdata = newdata, se_fit = FALSE, re_form_iid = NA,
          offset = rep.int(median(dat$effort), nrow(newdata)))
  
  #pred <- bind_rows(pred_list)
  group_id <- as.character(unique(pred$region))
  print(group_id)

  # print(head(pred))
  # pred_linkinv <- pred$est %>% m$family$linkinv()
  # print(head(pred_linkinv))

  #pred_ic[, group_id] 
  pred$pred_i <- pred$est %>% m_svc$family$linkinv() #pred_linkinv
  pred
  #se_pred_ic[, group_id] <- pred$est_se
})

#pred <- bind_rows(predict_list)

pred_ic <- array(NA, dim = c(nrow(grid_months)/nlevels(grid_months$region), 
                             nlevels(grid_months$region)),
                 dimnames = list(NULL, levels(grid_months$region)))
# se_pred_ic <- pred_ic
dim(pred_ic)

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
length(pred$prob_i)

pred %>%
  ggplot(data = .) +
  geom_histogram(aes(prob_i)) +
  facet_grid(region ~ month) 

glimpse(predict_list)
#glimpse(prob_i)

hist(pred$pred_i)
hist(pred$prob_i)

# pr1 <- predict(m, newdata = grid_months, 
#           offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 

saveRDS(pred, file = "data/fits/gsi-prediction-normalized-svc-sdmTMB.rds")
saveRDS(predict_list, file = "data/fits/gsi-prediction-normalized-svc-sdmTMB-list.rds")

# pr1 <- readRDS("data/fits/gsi-prediction.rds")

pred <- readRDS("data/fits/gsi-prediction-normalized-svc-sdmTMB.rds")


source("R/plot_map.R")

p <- pred %>%
  dplyr::select(utm_x_1000, utm_y_1000, prob_i, salmon_region = region, month_f) %>%
  ggplot(data = .) + # (dat, aes(X, Y, color = {{ column }})) +
  geom_raster( aes(x=utm_x_1000, y=utm_y_1000, fill = prob_i)) +
  scale_fill_viridis_c(name = "Probability of\nGSI assignment",#"Predicted\nnumber of\njuveniles\nin GSI samples",
                       labels = scales::comma,#) +#,
                       trans = ggsidekick::fourth_root_power_trans()) +
  ggsidekick::theme_sleek() +
  geom_sf(data = all_coast_km, color = "gray80", fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(pred$utm_x_1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(pred$utm_y_1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinoook GSI")

p

ggsave(p, filename = here::here("figs", "chinook_by_gsi_prob_i_svc_mapped_tauZ_sdmTMB_20240725.png"), width = 14, height = 12)
#ggsave(p, filename = here::here("figs", "chinook_by_gsi_prob_i_svc_int_no_tauZ_map.png"), width = 14, height = 12)

# +
#   ggtitle(paste0("Predicted distribution of juvenile ", species," by month (denser mesh)")) +
#   #scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
#   #scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
#   theme(legend.key.height = unit(0.06, 'npc'))

