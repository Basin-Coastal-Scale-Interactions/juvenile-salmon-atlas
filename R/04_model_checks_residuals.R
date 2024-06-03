### Model and residual checks
# Can't use normal residuals 

library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggsidekick)
library(sdmTMB)
library(here)
library(sdmTMBextra)
library(rstan)
library(TMB)

source(here("R", "plot_map.R"))
#mdmesh <- readRDS(here("data", "fits", "dmesh_allcoast.rds")) 

dat_tbl_dmesh <- readRDS(here::here("data", "fits", "all_dmesh_allcoast.rds"))
fits_list_mdmesh <- readRDS(here::here("data", "fits", "fits_list_mdmesh.rds"))

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

### TRYING TO USE MCMC RESIDUALS BUT SEEMS PROHIBITIVELY SLOW
### 
# mcmcs_samples <- sample_mle_mcmc(
#   dat_tbl_dmesh$fit[[2]],
#   mcmc_iter = 5,
#   mcmc_warmup = 2,
# )
# asd <- extract_mcmc(mcmcs_samples)
# dim(asd)
# head(asd)
# pred <- predict(obj_mle, mcmc_samples = extract_mcmc(samp), 
#                 model = model[[1]], nsim = nsim, offset = obj_mle$tmb_data$offset_i)
# pred
# 
# mcmcs_samples <- predict_mle_mcmc(
#   dat_tbl_dmesh$fit[[2]],
#   mcmc_iter = 22L,
#   mcmc_warmup = 20L,
#   print_stan_model = FALSE,
#   stan_args = NULL,
#   nsim = 1
# )

# Calculating residuals using MLE-MVN approach
resids <- purrr::map(dat_tbl_dmesh$fit, function (x) as.vector(residuals(x, type = "mle-mvn")))
names(resids) <- dat_tbl_dmesh$species

lapply(resids, function (x) is.infinite(x) |> which())

png(here::here("figs", "diagnostics", "mle_mlv_residuals.png"), height = 9, width = 6, units = "in", res = 200)
par(mfrow = c(3,2))
map2(resids, dat_tbl_dmesh$species, function (x, sp) {
  res <- x[is.finite(x)] # some Inf
  print(sp)
  qqnorm(res, main = paste0("Q-Q plot - mle-mvn residuals - ",sp)); qqline(res)
  })
dev.off()
par(mfrow = c(1,1))

# plotting residuals in space
purrr::pmap(
  list(resids, dat_tbl_dmesh$data, dat_tbl_dmesh$species),
  function(res, dat, sp) {
    df <- dat
    df$resids <- res
    
    ggplot(df) +
      geom_point(aes(x = utm_x, y = utm_y, colour = resids), shape = 1, size = 1.5) +
      facet_wrap(~month) +
      scale_colour_gradient2(low = "red", high = "blue") +
      ggsidekick::theme_sleek() +
      coord_fixed() +
      # geom_sf(data = all_coast, color = "gray80", fill = NA) +
      ggtitle(paste0("Randomized quantile residuals - ", sp))
    ggsave(here::here("figs", "diagnostics", paste0("spatial_residuals_", sp, ".png")), 
           width = 9, height = 7)
  }
)

# Looking at DHARMa residuals
dh_resids <- map(fits_list_mdmesh, function (x) {
  sim <- simulate(x, nsim = 100, type = "mle-mvn")
  sdmTMB::dharma_residuals(sim, x, return_DHARMa = TRUE)
})

# plotting DHARMa residuals
map2(dh_resids, dat_tbl_dmesh$species, function (x, sp) {
  png(here::here("figs", "diagnostics", paste0("dharma_resids_", sp, ".png")),
      height = 5, width = 10, units = "in", res = 300)
  plot(x, title = paste0("DHARMa residual - ", sp))
  dev.off()
})

saveRDS(dh_resids, file = here::here("data", "dharma-objects.rds"))
#dh_resids <- readRDS(file = here::here("data", "dharma-objects.rds"))

# DHARMa::testResiduals(dh_resids[[1]])

# DHARMa Moran's I test for distance-based autocorrelation
# need to remove duplicate coordinates
pmap(list(dh_resids, dat_tbl_dmesh$data, dat_tbl_dmesh$species), function (res, dat, sp) {
  dat1 <- dat %>%
    group_by(utm_x, utm_y) %>% 
    mutate(notdupe = n()==1)

  res2 <-  DHARMa::recalculateResiduals(res, sel = dat1$notdupe,
                                          aggregateBy = sum)
  dat1nd <- filter(dat1, notdupe)
  
  print(sp)
  png(here::here("figs", "diagnostics", paste0("dharma_distance_autocorr_", sp, ".png")),  
      height = 8, width = 8, units = "in", res = 300)
  print(DHARMa::testSpatialAutocorrelation(res2, 
                                   x = dat1nd$utm_x, 
                                   y = dat1nd$utm_y))
  mtext(sp, 3, adj = 0)
  dev.off()
})

# Moran's I autocorrelation test for residuals
moran_list <- purrr::pmap(
  list(dat_tbl_dmesh$fit, resids, dat_tbl_dmesh$data, dat_tbl_dmesh$species),
  function(fit, res, dat, sp) {
    #browser()
    df <- dat
    df$resids <- res
    
    inv_dist <- as.matrix(dist(df[,c("utm_x","utm_y")]))
    diag(inv_dist) <- 0
    moranI <- ape::Moran.I(df$resids, weight=inv_dist, na.rm=TRUE)
    #moranI$species <- sp
    moranI <- c("species" = sp, moranI)
    moranI
  }
)

moran_list %>% bind_rows()
# Seems like values for chum are positively autocorrelated

##########################################################

### Contribution of fixed vs random effects to predictions

# fixed_pred_list <- purrr::map2(
#   dat_tbl_mout3$fit, dat_tbl_mout3$data,
#   ~ predict(.x, newdata = .y)$est_non_rf %>% 
#     .x$family$linkinv(.)
# )
# 
# dat_tbl_mout3$fixed_pred_list <- fixed_pred_list


# # Monthly predictions for all species 
# all_pred_months <- purrr::map(
#   dat_tbl_dmesh$fit,
#   ~ predict(.x, newdata = grid_months, 
#             offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 
# )

# Fixed effect only predictions on PREDICTIVE GRID
fixed_pred_grid <- purrr::map(
  dat_tbl_dmesh$fit, 
  ~ predict(.x, newdata = grid_months, 
            offset = rep.int(-6.287114, nrow(grid_months)))$est_non_rf %>% 
    .x$family$linkinv(.)
)

fixed_pred_grid_1 <- fixed_pred_grid

# Fixed effect only predictions on PREDICTIVE GRID --- IN CHUNKS (to avoid std:bad_alloc errors)
fixed_pred_grid_4 <- purrr::map(
  dat_tbl_dmesh$fit, function (fit) {
    grid_chunks <- grid_months %>% group_split(group_id = row_number() %/% 500000) # 500000 is max size of chunk
    
    predict_list <- map(grid_chunks, 
                        ~predict(fit, newdata = .x,
                                 offset = rep.int(-6.287114, nrow(.x)))$est_non_rf %>% 
                          fit$family$linkinv(.))
    unlist(predict_list)
  })



fit1 <- fits_list_mdmesh[[1]]

grid_chunks <- grid_months %>% group_split(group_id = row_number() %/% 500000)
lapply(grid_chunks, nrow) %>% unlist

out <- map(grid_chunks, ~ predict(fit1, newdata = .x, offset = rep.int(-6.287114, nrow(.x)))$est_non_rf %>% 
      fit1$family$linkinv(.))

# fixed effect only predictions ON DATA POINTS
fixed_pred_list <- purrr::map2(
  dat_tbl_dmesh$fit, dat_tbl_dmesh$data,
  ~ predict(.x, newdata = .y, offset = rep.int(-6.287114, nrow(.y)))$est_non_rf %>% 
    .x$family$linkinv(.)
)

# Random field only predictions ON PREDICTIVE GRID
rf_pred_grid <- purrr::map(
  dat_tbl_dmesh$fit, 
  ~ predict(.x, newdata = grid_months, 
            offset = rep.int(-6.287114, nrow(grid_months)))$est_rf %>% 
    .x$family$linkinv(.)
)

for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
  p <- plot_map(grid_months, fixed_pred_grid[[i]], area = "BC", show_raw_data = TRUE, dat_tbl_dmesh$data[[i]]) + 
    ggtitle(paste0("NON random field - Predicted distribution of juvenile ", species," by month")) +
    # scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
    # scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  png(here::here("figs", "re_fe_contributions", paste0("nonRF_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}

for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
  p <- plot_map(grid_months, rf_pred_grid[[i]], area = "BC", show_raw_data = TRUE, dat_tbl_dmesh$data[[i]]) + 
    ggtitle(paste0("Random field - predicted distribution of juvenile ", species," by month")) +
    # scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
    # scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  png(here::here("figs", "re_fe_contributions", paste0("RF_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}




