# Diagnostic plots

library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggsidekick)
library(sdmTMB)
library(here)
library(sdmTMBextra)

source(here("R", "plot_map.R"))
mdmesh <- readRDS(here("data", "fits", "dmesh_allcoast.rds")) 

# mout3 <- readRDS(here("data", "fits", "mout3_allcoast.rds")) 
# dat_tbl_mout3 <- readRDS(here("data", "fits", "all_mout3_allcoast.rds"))
#dat_tbl_dmesh <- readRDS(here("data", "fits", "all_dmesh_allcoast.rds"))

index_grid <- readRDS(here("data", "pred_grid_all_coast.rds")) %>%
  mutate(year = 2019L, 
         year_f = as.factor(year),
         survey_f = as.factor("hss"),
         target_depth = 0,
         day_night = "DAY",
         scale_depth = -0.5427145)

index_grid_jun <- mutate(index_grid, month = 6,
                         month_f = month(month, label = TRUE),
                         month_adj = ifelse(month > 2, month - 2, month + 10)) # adjusting to start in March)

glimpse(index_grid_jun)
table(index_grid_jun$year)

grid_months <- index_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
  NULL)

dim(grid_months)


# grid_months_years <- grid_months %>%
#   select(-year, -year_f) %>%
#   replicate_df(., "year", c("1998", "2008", "2018")) %>%
#   mutate(year_f = as.factor(year), month_f = month(month, label = TRUE))

# Same
# grid_months_years2 <- index_grid %>%
#   replicate_df(., "month", 1:12) %>%
#   replicate_df(., "year", c("1998", "2008", "2018")) %>%
#   mutate(year_f = as.factor(year), month_f = month(month, label = TRUE))

dim(index_grid)
dim(index_grid_jun)
dim(grid_months)
#dim(grid_months_years)

# # Chinook only predictions
# pred_jun <- predict(mout3, newdata = index_grid_jun, 
#                     offset = rep.int(-6.287114, nrow(index_grid_jun)), re_form_iid = NA)
# pred_jun_dmesh <- predict(mdmesh, newdata = index_grid_jun, 
#                     offset = rep.int(-6.287114, nrow(index_grid_jun)), re_form_iid = NA)
# glimpse(pred_jun)
# 
# pred_months <- predict(mout3, newdata = grid_months, 
#                        offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 
# 
# 
# 
# plot_map(pred_jun, exp(est), area = "BC") +
#   scale_color_viridis_c(trans = fourth_root_power_trans()) +
#   ggtitle("Predicted juvenile chinook distribution in June")
# 
# png(here::here("figs", "chinook_juv_jun.png"),
#     height = 8, width = 8, units = "in", res = 200)
# plot_map(pred_jun_dmesh, exp(est), area = "BC") +
#   scale_color_viridis_c(trans = fourth_root_power_trans()) +
#   ggtitle("Predicted juvenile chinook distribution in June (dense mesh)")
# dev.off()

# Monthly predictions for all species
all_pred_months <- purrr::map(
  dat_tbl_dmesh$fit,
  ~ predict(.x, newdata = grid_months, 
            offset = rep.int(-6.287114, nrow(grid_months)), re_form_iid = NA) 
)

# saving to df ############# and adding adjusted month factor levels so maps start from May (except chum)
dat_tbl_dmesh$pred_grid <- all_pred_months

saveRDS(dat_tbl_dmesh, here::here("data", "fits", "all_dmesh_allcoast_withpred.rds"))

dat_tbl_dmesh <- readRDS(here::here("data", "fits", "all_dmesh_allcoast_withpred.rds"))

# png(here::here("figs", "seb", "chinook_juv_months_wdata4.png"),
#     height = 7, width = 10, units = "in", res = 300)
# plot_map(pred_months, exp(est), show_raw_data = TRUE, chinook_dat) + 
#   ggtitle("Predicted distribution of juvenile chinook by month") +
#   scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
#   scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
#   theme(legend.key.height = unit(0.06, 'npc'))
# dev.off()


### looping through to obtain distribution maps for all species

for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
  p <- plot_map(dat_tbl_dmesh$pred_grid[[i]], exp(est), area = "BC",
                show_raw_data = FALSE, dat_tbl_dmesh$data[[i]]) + 
    annotation_scale(location = "bl", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month (denser mesh)")) +
    #scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
    #scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  png(here::here("figs", paste0("BC_atlas_noraw_dmesh_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}


############

sims_list <- furrr::future_map(
  dat_tbl_dmesh$fit, function (x) {
    object <- x
    samp <- sample_mle_mcmc(object, mcmc_iter = 22L, mcmc_warmup = 20L, mcmc_chains = 5L,
                            stan_args = list(thin = 5L, cores = 5L))
    obj <- object$tmb_obj
    random <- unique(names(obj$env$par[obj$env$random]))
    pl <- as.list(object$sd_report, "Estimate")
    fixed <- !(names(pl) %in% random)
    map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
    obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
    obj_mle <- object
    obj_mle$tmb_obj <- obj
    obj_mle$tmb_map <- map
    simulate(obj_mle, mcmc_samples = sdmTMBextra::extract_mcmc(samp), nsim = 200L)
  }
)

saveRDS(sims_list,
        here::here("data", "fits", "nb_mcmc_draws_mout3.rds"))
}



sims_list <- readRDS(here::here("data", "fits", "nb_mcmc_draws_mout3.rds"))





all_fit_tbl$sims <- sims_list



### Plot distribution of sims vs observed:

for (i in 1:nrow(dat_tbl_mout3)) {
  species <- dat_tbl_mout3$species[i]
  sims <- all_fit_tbl$sims[[i]]
  dim(sims)
  
  obs <- all_fit_tbl$data[[i]]
  
  print(species)
  
  p <- ggplot(data = ,
              aes(+
    theme(legend.key.height = unit(0.06, 'npc'))
  png(here::here("figs", paste0(species, "_juv_months_loop.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}















######



dharma_res <- DHARMa::createDHARMa(
  simulatedResponse = dat_tbl_mout3$sims[[1]][1:(nrow(dat_tbl_mout3$sims[[1]])-2),],
  observedResponse = dat_tbl_mout3$data[[1]]$n_juv,
  fittedPredictedResponse = fixed_pred_list[[1]]
)
plot(dharma_res)
