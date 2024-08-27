# Cross-validaiton of stock compostion models

library(tidyverse)
library(sdmTMB)
library(RTMB)
library(ggsidekick)
library(here)

# loading sdmTMB model and data from model object
m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_sdmTMB.rds"))
dat <- m_svc$data
dat$offset <- m_svc$offset

sdmTMB_mesh <- m_svc$spde
sdmTMB_mesh2 <- readRDS("data/gsi-prop-mesh.rds")
unique(dat$year_adj)

# mapping of b_smooth
b_smooth_map <-  1:36
b_smooth_map[33:36] <- NA

# random folds by row
set.seed(42)
folds <-  sample(rep(seq(1L, 10), nrow(dat)), size = nrow(dat))                 

# assinging random folds to year
#folds <- sample(1:5, size = length(unique(dat$year_adj)), replace = TRUE)                 
#folds <-  sample(rep(seq(1L, 10), length(unique(dat$year_adj))), size = length(unique(dat$year_adj)))                 
#folds_key <- data.frame(year_adj = unique(dat$year_adj), fold = year_folds)

# joining folds to data
#dat <- left_join(dat, folds_key, by = "year_adj")
dat$fold <- folds

# subsetting data in samples and data out of samples
dat_in <- lapply(1:10, function(x) filter(dat, !fold %in% x))
dat_out <- lapply(1:10, function(x) filter(dat, fold %in% x))

cv_samples <- list("in" = dat_in, "out" = dat_out)
glimpse(cv_samples)

saveRDS(cv_samples, here::here("data", "stock-comp-cv-samples-seed42.R"))

lapply(dat_in, nrow)
lapply(dat_in, function (x) length(unique(x$region)))

# fitting sdmTMB models to resampled data with fold excluded
m_cv <- map(cv_samples[["in"]], function (x)  {
  new_mesh <- make_mesh(x, sdmTMB_mesh$xy_cols, mesh = sdmTMB_mesh$mesh)

  sdmTMB(
    m_svc$formula[[1]],
    family = tweedie(),
    spatial_varying = m_svc$spatial_varying_formula,
    offset = x$offset,
    data = x,
    time = "month_adj",
    extra_time = 1:12,
    spatial = "off",
    spatiotemporal = "rw",
    anisotropy = FALSE,
    mesh = new_mesh,
    silent = FALSE
  )
})

lapply(m_cv, sanity)
lapply(m_cv, function (x) sanity(x)$all_ok)

saveRDS(m_cv, here::here("data", "fits", "stock_composition_crossval_fits_sdmTMB_seed42.R"))

### Fitting RTMB models

# loading RTMB model function
nll <- readRDS(here::here("data", "stock-comp-nll-RTMB.rds"))

# functions to create TMB lists of data and pars
make_tmb_data <- function (dat, sdmTMB_mesh, m_svc) {
  list(
    observed = dat$stock_prop,
    covariate = dat$region,
    mesh = sdmTMB_mesh$mesh,
    offset = dat$effort,
    spde_m =  sdmTMB:::get_spde_matrices(sdmTMB_mesh), # change names for Q_spde functions # sdmTMB_mesh$spde,
    spde_aniso = sdmTMB:::make_anisotropy_spde(sdmTMB_mesh),
    spde =  sdmTMB_mesh$spde, # change names for Q_spde functions # sdmTMB_mesh$spde,
    year_f = as.factor(dat$year_adj_f),
    interpolator_data = fmesher::fm_basis(sdmTMB_mesh$mesh, loc = as.matrix(dat[, c("utm_x_1000", "utm_y_1000")])), 
    timesteps = sort(unique(dat$month_adj)),
    timesteps_et = 1:12,
    n_t = length(unique(dat$month_adj)), 
    n_tet = length(unique(dat$month_adj)) + 2, # +1 for extra_time
    t_i = as.numeric(dat$month_adj),
    nobs_RE = length(unique(dat$year_adj_f)),
    RE_indexes = as.numeric(as.factor(dat$year_adj_f)),
    mm_fixed = model.matrix(~ 0 + region, dat), # model matrix of fixed effects
    Zs         = m_svc$smoothers$Zs, # optional smoother basis function matrices
    Xs         = m_svc$smoothers$Xs, # optional smoother linear effect matrix
    b_smooth_start = m_svc$smoothers$b_smooth_start,
    z_i = m_svc$tmb_data$z_i, # model.matrix(~ 0 + region, dat)
    n_z = ncol(m_svc$tmb_data$z_i),
    A_spatial_index = sdmTMB_mesh$sdm_spatial_id,
    n_c = length(unique(dat$region)), # number of categories
    c_i = as.integer(dat$region)
  )
}

# this one use the smoothers instead of the sdmTMB fits
make_tmb_data_sm <- function (dat, sdmTMB_mesh, sm) {
  list(
    observed = dat$stock_prop,
    covariate = dat$region,
    mesh = sdmTMB_mesh$mesh,
    offset = dat$effort,
    spde_m =  sdmTMB:::get_spde_matrices(sdmTMB_mesh), # change names for Q_spde functions # sdmTMB_mesh$spde,
    spde_aniso = sdmTMB:::make_anisotropy_spde(sdmTMB_mesh),
    spde =  sdmTMB_mesh$spde, # change names for Q_spde functions # sdmTMB_mesh$spde,
    year_f = as.factor(dat$year_adj_f),
    interpolator_data = fmesher::fm_basis(sdmTMB_mesh$mesh, loc = as.matrix(dat[, c("utm_x_1000", "utm_y_1000")])), 
    timesteps = sort(unique(dat$month_adj)),
    timesteps_et = 1:12,
    n_t = length(unique(dat$month_adj)), 
    n_tet = length(unique(dat$month_adj)) + 2, # +1 for extra_time
    t_i = as.numeric(dat$month_adj),
    nobs_RE = length(unique(dat$year_adj_f)),
    RE_indexes = as.numeric(as.factor(dat$year_adj_f)),
    mm_fixed = model.matrix(~ 0 + region, dat), # model matrix of fixed effects
    Zs         = sm$Zs, # optional smoother basis function matrices
    Xs         = sm$Xs, # optional smoother linear effect matrix
    b_smooth_start = sm$b_smooth_start,
    z_i = m_svc$tmb_data$z_i, # model.matrix(~ 0 + region, dat)
    n_z = ncol(m_svc$tmb_data$z_i),
    A_spatial_index = sdmTMB_mesh$sdm_spatial_id,
    n_c = length(unique(dat$region)), # number of categories
    c_i = as.integer(dat$region)
  )
}

# Making parameter list
make_tmb_pars <- function (tmb_data, new_mesh, m_svc) {
  list(
    b_j = numeric(ncol(tmb_data$mm_fixed)),
    bs = if (m_svc$smoothers$has_smooths) numeric(ncol(m_svc$smoothers$Xs)) else numeric(0), # smoother linear effects
    log_tau_Z = numeric(tmb_data$n_z),
    log_tau_U = 0,
    log_kappa = 0,
    thetaf = 0,
    log_phi = 0,
    log_tau_G =  0,
    log_smooth_sigma = if (m_svc$smoothers$has_smooths) numeric(length(m_svc$smoothers$sm_dims)) else numeric(0), #
    # Random effect parameters, not returned
    RE = numeric(tmb_data$nobs_RE), # formula random effects
    b_smooth = if (m_svc$smoothers$has_smooths) numeric(sum(m_svc$smoothers$sm_dims)) else array(0), 
    rf = array(0, dim = c( new_mesh$mesh$n, tmb_data$n_tet, tmb_data$n_c)), # upsilon_stc
    rf2 = matrix(0, nrow = nrow(new_mesh$mesh$loc), ncol = tmb_data$n_z) # zeta_s
  )
}

# this one use the smoothers instead of the sdmTMB fits
make_tmb_pars_sm <- function (tmb_data, new_mesh, sm) {
  list(
    b_j = numeric(ncol(tmb_data$mm_fixed)),
    bs = if (sm$has_smooths) numeric(ncol(sm$Xs)) else numeric(0), # smoother linear effects
    log_tau_Z = numeric(tmb_data$n_z),
    log_tau_U = 0,
    log_kappa = 0,
    thetaf = 0,
    log_phi = 0,
    log_tau_G =  0,
    log_smooth_sigma = if (sm$has_smooths) numeric(length(sm$sm_dims)) else numeric(0), 
    RE = numeric(tmb_data$nobs_RE), # formula random effects
    b_smooth = if (sm$has_smooths) numeric(sum(sm$sm_dims)) else array(0), #  P-spline smooth parameters
    rf = array(0, dim = c( new_mesh$mesh$n, tmb_data$n_tet, tmb_data$n_c)), # upsilon_stc
    rf2 = matrix(0, nrow = nrow(new_mesh$mesh$loc), ncol = tmb_data$n_z) # zeta_s
  )
}

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

# Fitting RTMB models to kept data
obj_cv <- map2(cv_samples[["in"]], m_cv, function (x, m_cv)  {
  new_mesh <- make_mesh(x, sdmTMB_mesh$xy_cols, mesh = sdmTMB_mesh$mesh)
  
  # No need to recalculate smoothers as they have been recalculated in the sdmTMB
  # model fits and they are being fit to the exact same data
  tmb_data <<- make_tmb_data(x, new_mesh, m_cv) 
  params <- make_tmb_pars(tmb_data, new_mesh, m_cv) 
  
  obj <- MakeADFun(nll, params, random = c("RE", "rf", "rf2", "b_smooth"), 
                   map = list(log_tau_Z = factor(c(rep(1, 9), rep(2, 9))),
                              b_smooth = factor(b_smooth_map),
                              log_smooth_sigma = factor(c(1,2,3,4,5,6,7,8,NA))),
                   silent = FALSE)
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  obj
})
saveRDS(obj_cv, here::here("data", "fits", "stock_composition_crossval_fits_RTMB_seed.R"))

# calculating sdreport
sdrep_cv <- lapply(obj_cv, sdreport)
saveRDS(sdrep_cv, here::here("data", "fits", "stock_composition_crossval_sdrep_RTMB_seed.R"))

#saveRDS(sdrep_cv, here::here("data", "fits", "stock_composition_crossval_sdrep_RTMB_seed42.R"))

# Loading saved resamples and fitted models
obj_cv <- readRDS(here::here("data", "fits", "stock_composition_crossval_fits_RTMB_seed42.R"))

m_cv <- readRDS(here::here("data", "fits", "stock_composition_crossval.R"))
cv_samples <- readRDS(here::here("data", "stock-comp-cv-samples.R"))

# Predictions on withheld data for sdmTMB
pred_cv <- map2(m_cv, cv_samples[["out"]], function (m, d_out) {
  predict(m, newdata = d_out, 
            offset = d_out$offset, re_form_iid = NA)$est %>% m$family$linkinv()
})

glimpse(pred_cv)

# Predictions on withheld data for RTMB
pred_cv_RTMB <- map2(obj_cv, cv_samples[["out"]], function (obj, d_out) {
  sm_pred <- sdmTMB:::parse_smoothers(stock_prop ~ 0 + region + (1 | year_adj_f) + 
                                        s(scale_month_adj, by = region, bs = "tp", k = 6), 
                                      data = dat,
                                      newdata = d_out, 
                                      basis_prev = m_svc$smoothers$basis_out)
  
  new_mesh <- make_mesh(d_out, sdmTMB_mesh$xy_cols, mesh = sdmTMB_mesh$mesh)
  tmb_data <<- make_tmb_data_sm(d_out, new_mesh, sm_pred) # m_cv is on the kept data not withheld
  #params <- make_tmb_pars_sm(tmb_data, new_mesh, sm_pred) 
  lp <- obj$env$last.par.best
  obj$report(par = lp)$mu
})

glimpse(pred_cv_RTMB)

# calculating log likelihood based on dtweedie() for sdmTMB
cv_loglik_sdmTMB <- pmap(list(m_cv, cv_samples[["out"]], pred_cv), 
                  function (m_object, withheld_dat, withheld_mu) {
  sdmTMB:::ll_tweedie(m_object, withheld_dat$stock_prop, withheld_mu)
}) 

# calculating log likelihood based on dtweedie() for RTMB
cv_loglik_RTMB <- pmap(list(sdrep_cv, cv_samples[["out"]], pred_cv_RTMB), 
                       function (sdrep, withheld_dat, withheld_mu) {
                         tweedie_p <- sdrep$value["tweedie_p"]
                         phi <- sdrep$value["phi"]
                         #browser()
                         fishMod::dTweedie(y = withheld_dat$stock_prop, mu = withheld_mu, 
                                           p = tweedie_p, phi = phi, LOG = TRUE)
                       }) 

lapply(cv_loglik_sdmTMB, sum) %>% unlist() %>% sum()
# [1] -11803.91
lapply(cv_loglik_RTMB, sum) %>% unlist() %>% sum()
# [1] -11454.04
# RTMB model has a slightly higher nll than sdmTMB one hence it is slightly better

tibble(sdmTMB = lapply(cv_loglik_sdmTMB, sum) %>% unlist(),
           RTMB = lapply(cv_loglik_RTMB, sum) %>% unlist()) %>%
  mutate(diff = sdmTMB - RTMB)
# in 6 out of 10 folds the sdmTMB model has a higher nll (i.e. is better)

# in 6 out of 10 folds the sdmTMB model has a higher nll (i.e. is better) in seed 42 as well