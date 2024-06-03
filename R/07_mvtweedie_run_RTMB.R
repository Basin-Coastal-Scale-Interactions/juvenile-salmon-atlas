# RTMB 

library(tidyverse)
library(RTMB)
library(sdmTMB)
#library(Matrix)
# Load data & functions -----------------------------------------------------

dat <- readRDS(here::here("data", "chinook_gsi_counts.rds")) %>%
  mutate(season_n = as.numeric(season_f),
  scale_month_adj = scale(month_adj)[, 1],
  year_adj = ifelse(month > 2, year, year - 1), # adjusting year to represent fish cohorts
 # year_f = as.factor(year),
  year_adj_f = as.factor(year_adj))

Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}

# Mesh and SPDE matrix construction -----------------------------------------

sdmTMB_mesh <- readRDS("data/gsi-prop-mesh.rds")
spde <- sdmTMB_mesh$spde
mesh <- sdmTMB_mesh$mesh

# compute bilinear interpolation matrix from mesh to data:
# Identical to sdmTMB_mesh$A_st
interpolator_data <- fmesher::fm_basis(
  mesh,
  loc = as.matrix(dat[, c("utm_x_1000", "utm_y_1000")])
)

# extracting sdmTMB model object without fitting to simplify RTMB data entry

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

m_svc$spatial_varying_formula
m_svc$spatial_varying
m_svc$tmb_data$spatial_covariate
m_svc$tmb_data$spatial_only
m_svc$tmb_data$z_i %>% head()
table(dat$month_adj)

# RTMB model

nll <- function(par) {
  getAll(par, tmb_data)
  
  # parameter transformations
  #tau0 <- exp(log_tau0) # switch spatial = "off"
  #tau_E <- exp(log_tau_E)
  tau_U <- exp(log_tau_U)
  tau_Z <- exp(log_tau_Z)
  kappa <- exp(log_kappa)
  phi <- exp(log_phi)
  sigma_G <- exp(log_tau_G)
  smooth_sigma <- exp(log_smooth_sigma)
  tweedie_p <- plogis(thetaf) + 1
  
  # compute precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  #H <- MakeHv(log_H_input)
  #REPORT(H)
  
  # compute precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  # Two random fields, one for the intercept and one for the slope
  Q <- Q_spde(spde_m, kappa)
  #Q0 <- tau0^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  #Q_1 <- tau^2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  #Q <- tau^2 * Q_spde_aniso_v(spde_aniso, kappa, H)

  # Random walk of the spatiotemporal random field (intercept)
  # with separate walks for each group
  #rf0 %~% dgmrf(0, Q0) # switch spatial = "off"
  for (c in 1:n_c) {
    for (t in 1:n_tet) { # include extra_time timesteps in the RW
      if (t == 1) rf[, t, c] %~% dgmrf(0, tau_U^2 * Q)
      else {
        rf[, t, c] %~% dgmrf(rf[, t-1, c], tau_U^2 * Q)
      }
    }
  }
  
  # # Random walk of the spatiotemporal random field (intercept)
  # #rf0 %~% dgmrf(0, Q0) # switch spatial = "off"
  # rf[,1] %~% dgmrf(0, tau_E^2 * Q)
  # for (t in 2:n_tet) { # include extra_time timesteps in the RW
  #   rf[,t] %~% dgmrf(rf[,t-1], tau_E^2 * Q)
  # }
  
  # project random effects from vertices to data locations:
  #rf_at_observations0 <- interpolator_data %*% rf0  # switch spatial = "off"
  #rf_at_observations <- matrix(0, nrow = length(observed), ncol = length(timesteps_et)) # 
  rf_at_observations <- array(0, dim = c(length(observed), ncol = length(timesteps_et), n_c))
  rf2_at_observations <- matrix(0, nrow = length(observed), ncol = n_z)
  
  for (t in timesteps) { # `timesteps` doesn't include extra time values
    rf_at_observations[, t, ] <- as.vector(interpolator_data %*% rf[, t, ]) # rf column for extra time not used here
    #rf_at_observations[, t] <- as.vector(interpolator_data %*% rf[, t]) # rf column for extra time not used here
  }

  # projecting to data locations and drawing from separate distributions for
  # levels in covariate coefficients
  for (z in 1:n_z) {
    rf2_at_observations[, z] <- as.vector(interpolator_data %*% rf2[, z])
    rf2[, z] %~% dgmrf(0, tau_Z[z]^2 * Q)
  }
  
  log_mu <- numeric(length(observed))
  eta_iid_re_i <- numeric(length(observed))
  eta_fixed_i <- numeric(length(observed))
  eta_smooth_i <- numeric(length(observed))
  
  # Smoother effects
  for (s in 1:length(b_smooth_start)) { # iterate over # of smooth elements
    #beta_s <- b_smooth[(b_smooth_start[s]+1):(b_smooth_start[s]+ncol(Zs[[s]]))] #numeric(ncol(Zs[[s]]))
    beta_s <- numeric(ncol(Zs[[s]])) # reset beta_s for every smooth element

    for (j in 1:length(beta_s)) {
      beta_s[j] <- b_smooth[b_smooth_start[s] + j] # take out appropriate segment of b_smooth
    }
    
    # random effect part of smoother
    eta_smooth_i <- eta_smooth_i + rowSums(sweep(Zs[[s]], MARGIN = 2, beta_s, `*`)) 
    # Equivalent of Zs[[s]] %*% beta_s but this doesn't work on RTMB
    
    beta_s %~% dnorm(0, smooth_sigma[s], log = TRUE) # smoother as random effect distribution
  }
  
  # fixed effect part of smoother
  eta_smooth_i <- eta_smooth_i + rowSums(sweep(Xs, MARGIN = 2, bs, `*`)) 
  # Equivalent of Xs %*% bs but this doesn't work on RTMB

  # pick out the appropriate time slice for each row of data
  # and add on any other components:
  for (i in seq_along(observed)) {
    eta_iid_re_i[i] <- RE[RE_indexes[i]] # random effect
    eta_fixed_i[i] <- sum(mm_fixed[i,] * b_j) # fixed effects
    #eta_smooth_i[i] <- sum(Zs[i,] * beta_s) + Xs[i, 1] * bs # smoother effects, now estimated outside of this loop

    # for (int z = 0; z < n_z; z++)
    #   eta_i(i,m) += zeta_s_A(i,z,m) * z_i(i,z);
   #browser()
    log_mu[i] <- offset[i] + # rf_at_observations0[i] + # switch spatial = "off"
      eta_fixed_i[i] + eta_iid_re_i[i] + eta_smooth_i[i] + rf_at_observations[i, t_i[i], c_i[i]] +
      sum(rf2_at_observations[i, ] * z_i[i,]) # spatially varying coefficients
  }
  
  # For prediction
  mu <- exp(log_mu)
  ADREPORT(mu)
  
  #REPORT(rf_at_observations0) #< NEW
  REPORT(rf_at_observations)
  REPORT(rf2_at_observations)
  
  observed %~% dtweedie(exp(log_mu), phi, tweedie_p, log = TRUE)
  RE %~% dnorm(0, sigma_G, log = TRUE) # random effect distribution
  
  ADREPORT(tweedie_p)
  
  range <- sqrt(8) / kappa
  ADREPORT(range)
  #sigma0 <- 1 / sqrt(4 * pi * exp(2 * log_tau0 + 2 * log_kappa))
  sigma_U <- 1 / sqrt(4 * pi * exp(2 * log_tau_U + 2 * log_kappa))
  #sigma_E <- 1 / sqrt(4 * pi * exp(2 * log_tau_E + 2 * log_kappa))
  sigma_Z <- 1 / sqrt(4 * pi * exp(2 * log_tau_Z + 2 * log_kappa))
  #ADREPORT(sigma0) #< NEW
  ADREPORT(sigma_U)
  ADREPORT(sigma_Z)
  ADREPORT(phi)
  REPORT(sigma_G)
  
  # REPORT(b_smooth)  # smooth coefficients for penalized splines
  # REPORT(log_smooth_sigma)
  
  REPORT(b_smooth)         # smooth coefficients for penalized splines
  REPORT(log_smooth_sigma) # standard deviations of smooth random effects, in log-space
}

# parsing smoothers
sm <- m_svc$smoothers
# sm <- sdmTMB:::parse_smoothers(stock_prop ~  0 + region + (1 | year_f)
#                                + s(month_adj, by = region, bs = "tp", k = 6), dat)

tmb_data <- list(
  observed = dat$stock_prop,
  covariate = dat$region,
  mesh = sdmTMB_mesh$mesh,
  offset = dat$effort,
  spde_m =  sdmTMB:::get_spde_matrices(sdmTMB_mesh), # change names for Q_spde functions # sdmTMB_mesh$spde,
  spde_aniso = sdmTMB:::make_anisotropy_spde(sdmTMB_mesh),
  spde =  sdmTMB_mesh$spde, # change names for Q_spde functions # sdmTMB_mesh$spde,
  year_f = as.factor(dat$year_adj_f),
  interpolator_data = interpolator_data, # or A_st
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

glimpse(tmb_data)

# Order of parameters here is what determines order in sdrep output
par <- list(
  #log_tau0 = 0, #< NEW
  #log_H_input = numeric(2),
  b_j = numeric(ncol(tmb_data$mm_fixed)),
  bs = if (sm$has_smooths) numeric(ncol(sm$Xs)) else numeric(0), # smoother linear effects
  #zeta_s = matrix(0, nrow = nrow(sdmTMB_mesh$mesh$loc), ncol= tmb_data$n_z),
  log_tau_Z = numeric(tmb_data$n_z),
  #log_tau_E = 0,
  log_tau_U = 0,
  log_kappa = 0,
  thetaf = 0,
  log_phi = 0,
  log_tau_G =  0,
  log_smooth_sigma = if (sm$has_smooths) numeric(length(sm$sm_dims)) else numeric(0), # variances of spline REs if included
  # rf0 = rep(0, length = mesh$n), #< NEW
  # Random effect parameters, not returned
  RE = numeric(tmb_data$nobs_RE), # formula random effects
  b_smooth = if (sm$has_smooths) numeric(sum(sm$sm_dims)) else array(0), #  P-spline smooth parameters
  #rf = matrix(0, nrow = mesh$n, ncol = tmb_data$n_tet), # epsilon_st
  rf = array(0, dim = c( mesh$n, tmb_data$n_tet, tmb_data$n_c)), # upsilon_stc
  rf2 = matrix(0, nrow = nrow(sdmTMB_mesh$mesh$loc), ncol = tmb_data$n_z) # zeta_s
)
glimpse(par)


# fixing the last element of the smoother as it otherwise estimates a very high log_smooth_sigma
b_smooth_map <-  1:36
b_smooth_map[33:36] <- NA
b_smooth_map

# debugonce(nll)
obj <- MakeADFun(nll, par, random = c("RE", "rf", "rf2", "b_smooth"), 
                 map = list(log_tau_Z = factor(c(rep(1, 9), rep(2, 9))),
                            b_smooth = factor(b_smooth_map),
                            log_smooth_sigma = factor(c(1,2,3,4,5,6,7,8,NA))
                            ),
                # profile = c("b_j", "bs"),
                 silent = FALSE)  #< NEW: rf0 added

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdrep <- sdreport(obj)
sdrep

#saveRDS(obj, file = "data/fits/stock-prop-grouped-rw-RTMB.rds")
saveRDS(obj, file = "data/fits/stock-prop-grouped-rw-mapped-RTMB.rds")

m_svc <- update(m_svc, do_fit = TRUE, 
                control = sdmTMBcontrol(map = list(ln_tau_Z = factor(c(rep(1, 9), rep(2, 9))),
                                                   b_smooth = factor(b_smooth_map),
                                                   ln_smooth_sigma = factor(c(1,2,3,4,5,6,7,8,NA)))
                                        )
                )
                
m_svc$sd_report # compare with sdmTMB object
#saveRDS(m_svc, file = "data/fits/stock-prop-m-svc-sdmTMB.rds")

sanity(m_svc)

m_svc$tmb_random

cbind(summary(sdrep, "fixed"), summary(m_svc$sd_report, "fixed"))

# attempting to improve convergence with optimHess 
# to more closely match sdmTMB estimates and gradient
g <- as.numeric(obj$gr(opt$par))
h <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr)
opt$par <- opt$par - solve(h, g)
opt$objective <- obj$fn(opt$par)


# Predictive grid -----------------------------------------------------------

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
  

scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

unscale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x * sc_scale) + sc_center
}

scale_est(1:12, dat$month_adj)

pred_data_all <- index_grid %>%
  filter()
  replicate_df(index_grid, "month_adj", 1:12) %>%
  mutate(scale_month_adj = scale_est(month_adj, dat$month_adj)) %>%
  replicate_df("region", unique(dat$region))
glimpse(pred_data_all)


pred_data <- pred_data_all[1:10000,]

sm_pred <- sdmTMB:::parse_smoothers(stock_prop ~ 0 + region + (1 | year_adj_f) + 
                             s(scale_month_adj, by = region, bs = "tp", k = 6), 
                           data = dat,
                           newdata = pred_data, 
                           basis_prev = sm$basis_out)
# sm <- sdmTMB:::parse_smoothers(stock_prop ~  0 + region + (1 | year_f)
#                                + s(month_adj, by = region, bs = "tp", k = 6), dat)
tmb_newdata <- list(
  observed = numeric(nrow(pred_data)),
  covariate = pred_data$region,
  mesh = sdmTMB_mesh$mesh,
  offset = rep(median(dat$effort), nrow(pred_data)),
  spde_m =  sdmTMB:::get_spde_matrices(sdmTMB_mesh), # change names for Q_spde functions # sdmTMB_mesh$spde,
  spde_aniso = sdmTMB:::make_anisotropy_spde(sdmTMB_mesh),
  spde =  sdmTMB_mesh$spde, # change names for Q_spde functions # sdmTMB_mesh$spde,
  year_f = as.factor(pred_data$year_adj_f),
  interpolator_data = interpolator_data, # or A_st
  timesteps = sort(unique(pred_data$month_adj)),
  timesteps_et = 1:12,
  n_t = length(unique(pred_data$month_adj)), 
  n_tet = length(unique(pred_data$month_adj)), # +1 for extra_time
  t_i = as.numeric(pred_data$month_adj),
  nobs_RE = length(unique(pred_data$year_adj_f)),
  RE_indexes = as.numeric(as.factor(pred_data$year_adj_f)),
  mm_fixed = model.matrix(~ 0 + region, pred_data), # model matrix of fixed effects
  Zs         = sm_pred$Zs, # optional smoother basis function matrices
  Xs         = sm_pred$Xs, # optional s.moother linear effect matrix
  b_smooth_start = sm_pred$b_smooth_start,
  z_i = m_svc$tmb_data$z_i, # model.matrix(~ 0 + region, dat)
  n_z = ncol(m_svc$tmb_data$z_i),
  A_spatial_index = sdmTMB_mesh$sdm_spatial_id,
  n_c = length(unique(pred_data$region)), # number of categories
  c_i = as.integer(pred_data$region)
)

obj <- readRDS(file = "data/fits/stock-prop-grouped-rw-mapped-RTMB.rds")
obj$fn()


tmb_data <- tmb_newdata
lp <- obj$env$last.par.best
obj$report(par = lp)
summary(sdreport(obj))