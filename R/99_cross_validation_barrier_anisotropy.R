### Cross-validation

library(sdmTMB)
library(sdmTMBextra)
library(tidyverse)


dat_tbl_mout3 <- readRDS(here("data", "fits", "all_mout3_allcoast_withpred.rds"))
dat_tbl_barrier <- readRDS(here::here("data", "fits", "dat_tbl_barrier.rds"))


dat_tbl_mout3$fit[[1]]$formula


fit <- sdmTMB_cv(
  dat_tbl_mout3$fit[[1]]$formula[[1]],
  data = dat_tbl_mout3$data[[1]],
  mesh = spde2,
  family = sdmTMB::nbinom2(),
  k_folds = 8,
  fold_ids = NULL,
  parallel = TRUE,
  use_initial_fit = FALSE
)

fitb1 <- sdmTMB_cv(
  dat_tbl_barrier$fit[[1]]$formula[[1]],
  data = dat_tbl_barrier$data[[1]],
  mesh = spde2b,
  family = sdmTMB::nbinom2(),
  k_folds = 8,
  fold_ids = NULL,
  parallel = TRUE,
  use_initial_fit = FALSE
)

fitb2 <- sdmTMB_cv(
  dat_tbl_barrier$fit2[[1]]$formula[[1]],
  data = dat_tbl_barrier$data[[1]],
  mesh = spde2b2,
  family = sdmTMB::nbinom2(),
  k_folds = 8,
  fold_ids = NULL,
  parallel = TRUE,
  use_initial_fit = FALSE
)

fit$fold_loglik
fit$sum_loglik
fit$converged

fitb1$fold_loglik
fitb1$sum_loglik
fitb1$converged

fitb2$fold_loglik
fitb2$sum_loglik
fitb2$converged
