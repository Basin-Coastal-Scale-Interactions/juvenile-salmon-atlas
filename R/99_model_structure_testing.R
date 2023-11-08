

### Model structure testing (spline knots, spatiotemporal RF type)

mout3_noknots <- sdmTMB(n_juv ~ 0 + s(month, bs = "cc", k = 12) +
                          day_night + survey_f + scale_dist + scale_depth + (1 | year_f),
                        offset = chinook_dat$effort,
                        #knots = list(month = c(0,12)),
                        data = chinook_dat,
                        mesh = spde,
                        family = sdmTMB::nbinom2(),
                        spatial = "on",
                        time = "month",
                        extra_time = c(1,4), # no need to add 1 as RW isn't cyclical
                        spatiotemporal = "AR1",
                        anisotropy = TRUE,
                        silent = FALSE
)

mout3_rw <- sdmTMB(n_juv ~ 0 + s(month, bs = "cc") +
                  day_night + survey_f + scale_dist + scale_depth + (1 | year_f),
                offset = chinook_dat$effort,
                data = chinook_dat,
                mesh = spde,
                family = sdmTMB::nbinom2(),
                spatial = "on",
                time = "month",
                extra_time = c(4), # no need to add 1 as RW isn't cyclical
                spatiotemporal = "RW",
                anisotropy = TRUE,
                silent = FALSE
)

mout3_iid <- sdmTMB(n_juv ~ 0 + s(month, bs = "cc") +
                   day_night + survey_f + scale_dist + scale_depth + (1 | year_f),
                 offset = chinook_dat$effort,
                 data = chinook_dat,
                 mesh = spde,
                 family = sdmTMB::nbinom2(),
                 spatial = "on",
                 time = "month",
                 extra_time = c(1,4), # no need to add 1 as RW isn't cyclical
                 spatiotemporal = "iid",
                 anisotropy = TRUE,
                 silent = FALSE
)
sanity(mout3)
sanity(mout3_rw)
sanity(mout3_iid)



pred_months_05 <- predict(mout3_05, newdata = grid_months, re_form_iid = NA) 
pred_months_nk <- predict(mout3_noknots, newdata = grid_months, re_form_iid = NA) 

pred_months_iid <- predict(mout3_iid, newdata = grid_months, re_form_iid = NA) 
pred_months_rw <- predict(mout3_rw, newdata = filter(grid_months, month != 1), re_form_iid = NA) 

pred_months_years <- predict(mout3, newdata = grid_months_years, re_form_iid = NA)

pred_months_comb <- rbind(mutate(pred_months, spatiotemporal = "AR1"),
                          mutate(pred_months_iid, spatiotemporal = "IID"),
                          mutate(pred_months_rw, spatiotemporal = "RW"))

pred_months_knots <- rbind(mutate(pred_months_nk, knots = "NULL"),
                           mutate(pred_months, knots = "c(0,12)"),
                           mutate(pred_months_05, knots = "c(0.5, 12.5)"))



plot_map(pred_months_comb, exp(est)) +
  scale_color_viridis_c(trans = fourth_root_power_trans()) +
  facet_grid(spatiotemporal~month_f) +
  ggtitle("Predicted distribution of juvenile chinook by month with different spatiotemporal random fields (no data for Jan and Apr)")
ggsave(here::here("figs", "seb", "chinook_juv_months_strf.png"),
       height = 8, width = 20, units = "in")


plot_map(pred_months_knots, exp(est)) +
  scale_color_viridis_c(trans = fourth_root_power_trans()) +
  facet_grid(knots~month_f) +
  ggtitle("Predicted distribution of juvenile chinook by month with different spatiotemporal random fields (no data for Jan and Apr)")
ggsave(here::here("figs", "seb", "chinook_juv_months_knots.png"),
       height = 8, width = 20, units = "in")
