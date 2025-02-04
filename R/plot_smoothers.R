### plotting smoothers

# functions for quickly applying and unapplying scale()
scale_est <- function (x, unscaled_vec) {
  sc <- scale(unscaled_vec)
  sc_center <- attr(sc, "scaled:center")
  sc_scale <- attr(sc, "scaled:scale")
  (x - sc_center) / sc_scale
}

plot_smoother <- function (fit, regions = FALSE, species) {
  # create data frame with time steps being smoothed to predict with 
  nd <- data.frame(
    month = 1:12L,
    survey_f = "hss",
    day_night = "DAY",
    scale_depth = 0,
    scale_dist = 0,
    year_adj_f = "2010"
  ) %>%
    mutate(month_adj = case_when(species == "chum" ~
                                   ifelse(month > 2, # adjusting to start in April
                                          month - 2, 
                                          month + 10),
                                 species %in% c("coho", "pink","sockeye","chinook") ~ 
                                   ifelse(month > 3, # adjusting to start in April
                                          month - 3, 
                                          month + 9)),
           month_f =  month(month, label = T, abbr = T),
           scale_month_adj = scale_est(month_adj, fit$data$month_adj))
  
  
  # region needs to be a factor otherwise 
  if (regions) { 
    nd <- replicate_df(nd, "region", as.factor(levels(fit$data$region)))
  }
  
  # predicting using NA for both re_form and re_form_iid
  p <- predict(fit, 
               newdata = nd, se_fit = TRUE, offset = rep(-6.287114, nrow(nd)),
               re_form = NA, re_form_iid = NA)
  
  # plotting
  gp <- ggplot(p, aes(month_adj, exp(est),
                      ymin = exp(est - 1.96 * est_se),
                      ymax = exp(est + 1.96 * est_se))) +
    geom_line() +
    geom_ribbon(alpha = 0.4) +
    scale_x_continuous(labels = month(c(4:12, 1:3), label = T, abbr = T), breaks = 1:12) +
    coord_cartesian(expand = F) +
    ggsidekick::theme_sleek() + 
    labs(x = "Month", y = "Count")
  
  if (regions) {
    gp + facet_wrap(~region)
  } else {
    gp
  }
}
