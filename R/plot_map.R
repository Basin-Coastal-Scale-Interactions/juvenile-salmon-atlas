# shape file for coastline
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -129, ymin = 48.25, xmax = -123.25, ymax = 51.3) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))


plot_map <- function(dat, column, show_raw_data = FALSE, raw_dat) {
  pos_pt_fill = "#FFFFFF05"
  g <- ggplot() + # (dat, aes(X, Y, color = {{ column }})) +
    geom_raster(
      data = dat,
      aes(X, Y, fill = {{ column }})) +
    scale_fill_viridis_c(name = "Predicted\nnumber of\njuveniles\nper tow",
                         labels = comma,
                         trans = fourth_root_power_trans()) +
    ggsidekick::theme_sleek() +
    coord_fixed() +
    geom_sf(data = coast, color = "gray80", fill = "gray90")
  # geom_point(data = mutate(chinook_dat, month_f = month(month, label = TRUE)),
  #            aes(utm_x, utm_y, colour = if_else(n_juv == 0, "black", "red")), size =0.1) +
  
  if (show_raw_data) {
    ex <- 4
    g <- g +
      geom_point(
        data = filter(raw_dat, n_juv == 0),
        aes(utm_x, utm_y, pch = as.factor(ex)),
        pch = 4,
        col = "grey50",
        size = 0.5, alpha = 0.6
      ) + 
      scale_shape_manual(values = 0, name = "Observed\nnumber of\njuveniles\nper tow") +
      geom_point(
        data = filter(raw_dat, n_juv > 0),
        aes(utm_x, utm_y,
            size = n_juv,
        ), # fill = pos_pt_fill,
        col = "red", pch = 21, alpha = 0.6
      ) + 
      scale_size_binned_area(name = "Observed\nnumber of\njuveniles\nper tow",
                             labels = comma,
                             limits = c(1,2500),
                             breaks = c(1,10, 100, 500, 1500)) #,
    # trans = "log10")
    # scale_size_area(name = "Observed\nnumber of\njuveniles\nper tow",
    #                        labels = comma,
    #                 breaks = c(10, 100, 500, 1500)) #,
    #                 # trans = "log10")
  }
  
  g + 
    facet_wrap(~month_f) +
    scale_x_continuous(name = NULL, limits = c(462700, 814800), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
  
}
