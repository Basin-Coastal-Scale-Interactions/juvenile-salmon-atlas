# shape file for coastline
library(scales)
library(ggplot2)
library(ggspatial)
library(ggsidekick)
library(rnaturalearth)
library(rnaturalearthhires)

# coast outline crop
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                               returnclass = "sf"), 
                     rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -129, ymin = 48.25, xmax = -123.25, ymax = 51.3) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# as above but renamed for redundancy
south_coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -129, ymin = 48.25, xmax = -123.25, ymax = 51.3) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# new out
all_coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                               returnclass = "sf"), 
                     rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -138, ymin = 46.25, xmax = -119, ymax = 60) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# Need coast polygon in km for some graphs
all_coast_km <- rbind(rnaturalearth::ne_states( "United States of America", 
                                                returnclass = "sf"), 
                      rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -138, ymin = 46.25, xmax = -119, ymax = 60) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=km"))

# inlet specific
barkley <- rnaturalearth::ne_states( "Canada", returnclass = "sf") %>% 
  sf::st_crop(., xmin = -126, ymin = 48.6, xmax = -124.5, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# function to quickly provide plus and minus values around given value
# useful for plotting around specific UTM coordinates
pm <- function(x, pm) {
  return(c(x - pm, x + pm))
}

# Main map plotting function
plot_map <- function(dat, 
                     column, 
                     show_raw_data = FALSE,
                     raw_dat, 
                     area = "WCVI",
                     show_coast = TRUE) {
  pos_pt_fill = "#FFFFFF05"
  
  if (area == "Barkley") {
    maplayer <- barkley 
    } else {
      maplayer <- all_coast 
    }
  
  g <- ggplot() + # (dat, aes(X, Y, color = {{ column }})) +
    geom_raster(
      data = dat,
      aes(X, Y, fill = {{ column }})) +
    scale_fill_viridis_c(name = "Predicted\nnumber of\njuveniles\nper tow",
                         labels = comma,
                         trans = fourth_root_power_trans()) +
    ggsidekick::theme_sleek() +
    coord_fixed()
  
  if (show_coast == TRUE) {
    g <- g + geom_sf(data = maplayer, color = "gray80", fill = "gray90")
  }
  if (show_coast == "coastline") { 
    g <- g + geom_sf(data = maplayer, color = "gray80", fill = NA)
  }
  
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
  
  if (area == "WCVI") {
    g <- g + scale_x_continuous(name = NULL, 
                            limits = c(462700, 814800), expand = c(0, 0)) +
      scale_y_continuous(name = NULL, limits = c(5350050, 5681850), expand = c(0, 0))
  }
  
  if (area == "BC") {
    g <- g + scale_x_continuous(name = NULL, 
                                limits = c(100000, 875000), expand = c(0, 0)) +
      scale_y_continuous(name = NULL, limits = c(5350000, 6100000), expand = c(0, 0))
  }  
  
  if (area == "Barkley") {
    g <- g + scale_x_continuous(name = NULL, 
                                limits = c(718406, 812406), expand = c(0, 0)) +
      scale_y_continuous(name = NULL, limits = c(5388664, 5477664), expand = c(0, 0))
  }  
  
  if (area == "QCS") {
    g <- g + scale_x_continuous(name = NULL, 
                                limits = c(292406, 5824664), expand = c(0, 0)) +
      scale_y_continuous(name = NULL, limits = c(5538664, 5872664), expand = c(0, 0))
  }  
  
  
  g + 
    facet_wrap(~fct_reorder(month_f, month_adj)) +
    theme(legend.key.height = unit(0.06, 'npc'))
}


# from Philina
plot_mesh <- function(mesh_obj = spde2b,
                      data_obj = dat,
                      catch_var = "n_juv",
                      group = "number of juveniles") {
  ggplot() +
    inlabru::gg(mesh_obj$mesh) +
    coord_fixed() +
    geom_point(aes(utm_x_1000, utm_y_1000),
               shape = "x",
               size = 0.75,
               data = data_obj,
               inherit.aes = FALSE) +
    NULL +
    geom_point(
      aes(
        utm_x_1000, utm_y_1000,
        #size = .data[[catch_var]],
        #shape = survey_type,
        fill = .data[[catch_var]],
        #colour = .data[[catch_var]],
        NULL
      ),
      data = filter(data_obj, .data[[catch_var]] > 0)
    ) +
    facet_wrap( ~ species) +
    # scale_shape_discrete() +
    scale_fill_viridis_c(trans = "fourth_root_power") +
    scale_color_viridis_c(trans = "fourth_root_power") +
    #ggtitle(paste0(species, " (", group, ")")) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}

# plot_mesh(mesh_obj = bmesh01,
#           data_obj = dat,
#           catch_var = "n_juv",
#           group = "number of juveniles")
