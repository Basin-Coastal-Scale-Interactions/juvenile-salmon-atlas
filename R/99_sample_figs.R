### Plots of data summary

library(tidyverse)
library(cowplot)
library(sf)

chinook_dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))        
gsi_dat <- readRDS(here::here("data", "chinook_gsi_counts_20240722.rds")) %>%
  mutate(season_n = as.numeric(season_f),
         scale_month_adj = scale(month_adj)[, 1],
         year_adj = ifelse(month > 2, year, year - 1), # adjusting year to represent fish cohorts
         # year_f = as.factor(year),
         year_adj_f = as.factor(year_adj))


min_lat <- min(floor(c(chinook_dat$lat, gsi_dat$lat)) - 0.25)
max_lat <- max(c(chinook_dat$lat, gsi_dat$lat)) + 0.25
min_lon <- min(floor(c(chinook_dat$lon, gsi_dat$lon)) - 0.25)
max_lon <- max(c(chinook_dat$lon, gsi_dat$lon)) + 0.25

coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

set_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "gray90", size = 1.25) +
  geom_point(data = chinook_dat,
             aes(x = lon, y = lat, fill = survey_f), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13)) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15), 
           xlim = c(min_lon + 0.15, max_lon - 0.15))
set_map

set_map2 <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "gray90", size = 1.25) +
  geom_point(data = gsi_dat,
             aes(x = lon, y = lat, fill = survey_f), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13)) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15), 
           xlim = c(min_lon + 0.15, max_lon - 0.15))
set_map2

# inset
w_can <- map_data("world", region = c("usa", "canada")) %>%
  fortify(.)

state_prov <- rnaturalearth::ne_states(c("united states of america", "canada")) 

inset_map <- ggplot() +
  geom_sf(data = state_prov, color = "black", fill = "grey90") + 
  #geom_sf_text(data = state_prov, aes(label = postal)) +
  labs(x = "", y = "") +
  geom_rect(aes(xmin = min_lon, xmax = max_lon, ymin = min_lat, ymax = max_lat),
            fill = NA, lty = 2, col = "red") +
  coord_sf(ylim = c(40, 65), xlim = c(-150, -118), crs = 3005, # BC Albers
           default_crs = 4326) + # Defaults WGS 84 code
  ggsidekick::theme_sleek() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "top",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) 


# png(here::here("figs", "ms_figs_season_mvrw", "set_map.png"),
#     height = 4.5, width = 6, 
#     units = "in", res = 250)
p2 <- cowplot::ggdraw(set_map2) +
  cowplot::draw_plot(inset_map, x = 0.097, y = 0.065, #vjust = -0.2,
                     hjust = 0.1,
                     width = 0.4, height = 0.35)

ggsave(file = here::here("figs","set_maps.png"), 
       plot = plot_grid(set_map, p2,  vjust = 7, hjust = -1.6,
                        labels = c('a) Species-wide model', 'b) Stock-specific model')), 
       width = 12, height = 7.8, units = "in")

# dev.off()


### Bubble temporal coverage plots

# bubble plots of temporal coverage
bubble_temp_coverage <- chinook_dat %>% 
  group_by(year, week, survey_f) %>% 
  summarize(n_tows = length(unique_event), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(y = week, x = year, size = n_tows, 
                  #shape = season_f2,
                  color = survey_f),
              alpha = 0.3, width = 0.25) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number\nof tows") +
  scale_color_discrete(name = "Survey", labels = c("HSS", "IPES")) +
  labs(x = "Year", y = "Week") +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  )

bubble_combined <- bind_rows(chinook_dat %>% 
            group_by(year, week, survey_f) %>% 
            summarize(n_tows = length(unique(unique_event)), .groups = "drop") %>%
            mutate(model = "Species-wide model"),
            
          gsi_dat %>% 
            group_by(year, week, survey_f) %>% 
            summarize(n_tows = length(unique(unique_event)), .groups = "drop") %>%
            mutate(model = "Stock-specifc model")) %>%
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(y = week, x = year, size = n_tows, 
                  color = survey_f),
              alpha = 0.3, width = 0.25) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number\nof tows") +
  scale_color_discrete(name = "Survey", labels = c("HSS", "IPES")) +
  labs(x = "Year", y = "Week") +
  facet_wrap(~model, ncol = 1) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  )

ggsave(file = here::here("figs","temp_cov.png"), 
       plot = bubble_combined, 
       width = 7, height = 7, units = "in")
