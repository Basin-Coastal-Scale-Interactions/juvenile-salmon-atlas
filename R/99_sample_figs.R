### Plots of data summary

library(tidyverse)
library(cowplot)
library(sf)
library(knitr)
library(kableExtra)
options(knitr.kable.NA = '')


chinook_dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))        
gsi_dat <- readRDS(here::here("data", "chinook_gsi_counts_fitted.rds")) %>%
  mutate(scale_month_adj = scale(month_adj)[, 1],
         year_adj = ifelse(month > 2, year, year - 1), # adjusting year to represent fish cohorts
         # year_f = as.factor(year),
         year_adj_f = as.factor(year_adj))
m_svc <- readRDS(here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))
m_svc$version
gsidat <- m_svc$data
table(gsidat$region)


min_lat <- min(floor(c(chinook_dat$lat, gsi_dat$lat)) - 0.25)
max_lat <- max(c(chinook_dat$lat, gsi_dat$lat)) + 0.25
min_lon <- min(floor(c(chinook_dat$lon, gsi_dat$lon)) - 0.25)
max_lon <- max(c(chinook_dat$lon, gsi_dat$lon)) + 0.25

coast <- rbind(rnaturalearth::ne_states( "United States of America",
                                         returnclass = "sf"),
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>%
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

#scales::show_col(scales::hue_pal()(3))

# Species wide tows
set_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "gray90", size = 1.25) +
  geom_point(data = chinook_dat,
             aes(x = lon, y = lat), fill = "#F8766D", # fill = survey_f (in aes())
             shape = 21, alpha = 0.4) +
  #scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
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

# Stock-specific tows
set_map2 <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "gray90", size = 1.25) +
  geom_point(data = gsi_dat,
             aes(x = lon, y = lat), fill = "#F8766D", # fill = survey_f (in aes())
             shape = 21, alpha = 0.4) +
  #scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
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
  cowplot::draw_plot(inset_map, x = 0.097, y = 0.095, #vjust = -0.2,
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
  geom_jitter(aes(y = week, x = year, size = n_tows), 
                  #shape = season_f2,
                  color = "#F8766D",
              alpha = 0.3, width = 0.25) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number\nof tows") +
  #scale_color_discrete(name = "Survey", labels = c("HSS", "IPES")) +
  labs(x = "Year", y = "Week") +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  )

bubble_combined <- bind_rows(chinook_dat %>% 
            group_by(year, week) %>% 
            summarize(n_tows = length(unique(unique_event)), .groups = "drop") %>%
            mutate(model = "Species-wide model"),
            
          gsi_dat %>% 
            group_by(year, week) %>% 
            summarize(n_tows = length(unique(unique_event)), .groups = "drop") %>%
            mutate(model = "Stock-specifc model")) %>%
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(y = week, x = year, size = n_tows), 
                  color = "#F8766D",
              alpha = 0.3, width = 0.25) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number\nof tows") +
 # scale_color_discrete(name = "Survey", labels = c("HSS", "IPES")) +
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

### Figures of spatial tow coverage by month

#### Chinook abundance data

# adding missing month factor levels
cd2 <- chinook_dat %>%
  mutate(month_f = ordered(month_f, levels = month(1:12, abbr = TRUE, label = TRUE)))

nodata_text <- data.frame(
  label = c("No data", rep("", 11)),
  month_f   = ordered(levels(cd2$month_f))
)

ch_monthly <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "gray90", size = 1.25) +
  geom_point(data = cd2,
             aes(x = lon, y = lat), fill = "#F8766D", 
             shape = 21, alpha = 0.4) +
  #scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0), labels = c("136°W", "","132°W", "","128°W", "","124°W") ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(#legend.position = "outside", #legend.position.inside = c(0.9, 0.8),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13)) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15), 
           xlim = c(min_lon + 0.15, max_lon - 0.15)) +
  facet_wrap(~month_f, nrow = 3, drop = FALSE)   +
  geom_text(data = nodata_text, mapping = aes(x = -130, y = 48, label = label), size = 6)

ggsave(here::here("figs", "chinook_data_by_month.png"), 
       ch_monthly, width = 10, height = 8, units = "in")

#### GSI sampling monthly coverage

# adding missing month factor levels
gd2 <- gsi_dat %>%
  mutate(month_f = ordered(month_f, levels = month(1:12, abbr = TRUE, label = TRUE)))

nodatagsi_text <- data.frame(
  label = c("No data", "", "", "No data", rep("", 8)),
  month_f   = ordered(levels(cd2$month_f))
)

gsi_monthly <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "gray90", size = 1.25) +
  geom_point(data = gd2,
             aes(x = lon, y = lat), fill = "#F8766D", 
             shape = 21, alpha = 0.4) +
 # scale_fill_discrete(name = "Survey", labels=c("HSS", "IPES")) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0), labels = c("136°W", "","132°W", "","128°W", "","124°W") ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(#legend.position = "outside", #legend.position.inside = c(0.9, 0.8),
    legend.title = element_text(size=15),
    legend.text = element_text(size=13)) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15), 
           xlim = c(min_lon + 0.15, max_lon - 0.15)) +
  facet_wrap(~month_f, nrow = 3, drop = FALSE)   +
  geom_text(data = nodatagsi_text, mapping = aes(x = -130, y = 48, label = label), size = 6)

ggsave(here::here("figs", "gsi_data_by_month.png"), 
       gsi_monthly, width = 10, height = 8, units = "in")


### Summary table of summed GSI proportions by stock grouping and months

gsi_tb <- gsidat %>%
  group_by(month_f, region) %>%
  summarise(sum = sum(stock_prop)) %>%
  pivot_wider(names_from = month_f, values_from = sum) %>%
  ungroup() %>%
  group_by(region) %>%
  mutate(Total = sum(Feb+Mar+Jun+Jul+Aug+Sep+Oct+Nov+Dec)) %>% 
  mutate(Jan = 0, Apr = 0, May = 0) %>%
  select(Region = region, month(1:12, abbr = TRUE, label = TRUE), Total) %>%
  janitor::adorn_totals("row") %>% 
  mutate_if(is.numeric, ~round(., 1))

gsi_kb <- kable(gsi_tb, format = "latex", align = "lrrrrrrrrrrrrr",
                label = "gsi-summary", , booktabs = TRUE,
                linesep = c(rep("", 9), "\\addlinespace"),
                caption = "Summary table of summed proportion of GSI assigment
                by region and month.") %>%
  # kable_styling(font_size = 6) %>%
  kable_styling(latex_options = "hold_position")

writeLines(gsi_kb, here("tables","gsi_summary.tex"))
