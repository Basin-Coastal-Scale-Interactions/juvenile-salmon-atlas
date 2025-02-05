### Looking at Genetic Stock Identification (GSI) data for juvenile chinook
### and joining it with spatial data

library(readxl)
library(tidyverse)
library(ggsidekick) # for theme_sleek()
library(here)
library(janitor)

# reading unique event data (date, location, etc.) in catch data set
dat <- readRDS(here::here("data", "cleaned_all_bridge_catch_dat.rds"))

range(dat$date)

# extracting bridge data from juvenile salmon catch data set
bridge <- dat %>%
  filter(species == "chinook") %>%
  # select(unique_event:utm_y) %>%
  # unique() %>%
  as_tibble() #%>%
# # select(-ends_with("cycle")) %>%
# mutate(utm_x_1000 = utm_x / 1000,
#        utm_y_1000 = utm_y / 1000)

glimpse(bridge)

# reading GSI data (import has data from 1998 to 2021)
# gsi <- read_xlsx("data-raw/CHINOOK_event_catch_top_gsi_20240212.xlsx") %>%
#   janitor::clean_names()
# 
# gsi_mu <- read_xlsx("data-raw/CHINOOK_event_catch_top_gsi_20240212.xlsx",
#                             sheet = "catch_top_mu") %>%
#   janitor::clean_names() %>%
#   mutate(unique_event = gsub("-124J", "", unique_catch))
#          #year = as.numeric(sub("\\D*(\\d{4}).*", "\\1", unique_event)))
# 
# # same data set but arranged by fish rather than by tow
# gsi_by_fish <- read_xlsx("data-raw/CHINOOK_fish_top_gsi_20240212.xlsx",
#                          sheet = "top_mu") %>%
#   janitor::clean_names() %>%
#   mutate(unique_event = gsub("-124[JA]?-.+$", "", unique_fish)) 
# 
# gsi_by_fish[grep("124A", gsi_by_fish$unique_fish),] 
# checked this fish and this is a juvenile based on size that was originally
# labeled as an adult
#
# All percentages for all fish at MU level
# gsi_all_mu_in <- read_xlsx(here("data-raw", "CHINOOK_fish_all_gsi_20250122.xlsx"), sheet = "all_mu") %>%
#   janitor::clean_names()
# glimpse(gsi_all_mu_in)
# 
# grep("2023", perl = TRUE, gsi_all_mu_in$unique_fish) %>%
#   gsi_all_mu_in[., ]
# 
# All percentages for all fish at stock level
gsi_all_stock_in <- read_xlsx(here("data-raw", "CHINOOK_fish_all_gsi_20250122.xlsx"), sheet = "all_stock") %>%
  janitor::clean_names()
glimpse(gsi_all_stock_in)

# Key for unique_event
event_key <- read_xlsx(here("data-raw", "CHINOOK_fish_all_gsi_20250122.xlsx"), sheet = "fish_catch_event") %>%
  janitor::clean_names()

# Key for stock to Cam's regions

stock_key <- read_csv("data-raw/seb_key_2024-12-13.csv") 

filter(stock_key, is.na(juv_marine))

stock_key %>% filter(Region3Name == "North/Central BC") %>%
  pull(juv_marine) %>% unique()

filter(stock_key, grepl("^LEWIS", stock))

#gsi_all <- left_join(gsi_all_mu_in, event_key, by = "unique_fish") 

stock_key %>% filter(stock == "TETE_JAUNE") %>%
  pull(juv_marine) %>% unique()


# Joining GSI data with event key and then joining with stock key
gsi_all_stock <- left_join(gsi_all_stock_in, event_key, by = "unique_fish") %>%
  select(-stock) %>%
  rename(stock = stock_final) %>%
  mutate(stock = if_else(stock == "LEWIS_ LAKE", "LEWIS_LAKE", stock)) %>% # Renaming misspelled stock
  left_join(select(stock_key, stock, juv_marine), by = "stock") %>% # SELECTING REGION3 HERE
  select(unique_event, unique_fish, stock, region = juv_marine, stock_prob)
  # Need to define "region" to either Region2Name or Region3Name if using other stock key
  # left_join(select(stock_key, stock, Region1Name, Region2Name, Region3Name), by = "stock") %>% # SELECTING REGION3 HERE
  # select(unique_event, unique_fish, region = Region3Name, stock_prob)

# checking for NA regions that aren't unknowns or NAs
filter(gsi_all_stock, is.na(region)) %>% pull(stock) %>% unique()

filter(gsi_all_stock, unique_fish == "HS200215-VI01-124-025")


barkley_ue <- readRDS(here("data", "barkley_unique_events.rds"))
barkley_stock <- left_join(barkley_ue, gsi_all_stock, by = "unique_event") %>%
  mutate(month = month(date)) %>%
  filter(month %in% 6:7)

barkley_fish <- filter(barkley_stock, unique_fish %in% filter(barkley_stock, region == "Salish Coastal")$unique_fish)

barkley_fish %>% filter(region == "Salish Coastal")

filter(barkley_fish, unique_fish %in% 
         filter(barkley_fish, region == "Salish Coastal" & 
                  stock_prob < 30 & 
                  stock_prob > 10)$unique_fish) %>%
  print(n = 25)

# table(gsi_all_stock$Region1Name)
# table(gsi_all_stock$Region2Name)
# table(gsi_all_stock$Region3Name)

table(gsi_all_stock$region)


baselines_tb <- gsi_all_stock %>% 
  drop_na() %>%
  filter(region != "Cali", stock_prob >= 20) %>%
  # mutate(Region = fct_recode(region,
  #                     "Fraser Summer 4.1" = "Fraser_Summer_4.1",
  #                     "Northern BC/AK" = "NBC_SEAK"           
  # )) %>%
  mutate(Region = fct_relevel(region,
                              c("NBC/SEAK", "Central BC", "Salish Coastal", 
                                "WCVI", "Fraser Yearling",   "Fraser Summer 4.1", 
                                "WA/OR Coastal", "Col. Coastal", "Col. North", 
                              "Upper Col. Yearling"))) %>%
  group_by(Region) %>%
  reframe(Baselines = str_to_title(str_replace_all(unique(stock), "_", " "))) %>% 
  group_split(group_id = c(rep(1:3, times = 70),1)) %>%  # nrow 71, 70, 70
  map(function (x) select(x, -group_id)) %>% 
  map(function (x) if (nrow(x) == 70) {
    bind_rows(x, tibble(Region = " ", Baselines = " "))
    } else {
      return(x) 
  })  %>%
   bind_cols(.name_repair = "minimal")



library(kableExtra)
baselines_kb <- knitr::kable(list(baselines_tb), format = "latex", align = "|ll|ll|ll|",
             caption = "Genetic stock identification (GSI) assignment baselines detected\
             in all sampled juvenile chinook salmon with values higher than 0.2.") %>%
  kable_styling(font_size = 6)


writeLines(baselines_kb, here("tables","baselines.tex"))


library(xtable)
xtable(baselines_tb)

# gsi_all_stock %>%
#   select(Region1Name, Region2Name, Region3Name) %>%
#   distinct() %>% print(n = 46)
# 
# gsi_all_stock %>%
#   select(Region1Name, Region2Name, Region3Name) %>%
#   group_by(Region1Name, Region2Name) %>%
#   reframe(count = n(), 
#             Region3Name = unique(Region3Name) )%>% 
#   print(n = 46)

#nas_mu <- gsi_all[is.na(gsi_all$mu_name),]
nas_stock <- gsi_all_stock[is.na(gsi_all_stock$region), ]

unique(nas_stock$stock)
table(nas_stock$stock)
filter(nas_stock, is.na(stock))

# table(gsi_all$mu_name)
# levels(factor(gsi_all$mu_name))

table(gsi_all_stock$region)
levels(factor(gsi_all_stock$region))

# Only 9 fish are from Cali region and have a GSI probability greater than 40%
# of which half are from hatchery origins
gsi_all_stock %>%
  filter(region == "Cali" & stock_prob > 40) %>% print(n = 68)

##########

# get unique_fish values for GSIed fish that have a max stock_prop less than 0.5
low_gsi <- gsi_all_stock %>% filter(!is.na(region)) %>% 
  select(-stock) %>%
  #mutate(region = fct_other(region, drop = "Fraser River", other_level = "SOG")) %>% # Combining Fraser into SOG
  complete(region, nesting(unique_event, unique_fish)) %>%
  mutate(stock_prop = if_else(is.na(stock_prob), 0, stock_prob/100)) %>%
  group_by(unique_fish) %>%
  summarize(max_prop = max(stock_prop)) %>%
  filter(max_prop < 0.5) %>%
  pull(unique_fish)

# Wrangling data so that there are columns for all regions  with zeroes for 
# each unique event where GSI regions are absent
gsi_all_stock_grouped <- gsi_all_stock %>% filter(!is.na(region)) %>% 
  select(-stock) %>%
  #mutate(region = fct_other(region, drop = "Fraser River", other_level = "SOG")) %>% # Combining Fraser into SOG
  complete(region, nesting(unique_event, unique_fish)) %>%
  mutate(stock_prop = if_else(is.na(stock_prob), 0, stock_prob/100)) %>% # converting NAs to zeros
  filter(!unique_fish %in% low_gsi) %>% # removing fish with max stock prop lower than 0.5
  group_by(unique_event, region) %>% # adding probabilities for across all fish in each unique_event
  summarize(stock_prop = sum(stock_prop), n_fish = length(unique(unique_fish))) %>%
  ungroup() %>%
  arrange(unique_event, region) 

gsi_all_stock %>% filter(!is.na(region)) %>% 
  select(-stock) %>%
  #mutate(region = fct_other(region, drop = "Fraser River", other_level = "SOG")) %>% # Combining Fraser into SOG
  complete(region, nesting(unique_event, unique_fish)) %>%
  mutate(stock_prop = if_else(is.na(stock_prob), 0, stock_prob/100)) %>%
  group_by(unique_fish) %>%
  slice(which.max(stock_prob)) %>%
  #filter(stock_prop >= .5) %>%
  ungroup() %>%
  group_by(unique_event) %>%
  tally()

gsi_all_stock %>% filter(!is.na(region)) %>% 
  select(-stock) %>%
  #mutate(region = fct_other(region, drop = "Fraser River", other_level = "SOG")) %>% # Combining Fraser into SOG
  complete(region, nesting(unique_event, unique_fish)) %>%
  mutate(stock_prop = if_else(is.na(stock_prob), 0, stock_prob/100)) %>%
  group_by(unique_fish) %>%
  slice(which.max(stock_prop)) %>%
  ggplot() +
  geom_histogram(aes(stock_prop), binwidth = 0.05) +
  theme_sleek() +
  geom_vline(aes(xintercept = .50))

gsi_all_stock_grouped %>%
  group_by(unique_event) %>%
  summarise(ratio = sum(stock_prop)/unique(n_fish)) %>%
  pull(ratio) %>% hist()

lowratios <- gsi_all_stock_grouped %>%
  group_by(unique_event) %>%
  summarise(ratio = sum(stock_prop)/unique(n_fish)) %>%
  filter(ratio < 0.4) %>% pull(unique_event)

highratios <- gsi_all_stock_grouped %>%
  group_by(unique_event) %>%
  summarise(ratio = sum(stock_prop)/unique(n_fish)) %>%
  filter(ratio > 1) %>% pull(unique_event)

# Looking at unique fish with total GSI probabilities higher than 1 and
# lower than 0.4
filter(gsi_all_stock, unique_event %in% lowratios) %>% print(n = 90)
filter(gsi_all_stock, unique_event %in% highratios) %>% print(n = 90)


filter(gsi_all_stock, unique_event == "BCSI-201731-HW01")
#filter(gsi_all, unique_event == "BCSI-201731-HW01")

# Removing "Cali" region here
chinook_gsi_counts <- inner_join(gsi_all_stock_grouped, bridge, by = "unique_event") %>%
  filter(region != "Cali") %>%
  mutate(
    region = as.factor(region),
    year_f = as.factor(year),
    month_f = month(date, label = TRUE),
    month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3),
   # scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night))

glimpse(chinook_gsi_counts)

table(chinook_gsi_counts$region)

saveRDS(chinook_gsi_counts, 
        file = here("data", "chinook_gsi_counts.rds"))

filter(chinook_gsi_counts, year %in% 2022:2024)

##########################


gsi_all_stock_grouped %>%
  group_by(unique_event) %>%
  summarise(n_fish = unique(n_fish)) %>%
  ggplot(.) + geom_histogram(aes(n_fish), binwidth = 1) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(breaks = seq(0, 55, by = 5)) +
  ggtitle("Distribution of number of fish in GSI samples by tow in catch dataset")
ggsave(here("figs", "gsi", "hist-gsi-samples-by-tow.png"), width = 9, height = 5)

chinook_gsi_counts %>%
  ggplot(data = .) +
  geom_col(aes(x = week,  y = stock_prop, fill = region), position = "stack") +
  #coord_flip() +
  #facet_wrap(~month, scales = "free_y") +
  theme_sleek() +
  NULL
ggsave(here("figs", "gsi", "gsi-chinook-by-week.png"), width = 9, height = 5)

month_labels <- c('J','F','M','A','M','J','J','A','S','O','N','D')

# faceting by year
chinook_gsi_counts %>%
  mutate(month_f =  month(month, label = TRUE)) %>%
  ggplot(data = .) +
  geom_col(aes(x = month_f, y = stock_prop, fill =region), position = "stack") +
  #coord_flip() +
  facet_wrap(~year, scales = "free_y") +
  scale_x_discrete(labels = month_labels, drop = FALSE) +
  theme_sleek() +
  NULL
ggsave(here("figs", "gsi", "gsi-chinook-by-month-year.png"), width = 9, height = 8)



# faceting by year but with fixed y scale
chinook_gsi_counts %>%
  mutate(month_f =  month(month, label = TRUE)) %>%
  ggplot(data = .) +
  geom_col(aes(x = month_f, y = stock_prop, fill =region), position = "stack") +
  #coord_flip() +
  facet_wrap(~year, scales = "fixed") +
  scale_x_discrete(labels = month_labels, drop = FALSE) +
  theme_sleek() +
  NULL
ggsave(here("figs", "gsi", "gsi-chinook-by-month-year-fixed-y.png"), width = 9, height = 8)


gsi_all_stock_grouped %>%
  group_by(unique_event) %>%
  summarise(n_fish = unique(n_fish)) %>%
  group_by(n_fish) %>% count(name = "ncount") %>%
  ungroup() %>%
  mutate(ntotal = n_fish * ncount,
         cumcount = cumsum(ntotal),
         cumprop = cumcount/sum(ntotal)) %>%
  ggplot(.) + geom_col(aes(x = n_fish, y = ntotal)) +
  geom_line(aes(x = n_fish, y = cumprop*1300)) + geom_point(aes(x = n_fish, y = cumprop*1300)) +
  scale_x_continuous(breaks = seq(0, 55, by = 5)) +
  scale_y_continuous("Fish count per number of samples", breaks = seq(0, 1400, by = 200),
                     sec.axis = sec_axis(~./1300, name = "Cumulative proportion")) +
  ggtitle("Distribution of number of fish in GSI samples by tow in catch dataset") +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "hist-gsi-total-fish-by-tow.png"), width = 9, height = 5)


# Plotting locations of individual data points by month color coded by MU
source(here("R", "plot_map.R"))

chinook_gsi_counts %>%
    ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = stock_prop)) +
  xlim(range(chinook_gsi_counts$utm_x)) +
  ylim(range(chinook_gsi_counts$utm_y)) +
  facet_wrap(~season_f, nrow = 2) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek()

ggsave(here("figs", "gsi", "gsi-locations-by-season.png"), width = 9, height = 6)

chinook_gsi_counts %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = stock_prop)) +
  xlim(range(chinook_gsi_counts$utm_x)) +
  ylim(range(chinook_gsi_counts$utm_y)) +
  xlab("")  + ylab("") +
  facet_wrap(~month_f, drop = FALSE) +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "gsi-locations-by-month.png"), width = 9, height = 6)


chinook_gsi_counts %>%
  full_join(crossing(season_f = unique(chinook_gsi_counts$season_f), year=unique(chinook_gsi_counts$year))) %>% 
  mutate(empty=ifelse(is.na(stock_prop), "No samples", NA_character_),
         x=mean(range(utm_x, na.rm=TRUE)), 
         y=mean(range(utm_y, na.rm=TRUE))) %>%
  ungroup() %>%
  select(-month_f) %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = stock_prop)) +
  facet_grid(season_f~year) + # , margins = "year") + # margins doesn't work with geom_sf()
  geom_text(aes(x, y, label=empty), colour="red", size=3) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  xlim(range(chinook_gsi_counts$utm_x)) +
  ylim(range(chinook_gsi_counts$utm_y)) +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "gsi-locations-season-by-year.png"), width = 40, height = 5)

library(gridExtra)

png(here("figs", "gsi", "gsi-locations-season-by-year-2.png"), width = 20, height = 9, units = "in", res = 200)

grid.arrange(

chinook_gsi_counts %>%
  full_join(crossing(season_f = unique(chinook_gsi_counts$season_f), year=unique(chinook_gsi_counts$year))) %>% 
  mutate(empty=ifelse(is.na(stock_prop), "No samples", NA_character_),
         x=mean(range(utm_x, na.rm=TRUE)), 
         y=mean(range(utm_y, na.rm=TRUE))) %>%
  ungroup() %>%
  filter(year %in% 1998:2009) %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = stock_prop)) +
  xlim(range(chinook_gsi_counts$utm_x)) +
  ylim(range(chinook_gsi_counts$utm_y)) +
  facet_grid(season_f~year) +
  geom_text(aes(x, y, label=empty), colour="red", size=3) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek(),

chinook_gsi_counts %>%
  full_join(crossing(season_f = unique(chinook_gsi_counts$season_f), year=unique(chinook_gsi_counts$year))) %>% 
  mutate(empty=ifelse(is.na(stock_prop), "No samples", NA_character_),
         x=mean(range(utm_x, na.rm=TRUE)), 
         y=mean(range(utm_y, na.rm=TRUE))) %>%
  ungroup() %>%
  filter(year %in% 2010:2021) %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = stock_prop)) +
  xlim(range(chinook_gsi_counts$utm_x)) +
  ylim(range(chinook_gsi_counts$utm_y)) +
  facet_grid(season_f~year) +
  geom_text(aes(x, y, label=empty), colour="red", size=3) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek(),
nrow=2)

dev.off()

