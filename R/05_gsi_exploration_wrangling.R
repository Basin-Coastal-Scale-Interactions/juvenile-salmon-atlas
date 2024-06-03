### Looking at Genetic Stock Identification (GSI) data for juvenile chinook
### and joining it with spatial data

library(readxl)
library(tidyverse)
library(ggsidekick) # for theme_sleek()
library(here)
library(janitor)

# reading unique event data (date, location, etc.) in catch data set
#dat_in <- readRDS(here::here("data", "catch_survey_sbc_withnorth.rds")) 
dat_trim <- readRDS(here::here("data", "cleaned_all_bridge_dat_apr_2024.rds"))
#dat_trim <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))
glimpse(dat_trim)

# extracting bridge data from juvenile salmon catch data set
bridge <- dat_trim %>%
  filter(species == "chinook") %>%
  # select(unique_event:utm_y) %>%
  # unique() %>%
  as_tibble() #%>%
# # select(-ends_with("cycle")) %>%
# mutate(utm_x_1000 = utm_x / 1000,
#        utm_y_1000 = utm_y / 1000)

glimpse(bridge)


# reading GSI data (import has data from 1998 to 2021)
gsi <- read_xlsx("data-raw/CHINOOK_event_catch_top_gsi_20240212.xlsx") %>%
  janitor::clean_names()

gsi_mu <- read_xlsx("data-raw/CHINOOK_event_catch_top_gsi_20240212.xlsx",
                            sheet = "catch_top_mu") %>%
  janitor::clean_names() %>%
  mutate(unique_event = gsub("-124J", "", unique_catch))
         #year = as.numeric(sub("\\D*(\\d{4}).*", "\\1", unique_event)))

# same data set but arranged by fish rather than by tow
gsi_by_fish <- read_xlsx("data-raw/CHINOOK_fish_top_gsi_20240212.xlsx",
                         sheet = "top_mu") %>%
  janitor::clean_names() %>%
  mutate(unique_event = gsub("-124[JA]?-.+$", "", unique_fish)) 

gsi_by_fish[grep("124A", gsi_by_fish$unique_fish),] 
# checked this fish and this is a juvenile based on size that was originally
# labeled as an adult

# All percentages for all fish at MU level
gsi_all_mu_in <- read_xlsx(here("data-raw", "CHINOOK_fish_all_gsi_20240424.xlsx"), sheet = "all_mu") %>%
  janitor::clean_names()
glimpse(gsi_all_mu_in)

# All percentages for all fish at MU level
gsi_all_stock_in <- read_xlsx(here("data-raw", "CHINOOK_fish_all_gsi_20240424.xlsx"), sheet = "all_stock") %>%
  janitor::clean_names()
glimpse(gsi_all_stock_in)

# Key for unique_event
event_key <- read_xlsx(here("data-raw", "CHINOOK_fish_all_gsi_20240424.xlsx"), sheet = "fish_catch_event") %>%
  janitor::clean_names()

# Key for stock to Cam's regions

#stock_key <-readRDS("data-raw/finalStockList_Mar2024.rds")
#stock_key <-read_csv("data-raw/finalStockList_Mar2024_updatedSeb.csv")
stock_key <-read_csv("data-raw/finalStockList_Apr2024.csv")

filter(stock_key, grepl("^LEWIS", stock))


gsi_all <- left_join(gsi_all_mu_in, event_key, by = "unique_fish") 

gsi_all_stock <- left_join(gsi_all_stock_in, event_key, by = "unique_fish") %>%
  rename(stock = stock_final) %>%
  mutate(stock = if_else(stock == "LEWIS_ LAKE", "LEWIS_LAKE", stock)) %>% # Renaming misspelled stock
  left_join(select(stock_key, stock, Region3Name), by = "stock") %>% # SELECTING REGION3 HERE
  select(unique_event, unique_fish, stock, region = Region3Name, stock_prob)  

table(gsi_all_stock$region)

nas <- gsi_all[is.na(gsi_all$mu_name),]
nas2 <- gsi_all_stock[is.na(gsi_all_stock$region), ]

unique(nas2$stock)

table(nas2$stock)


table(gsi_all$mu_name)
levels(factor(gsi_all$mu_name))


table(gsi_all_stock$region)
levels(factor(gsi_all_stock$region))

##########

gsi_all_stock_grouped <- gsi_all_stock %>% filter(!is.na(region)) %>% 
  select(-stock) %>%
  mutate(region = fct_other(region, drop = "Fraser River", other_level = "SOG")) %>% # Combining Fraser into SOG
  complete(region, nesting(unique_event, unique_fish)) %>%
  mutate(stock_prop = if_else(is.na(stock_prob), 0, stock_prob/100)) %>% # converting NAs to zeros
  group_by(unique_event, region) %>% # adding probabilities for across all fish in each unique_event
  summarize(stock_prop = sum(stock_prop), n_fish = length(unique(unique_fish))) %>%
  ungroup() %>%
  arrange(unique_event, region) 


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


filter(gsi_all_stock, unique_event %in% lowratios) %>% print(n = 90)
filter(gsi_all_stock, unique_event %in% highratios) %>% print(n = 90)

filter(gsi_all_stock)

filter(gsi_all_stock, unique_event == "BCSI-201731-HW01")
filter(gsi_all, unique_event == "BCSI-201731-HW01")


chinook_gsi_counts <- inner_join(gsi_all_stock_grouped, bridge, by = "unique_event") %>%
  mutate(
    region = as.factor(region),
    year_f = as.factor(year),
    month_f = month(date, label = TRUE),
    month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night))
glimpse(chinook_gsi_counts)

saveRDS(chinook_gsi_counts, file = here("data", "chinook_gsi_counts.rds"))



gsi_all %>% filter(!is.na(mu_name)) %>% 
  group_by(mu_name) %>% summarize(count = n()) %>% print(n = 22)

gsi_all %>% filter(!is.na(mu_name)) %>% 
mutate(mu_name = fct_lump_min(mu_name, min = 400)) %>% 
  group_by(mu_name) %>% summarize(count = n()) %>% print(n = 22)

gsi_all_grouped <- gsi_all %>% filter(!is.na(mu_name)) %>% 
  mutate(mu_name = fct_lump_min(mu_name, min = 400)) %>% select(unique_event, unique_fish, starts_with("mu")) 

gsi_all_grouped %>% 
  group_by(unique_event, unique_fish) %>% 
  pivot_wider(id_cols = -unique_event, names_from = mu_name, values_from = mu_prob)

gsi_all_grouped %>% complete(mu_name, nesting(unique_event, unique_fish)) %>% pull(mu_name) %>% table()
gsi_all_grouped %>% group_by(unique_event) %>%
  complete(mu_name, unique_fish) %>% pull(mu_name) %>% table()

levels(gsi_all_grouped$mu_name)

filter(gsi_all_grouped, unique_event == "BCSI-201731-HW03")
filter(gsi_all_grouped, unique_event == "BCSI-201731-HW15")

gsi_all_expanded <- gsi_all_grouped %>% 
  group_by(unique_event, unique_fish, mu_name) %>%
  summarise(mu_prob = sum(mu_prob)) %>% # adding "Other" mu_names as there can be multiple per unique_fish
  arrange(mu_name, .by_group = TRUE) %>% 
  ungroup() %>%
  arrange(unique_event, unique_fish, mu_name) %>%
  complete(mu_name, nesting(unique_event, unique_fish)) %>% # expanding grid 
  group_by(unique_event, unique_fish) %>%
  mutate(mu_prob = if_else(is.na(mu_prob), 0, mu_prob)) %>% # converting NAs to zeros
  group_by(unique_event, mu_name) %>% # adding probabilities for across all fish in each unique_event
  summarize(mu_prop = sum(mu_prob)/100, n_fish = n()) %>% 
  ungroup() %>%
  arrange(unique_event, mu_name) 

gsi_all_expanded %>% print(n =26)

table(gsi_all_expanded$mu_name)

gsi_all_expanded %>% ungroup() %>% group_by(unique_event) %>%
  summarise(propsum = sum(mu_prop), propratio = sum(mu_prop)/unique(n_fish)) %>%
  pull(propratio) %>% hist(breaks = seq(0.5,1.1, by=0.01))


# chinook_gsi_counts <-inner_join(gsi_all_expanded, bridge, by = "unique_event") %>%
#   mutate(
#     year_f = as.factor(year),
#     month_f = month(date, label = TRUE),
#     month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
#     yday = lubridate::yday(date),
#     utm_x_1000 = utm_x / 1000,
#     utm_y_1000 = utm_y / 1000,
#     effort = log(volume_km3),
#     scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
#     scale_depth = scale(as.numeric(target_depth))[ , 1],
#     year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
#     day_night = as.factor(day_night))
# glimpse(chinook_gsi_counts)
# 
# saveRDS(chinook_gsi_counts, file = here("data", "chinook_gsi_counts.rds"))


hist(chinook_gsi_counts$scale_dist)
chinook_gsi_counts%>% pull(mu_name) %>% table()

#%>% #%>%
  filter(!is.na(mu_pl_name) & !is.na(lat))
glimpse(gsi_trim)




##########################

head(gsi$unique_catch)
head(gsi_mu$unique_catch)
head(gsi_by_fish$unique_fish)

#find which strings in vector1 are NOT in vector2
gsi_mu$unique_event[!(gsi_mu$unique_event %in% bridge$unique_event)]# lots!
base::setdiff(unique(dat_trim$unique_event), unique(gsi_mu$unique_event))

# this GSI tow is not in the catch data
filter(gsi_mu, unique_event == "HS201637-HW08") 
filter(dat_trim, unique_event == "HS201637-HW08")
filter(bridge, unique_event == "HS201637-HW08")


dat_trim[grep("^HS20", bridge$unique_event),"unique_event"]
#dat_in[grep("^HS2016", dat_in$unique_event),"unique_event"]

gsi_trim <- left_join(gsi_mu, bridge, by = "unique_event") %>% #%>%
 filter(!is.na(mu_pl_name) & !is.na(lat))
glimpse(gsi_trim)


saveRDS(gsi_trim, here::here("data", "gsi-trim.rds"))

gsi_trim_na <- filter(gsi_trim, is.na(mu_pl_name))





gsi_fish_trim <- left_join(gsi_by_fish, bridge, by = "unique_event") %>%
 # drop_na() %>%
  mutate(MU = fct_lump_min(mu_pl_name, min = 100),
         season_f = case_when(
           month %in% c("2", "3", "4") ~ "spring",
           month %in% c("5", "6", "7", "8") ~ "summer",
           month %in% c("9", "10", "11", "12") ~ "winter"
         ),
         month_f =  month(month, label = TRUE)
        )

# Range of years 
range(gsi_fish_trim$year)

# looking at the GSI tows that aren't in list of salmon catch tows
gsi_out <- left_join(gsi_mu, bridge, by = "unique_event")  %>%
  filter(is.na(year)) %>% 
  mutate(year2 = as.numeric(sub("\\D*(\\d{4}).*", "\\1", unique_event))) %>%
  select(year2, everything()) 

gsi_out

gsi_out %>%
  group_by(year2) %>%
  count(unique_event) %>%
  print(n = 700)
  

gsi_trim %>%
  group_by(unique_event) %>%
  summarise(n = n())


#(gsi_mu$unique_event %in% dat_in$unique_event)


dat_trim %>%
  filter(str_detect(unique_event, "IPES2019-124"))


# do
gsi_mu %>% group_by(unique_catch) %>%
summarise(n = n()) %>% pull(n) %>% hist()


gsi_mu %>% group_by(mu_pl_name) %>%
  summarise(n = n()) %>%
  print(n = 23)

gsi_mu %>% filter(is.na(mu_pl_name))

gsi_mu %>%
  group_by(unique_event) %>%
  summarise(n = sum(mu_n)) %>%
  pull(n) %>% hist(breaks = 0:55, main = "Distribution of GSI samples per tow")

filter(gsi_trim) %>%
ggplot(data = .) +
  geom_col(aes(y = mu_n, x = unique_event, fill = mu_pl_name), position = "stack") +
  coord_flip() +
  facet_grid(year~month, scales = "free_y")

filter(gsi_trim) %>%
  ggplot(data = .) +
  geom_col(aes(y = mu_n, x = week, fill = mu_pl_name), position = "stack") +
  coord_flip() +
  facet_wrap(~month, scales = "free_y")


filter(gsi_trim) %>%
  ggplot(data = .) +
  geom_col(aes(y = mu_n, x = week, fill = mu_pl_name), position = "stack") +
  coord_flip() +
  facet_wrap(~month, scales = "free_y")


filter(gsi_trim) %>%
  ggplot(data = .) +
  geom_col(aes(y = mu_n, x = week, fill = mu_pl_name), position = "stack") 

gsi_trim %>%
  group_by(month, mu_pl_name) %>%
  summarise(n = sum(mu_n)) %>%
  print(n = 90)

# plots by fish as this way MUs with few samples can more easily
# be lumped into "Other" category

table(gsi_fish_trim$mu_pl_name)

gsi_fish_trim %>%
  group_by(mu_pl_name) %>%
  count()

gsi_fish_trim %>%
  count(mu_pl_name)

gsi_fish_trim %>%
  count(unique_event) %>%
  ggplot(.) + geom_histogram(aes(n), binwidth = 1) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(breaks = seq(0, 55, by = 5)) +
  ggtitle("Distribution of number of fish in GSI samples by tow in catch dataset")
ggsave(here("figs", "gsi", "hist-gsi-samples-by-tow.png"), width = 9, height = 5)

gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  mutate(MU = fct_lump_min(mu_pl_name, min = 100)) %>%
  ggplot(data = .) +
  geom_bar(aes(x = week, fill = MU), position = "stack") +
  #coord_flip() +
  #facet_wrap(~month, scales = "free_y") +
  theme_sleek() +
  NULL
ggsave(here("figs", "gsi", "gsi-chinook-by-week.png"), width = 9, height = 5)

month_labels <- c('J','F','M','A','M','J','J','A','S','O','N','D')

# faceting by year
gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  mutate(MU = fct_lump_min(mu_pl_name, min = 100)) %>%
  mutate(month_f =  month(month, label = TRUE)) %>%
  ggplot(data = .) +
  geom_bar(aes(x = month_f, fill = MU), position = "stack") +
  #coord_flip() +
  facet_wrap(~year, scales = "free_y") +
  scale_x_discrete(labels = month_labels, drop = FALSE) +
  theme_sleek() +
  NULL
ggsave(here("figs", "gsi", "gsi-chinook-by-month-year.png"), width = 9, height = 8)



# faceting by year but with fixed y scale
gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  mutate(MU = fct_lump_min(mu_pl_name, min = 100)) %>%
  mutate(month_f =  month(month, label = TRUE)) %>%
  ggplot(data = .) +
  geom_bar(aes(x = month_f, fill = MU), position = "stack") +
  #coord_flip() +
  facet_wrap(~year, scales = "fixed") +
  scale_x_discrete(labels = month_labels, drop = FALSE) +

  theme_sleek() +
  NULL
ggsave(here("figs", "gsi", "gsi-chinook-by-month-year-fixed-y.png"), width = 9, height = 8)


# Need to figure out method for estimating TOTAL number of fish in each sample 
# size bin

gsi_fish_trim %>%
  count(unique_event) %>%
  group_by(n) %>% count(name = "ncount") %>%
  ungroup() %>%
  mutate(ntotal = n * ncount,
         cumcount = cumsum(ntotal),
         cumprop = cumcount/sum(ntotal)) %>%
  ggplot(.) + geom_col(aes(x = n, y = ntotal)) +
  geom_line(aes(x = n, y = cumprop*1300)) + geom_point(aes(x = n, y = cumprop*1300)) +
  scale_x_continuous(breaks = seq(0, 55, by = 5)) +
  scale_y_continuous("Fish count per number of samples", breaks = seq(0, 1400, by = 200),
                     sec.axis = sec_axis(~./1300, name = "Cumulative proportion")) +
  ggtitle("Distribution of number of fish in GSI samples by tow in catch dataset") +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "hist-gsi-total-fish-by-tow.png"), width = 9, height = 5)


# Plotting locations of individual data points by month color coded by MU
source(here("R", "plot_map.R"))

gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = MU)) +
  xlim(range(gsi_fish_trim$utm_x)) +
  ylim(range(gsi_fish_trim$utm_y)) +
  facet_wrap(~season_f, nrow = 2) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "gsi-locations-by-season.png"), width = 9, height = 6)

gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = MU)) +
  xlim(range(gsi_fish_trim$utm_x)) +
  ylim(range(gsi_fish_trim$utm_y)) +
  xlab("")  + ylab("") +
  facet_wrap(~month_f, drop = FALSE) +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "gsi-locations-by-month.png"), width = 9, height = 6)


gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  full_join(crossing(season_f = unique(gsi_fish_trim$season_f), year=unique(gsi_fish_trim$year))) %>% 
  mutate(empty=ifelse(is.na(MU), "No samples", NA_character_),
         x=mean(range(utm_x, na.rm=TRUE)), 
         y=mean(range(utm_y, na.rm=TRUE))) %>%
  ungroup() %>%
  select(-month_f) %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = MU)) +
  facet_grid(season_f~year) + # , margins = "year") + # margins doesn't work with geom_sf()
  geom_text(aes(x, y, label=empty), colour="red", size=3) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  xlim(range(gsi_fish_trim$utm_x)) +
  ylim(range(gsi_fish_trim$utm_y)) +
  ggsidekick::theme_sleek()
ggsave(here("figs", "gsi", "gsi-locations-season-by-year.png"), width = 40, height = 5)


library(gridExtra)

png(here("figs", "gsi", "gsi-locations-season-by-year-2.png"), width = 20, height = 9, units = "in", res = 200)

grid.arrange(

gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  full_join(crossing(season_f = unique(gsi_fish_trim$season_f), year=unique(gsi_fish_trim$year))) %>% 
  mutate(empty=ifelse(is.na(MU), "No samples", NA_character_),
         x=mean(range(utm_x, na.rm=TRUE)), 
         y=mean(range(utm_y, na.rm=TRUE))) %>%
  ungroup() %>%
  filter(year %in% 1998:2009) %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = MU)) +
  xlim(range(gsi_fish_trim$utm_x)) +
  ylim(range(gsi_fish_trim$utm_y)) +
  facet_grid(season_f~year) +
  geom_text(aes(x, y, label=empty), colour="red", size=3) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek(),

gsi_fish_trim %>%
  filter(mu_pl_name != "Below_Limit") %>%
  full_join(crossing(season_f = unique(gsi_fish_trim$season_f), year=unique(gsi_fish_trim$year))) %>% 
  mutate(empty=ifelse(is.na(MU), "No samples", NA_character_),
         x=mean(range(utm_x, na.rm=TRUE)), 
         y=mean(range(utm_y, na.rm=TRUE))) %>%
  ungroup() %>%
  filter(year %in% 2010:2021) %>%
  #mutate(MU = fct_lump_n(mu_pl_name, n = 9)) %>%
  ggplot(data = .) +
  geom_sf(data = all_coast, color = "gray80", fill = NA) +
  geom_point(aes(x = utm_x, y = utm_y, color = MU)) +
  xlim(range(gsi_fish_trim$utm_x)) +
  ylim(range(gsi_fish_trim$utm_y)) +
  facet_grid(season_f~year) +
  geom_text(aes(x, y, label=empty), colour="red", size=3) +
  xlab("")  + ylab("") +
  #facet_grid(month_f ~ year) +
  ggsidekick::theme_sleek(),
nrow=2)

dev.off()

