### Data wrangling prior to model fitting

library(tidyverse)
library(sf)
library(sp)
library(hms)


# add lat/lon
get_ll <- function(x, y, zone, loc){
  points_utm <- SpatialPoints(
    cbind(x, y), 
    proj4string = CRS(paste0("+proj=utm +zone=",zone[1]," +units=m"))
  )
  points_ll <- spTransform(
    points_utm, 
    CRS("+proj=longlat +datum=WGS84")
  )
  
  
  if (loc == "x") {
    return(coordinates(points_ll)[,1])
  } else if (loc == "y") {
    return(coordinates(points_ll)[,2])
  }
}

## add UTM to bridge

## This code uses the start longitude to get the UTM zone. 
##  NOTE that for some modelling consistent zone of 9 more appropriate 

## use zone = 9 for bc coastwide modelling as per Sean Anderfson suggestion. 
##  use zone = "vary" if you want the correct utm zone for the longitude to be used
##  note that data transcribed from multiple zone transformation should not be combined

# get_utm <- function(x, y, zone, loc){
#   
#   if (zone == "9") {
#     utm_zone <- zone
#   } else if (zone == "vary") {
#     utm_zone <- floor((x/6)+31)
#   }
#   
#   epsg = paste0("+init=epsg:326",(str_pad(utm_zone, 2,
#                                           side = "left", pad = "0")))
#   
#   points = SpatialPoints(cbind(x, y),
#                          proj4string = CRS("+proj=longlat +datum=WGS84"))
#   
#   points_utm = spTransform(points, CRSobj = CRS(epsg))
#   
#   if (loc == "x") {
#     return(coordinates(points_utm)[,1])
#   } else if (loc == "y") {
#     return(coordinates(points_utm)[,2])
#   } else if (loc == "z") {
#     return(utm_zone)
#   }
# }

get_utm <- function(x, y, zone, loc){
  points = SpatialPoints(cbind(x, y),
                         proj4string = CRS("+proj=longlat +datum=WGS84"))
  points_utm = spTransform(
    points, CRS(paste0("+proj=utm +zone=",zone[1]," +units=m"))
  )
  if (loc == "x") {
    return(coordinates(points_utm)[,1])
  } else if (loc == "y") {
    return(coordinates(points_utm)[,2])
  }
}

# testutm <- test1 %>% 
#   mutate(
#     utm_x = get_utm(START_LONGITUDE, START_LATITUDE, zone = 9, loc = "x"),
#     utm_y = get_utm(START_LONGITUDE, START_LATITUDE, zone = 9, loc = "y"),
#     utm_zone = get_utm(START_LONGITUDE, START_LATITUDE, zone = 9, loc = "z")
#   ) 

### Data pull done by myself using Amy's code and local db version BCSI_be_20240918
## Nov 20 2024 update
## This data set does not have missing species in tows with zero salmon
## so no need to wrangle as below
dat_in <- read.csv(here::here("data-raw", "Seb_SalmonCounts_all_20250121.csv"),
                    stringsAsFactors = FALSE) %>% 
  janitor::clean_names() %>% 
  mutate(
         utm_x = get_utm(lon, lat, zone = 9, loc = "x"),
         utm_y = get_utm(lon, lat, zone = 9, loc = "y"),
  )
glimpse(dat_in)



## MISSING EVENTS WITH ZERO TOTAL CATCHES --------------------------------------

dat_in %>% group_by(unique_event) %>%
  summarise(n = n()) %>%
  pull(n) %>% range()

# catch data excludes stations w/ 0 catches for all species; correct and check
dat_in %>% group_by(unique_event) %>%
  summarise(n = sum(n_total)) %>% 
  filter(n == 0)
# if no rows returned then there are no unique_events with zero salmon catches,
# these are being likely shown as NAs

dat_in %>% group_by(unique_event) %>%
  summarise(n = sum(n_total)) %>% 
  filter(is.na(n)) %>% dim()

filter(dat_in, scientific_name == "") %>% dim() # same as above

missing_catches <- dat_in %>% 
  filter(is.na(n_total))

dim(missing_catches) # must be same row number as is.na(n) test above

sp_vec <- filter(dat_in, scientific_name != "") %>% 
  pull(scientific_name) %>% unique()

catch_list <- vector(mode = "list", length = length(sp_vec))

for (i in seq_along(catch_list)) {
  catch_list[[i]] <- missing_catches %>% 
    mutate(
      n_juv = 0,
      n_ad = 0,
      n_total = 0,
      scientific_name = sp_vec[i]
    )
} 

catch_list_flat <- bind_rows(catch_list) %>% 
  arrange(unique_event, scientific_name)

catch_out <- dat_in %>% 
  filter(!is.na(n_total)) %>% 
  rbind(., bind_rows(catch_list))

identical(catch_out, dat_in)

# This illustrates the before/after of the checklists with zero catches
dat_in %>% filter(unique_event == "BCSI-201778-QCSD01")
catch_out %>% filter(unique_event == "BCSI-201778-QCSD01")

dim(dat_in)
dim(catch_out)
#dat_in <- catch_out

### checking for samples with empty day_night field 

dn <- filter(dat_in, day_night == "")
dn2023 <- filter(dat_in, year == "2023")

# creating new column with new day_night value based on time of date field
# however this can be problematic as sometimes there can be a lag between
# the start of a unique_event and the time fishing begins.
# Need tow start or end time instead.
# Also, can use the oce package to quickly calculate sunset/sunrise based on 
# local time and lat/lon using sunAngle()
dat_in %>%
  mutate(time = as_hms(ymd_hms(date)),
         diff = Mod(difftime(time, as_hms("13:00:00"), units = "hours")),
         dn = case_when(
           diff > 8 & time != as_hms("00:00:00") ~ "NIGHT",
           diff <= 8 & time != as_hms("00:00:00") ~ "DAY",
           time == as_hms("00:00:00") ~ "notime")
  ) %>%
  filter(scientific_name == "ONCORHYNCHUS TSHAWYTSCHA" & dn != day_night & dn != "notime") %>%
  select(unique_event, date, time, day_night, dn) 

## FLAG IPES -------------------------------------------------------------------

## import and transform IPES grid
# ipes_sf_deg <- 
ipes_sf_poly <- st_read(
  here::here("data-raw", "spatial", "ipes_shapefiles", "IPES_Grid_UTM9.shp")) %>% 
  #convert into single polygon
  st_union(., by_feature = FALSE)


dat_sf <- dat_in %>% 
  select(unique_event, utm_y, utm_x) %>% 
  st_as_sf(., coords = c("utm_x", "utm_y"), 
           crs = st_crs(ipes_sf_poly))

# check
# ggplot() +
#   geom_sf(data = st_intersection(dat_sf, ipes_sf_poly))

# extract events within ipes survey grid
ipes_grid_events <- st_intersection(dat_sf, ipes_sf_poly) %>% 
  pull(., unique_event)

##------------------------------------------------------------------------------
# make year key representing pink/sockeye cycle lines
yr_key <- data.frame(
  year = unique(dat_in$year)
) %>% 
  arrange(year) %>% 
  mutate(
    sox_cycle = rep(1:4, length.out = length(unique(dat_in$year))) %>% 
      as.factor(),
    pink_cycle = rep(1:2, length.out = length(unique(dat_in$year))) %>% 
      as.factor()
  )



## MISSING BATHY DATA ----------------------------------------------------------

missing_bathy <- dat_in # %>%
 # no depth_mean_m column at the moment as don't have that db table
 # filter(is.na(depth_mean_m) | depth_mean_m < 0) 

dim(missing_bathy)

# Check for local NOAA bathymetric grid file before downloading
bathy_grid_filename <- paste("NOAA_bathy_grid",
                             round(min(missing_bathy$lat), 3),
                             round(max(missing_bathy$lat), 3),
                             round(min(missing_bathy$lon), 3),
                             round(max(missing_bathy$lon), 3), 
                             "res_0.1.rds", sep = "_")

if (file.exists(here::here("data-raw", "spatial", bathy_grid_filename))) {
  bathy_grid <- readRDS(file = here::here("data-raw", "spatial",
                                          bathy_grid_filename))
} else {
  
  bathy_grid <- marmap::getNOAA.bathy(
    lon1 = min(missing_bathy$lon),
    lon2 = max(missing_bathy$lon),
    lat1 = min(missing_bathy$lat),
    lat2 = max(missing_bathy$lat),
    resolution = 0.1
  )
  
  saveRDS(bathy_grid, 
          file = here::here("data-raw", "spatial", 
                            bathy_grid_filename))
}

bathy_matrix <- matrix(data = NA, 
                       nrow = nrow(bathy_grid), ncol = ncol(bathy_grid))

for (i in 1:nrow(bathy_matrix)) {
  bathy_matrix[i, ] <- bathy_grid[i, ] %>% as.numeric()
}

dimnames(bathy_matrix) <- dimnames(bathy_grid)
bathy_df <- data.frame(
  lon = rownames(bathy_matrix),
  bathy_matrix) %>% 
  pivot_longer(., cols = -(lon), names_to = "lat",
               names_prefix = "X", values_to = "depth") %>%
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat)) 

# subset to remove missing depths that are on land or inshore of 20 m isobath
bathy_df_trim <- bathy_df %>%
  filter(!is.na(depth),
         !depth > -20) 

# identify closest depth data
tt <- hutilscpp::match_nrst_haversine(
  lat = missing_bathy$lat,
  lon = missing_bathy$lon,
  addresses_lat = as.numeric(bathy_df_trim$lat),
  addresses_lon = as.numeric(bathy_df_trim$lon)
)
missing_bathy$pred_depth_mean_m <- -1 * bathy_df_trim$depth[tt$pos]
hist(missing_bathy$pred_depth_mean_m)

## MISSING BRIDGE DATA ---------------------------------------------------------

imp_dat <- dat_in %>% 
  mutate(
    # replace 0s
    volume_km3 = ifelse(volume_km3 == "0", NaN, volume_km3),
   # dist_to_coast_km = ifelse(dist_to_coast_km == "0", NaN, dist_to_coast_km), # No dist_to_coast_km
   NULL
  ) 



# small number so probably fine but should ultimately redo distance coast calc
# imp_dist <- imp_dat %>% 
#   select(unique_event, utm_x, utm_y, dist_to_coast_km) %>%
#   distinct() %>% 
#   VIM::kNN(.) 

# NOTE THIS SHOULD BE REPEATED WITH MORE BRIDGE VARIABLES
imp_eff <- imp_dat %>% 
  select(unique_event, year, week, utm_x, utm_y, target_depth,  
         volume_km3) %>%
  distinct() %>% 
  VIM::kNN(.) 




## REJOIN ----------------------------------------------------------------------

dat <- dat_in %>%
  #remove data imputed above
  select(-c(utm_x, utm_y, year, week, target_depth, volume_km3#,  dist_to_coast_km
            )) %>%
  # add missing bathy data
  left_join(
    ., 
    missing_bathy %>% select(unique_event, scientific_name, pred_depth_mean_m), 
    by = c("unique_event", "scientific_name")
  ) %>% 
  # add imputed effort data
  left_join(
    ., 
    imp_eff %>% select(-ends_with("imp")),
    by = "unique_event"
  ) %>%
  # add imputed distance data
  # left_join(
  #   ., 
  #   imp_dist %>% select(-ends_with("imp"), -utm_x, -utm_y),
  #   by = "unique_event"
  # ) %>%
  # add cycle line IDs
  left_join(
    ., 
    yr_key,
    by = "year"
  ) %>%
  mutate(
    # use predicted depth if actual is missing
    # or in this case always use predicted depth as no depth available from db
    # depth_mean_m = ifelse(
    #   is.na(pred_depth_mean_m), depth_mean_m, pred_depth_mean_m
    # ),
    depth_mean_m = pred_depth_mean_m,
    #add common names
    species = case_when(
      grepl("GORBUSCHA", scientific_name) ~ "pink",
      grepl("KETA", scientific_name) ~ "chum",
      grepl("KISUTCH", scientific_name) ~ "coho",
      grepl("NERKA", scientific_name) ~ "sockeye",
      grepl("TSHAWYTSCHA", scientific_name) ~ "chinook",
      TRUE ~ "other"
    ),
    # adjust season
    season_f = case_when(
      month %in% c("2", "3", "4") ~ "sp",
      month %in% c("5", "6", "7", "8") ~ "su",
      month %in% c("9", "10", "11", "12") ~ "wi"
    ) %>% 
      as.factor(.),
    # define core area as square around BC (with a little of AK and WA)
    synoptic_station = ifelse(
      lat > 47 & lat < 57 &  # #!grepl("GS", unique_event) &  
        # Originally cutting at 56 N but if we want to include Sumner St in AK then changing to 57 N
        lon > -136, 
      TRUE,
      FALSE
    ),
    # define IPES based on intersections with grid
    ipes_grid = ifelse(
      unique_event %in% ipes_grid_events,
      TRUE,
      FALSE
    ),
    survey_f = ifelse(
      year > 2016 & season_f == "su" & ipes_grid == TRUE,
      "ipes", 
      "hss"
    ) %>% 
      as.factor(),
    year_f = as.factor(year),
    yday = lubridate::yday(date)
  ) 


# subset to core area and remove rows with missing data
dat_trim <- dat %>%
  filter(
    # remove steelhead and mystery species
    !species == "other",
    # exclude outside of box around BC
    synoptic_station == TRUE,
    # previously excluding years prior to 1995 due to small sample sizes,
    # however these older surveys cover less typical months, which is better
    # for modelling seasonal abundances (atlas) rather than annual abundances (index).
    # Importantly, 1995 has data for April, which no other year has
    # year >= 1998,
    # exclude SoG and Puget Sound PFMAs
    # !pfma %in% c("13", "14", "15", "16", "17", "18", "19","28", "29", # This line is all SoG PFMAs
    #              "AKPEN", "KI", "OR", "PS", "PWS", "SCA"), # PS = Puget Sound
    # OFF (offshore), SEA(southeast Alaska) and ISEA (Inner southeast Alaska)
    # have points that are within synoptic station that are useful for modelling 
    # So they aren't excluded.
    # remove one station clearly on land (CHECK FUTURE DATA PULLS)
    
    # !pfma %in% c("13", "14", "15", "16", "17", "18", "19","28", "29",
    #              "AKPEN", "ISEA", "KI", "OFF", "OR", "PS", "PWS", "SCA", "SEA") # PS = Puget Sound
  ) %>%
  select(unique_event, date, year, year_f, month, week, day, yday, day_night,
         pink_cycle, sox_cycle, pfma,
         season_f, survey_f, lat, lon, utm_x, utm_y, target_depth, depth_mean_m,# dist_to_coast_km, 
         volume_km3, n_juv, n_ad, n_total, species) %>% 
  droplevels()

# filedate <- sprintf("%04d%02d%02d",year(today()),month(today()),day(today()))
# ymd(Sys.Date())
# today()

saveRDS(dat_trim, here("data", "cleaned_atlas_catch_dat.rds"))
saveRDS(dat, here("data", "cleaned_all_bridge_catch_dat.rds"))

dim(dat_trim)

## MISC CHECKS -----------------------------------------------------------------

# dat_trim <- readRDS(here("data", "cleaned_atlas_catch_dat.rds"))
# dat <- readRDS(here("data", "cleaned_all_bridge_catch_dat.rds"))

min_lat <- min(floor(dat_trim$lat) - 0.1)
max_lat <- max(dat_trim$lat) + 0.1
min_lon <- min(floor(dat_trim$lon) - 0.1)
max_lon <- max(dat_trim$lon) + 0.1

bc_coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = min_lon-2, ymin = min_lat-1, xmax = max_lon+2, ymax = max_lat+1) %>%
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

dt_chinook <- dat_trim %>% filter(species == "chinook")

all_coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                            returnclass = "sf"), 
                  rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = -163, ymin = 40, xmax = max_lon+2, ymax = 61) %>%
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

dat %>%
  group_by(pfma) %>% tally() %>% print(n = 70)

dat_trim %>%
  group_by(pfma) %>% tally() %>% print(n = 70)

# All points (inc outside synoptic) from pfmas with letters
dat %>% filter(!grepl("[0-9]", pfma)) %>% 
  ggplot() +
  geom_sf(data = all_coast) + 
  geom_point(aes(x = lon, y = lat, colour = pfma)) 

# trimmed data map
dat_trim %>% filter(species == "chinook") %>%
  ggplot() +
  geom_point(aes(x = lon, y = lat, colour = pfma)) +
  geom_sf(data = bc_coast) +
  ggsidekick::theme_sleek() +
  ylim(c(min_lat, max_lat)) + 
  xlim(c(min_lon, max_lon)) 

ggsave(here::here("figs", "dataset", "trimmed-data-by-pfma-nov-2024.png"), height = 9, width = 9)

dt_chinook <- dat_trim %>% filter(species == "chinook")

dat_trim %>% filter(species == "chinook") %>%
  filter(pfma %in% c("13", "14", "15", "16", "17", "18", "19","28", "29"))  %>%
  ggplot() +
  geom_point(aes(x = lon, y = lat, colour = pfma)) +
  geom_sf(data = bc_coast) +
  ggsidekick::theme_sleek() +
  ylim(c(48, 51)) + 
  xlim(c(-127, -122)) +
  facet_grid(month~year)
         
# Only missing data for January
table(dt_chinook$year, dt_chinook$month)
table(dt_chinook$month)


filter(dat_trim, pfma == "IBC")

table(dat_trim$month) # observations per month 
table(dat$month)

ggplot(dat_trim %>% filter(species == "sockeye")) +
  geom_point(aes(x = lon, y = lat, colour = survey_f)) +
  geom_sf(data = bc_coast) +
  #facet_wrap(year_f ~ season_f) +
  xlim(c(min_lon, max_lon)) +
  ylim(c(min_lat, max_lat)) +
  ggsidekick::theme_sleek()
ggsave(here::here("figs", "dataset", "trimmed-data-by-survey-apr-2024.png"), height = 9, width = 8)

dat_trim$season_f2 <- fct_recode(
  dat_trim$season_f, "spring" = "sp", "summer" = "su", "fall" = "wi"
)

shape_pal <- c(21, 22, 23)
names(shape_pal) <- unique(dat_trim$season_f2)



# bubble plots of temporal coverage
dt_chinook %>% 
  group_by(year, week, survey_f, season_f2) %>% 
  summarize(n_tows = length(unique_event), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(y = week, x = year, size = n_tows, fill = survey_f, 
                  shape = season_f2),
              alpha = 0.3, width = 0.25) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number\nof\nTows") +
  scale_fill_discrete(name = "Survey") +
  scale_shape_manual(values = shape_pal, name = "Season") +
  labs(x = "Year", y = "Week") +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  )
ggsave(here::here("figs", "dataset", "temporal-coverage-by-year-week-apr-2024.png"), height = 6, width = 10)

dt_chinook %>% 
  mutate(month_f = month(month, label = TRUE)) %>%
  group_by(year, month_f, survey_f) %>% 
  summarize(n_tows = length(unique_event), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_point(aes(y = month_f, x = year, size = n_tows, color = survey_f),
              alpha = 0.5) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number\nof\nTows") +
  scale_fill_discrete(name = "Survey") +
  labs(x = "Year", y = "Month") 
ggsave(here::here("figs", "dataset", "temporal-coverage-by-year-month-apr-2024.png"), height = 4, width = 9)

dt_chinook$CommonDate <- as_date(dt_chinook$date)
year(dt_chinook$CommonDate) <- 2000
dt_chinook$date2 <- as_datetime(dt_chinook$date)

dt_chinook %>% 
  mutate(month_f = month(month, label = TRUE)) %>%
  group_by(year, CommonDate) %>% 
  summarize(n_tows = length(unique_event), date2 = first(date2), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_point(aes(y = CommonDate, x =date2, size = n_tows),
             alpha = 0.5) +
  ggsidekick::theme_sleek() +
  scale_size_area(name = "Number of\nTows per day") +
  scale_fill_discrete(name = "Survey") +
  labs(x = "Year", y = "Month") 
ggsave(here::here("figs", "dataset", "temporal-coverage-by-year-day-apr-2024.png"), height = 4, width = 9)


dt_chinook %>% 
  mutate(month_f = month(month, label = TRUE)) %>%
  group_by(year, CommonDate) %>% 
  summarize(n_tows = length(unique_event), date2 = first(date2), .groups = "drop") %>%
  filter(year %in% 2023)



# on dat, not dat_trim, figuring put pfsa = "PS" is Puget Sound
dat %>% filter(species == "sockeye" & pfma == "PS") %>%
  ggplot() +
  geom_point(aes(x = lon, y = lat, colour = survey_f)) +
  geom_sf(data = coast) +
  xlim(c(min_lon, max_lon+1.2)) +
  ylim(c(min_lat, max_lat))

dat %>% filter(species == "sockeye" & synoptic_station != TRUE &
                 !pfma %in% c("13", "14", "15", "16", "17", "18", "19", "PS") &
                 lat < 55 & lat > 46 &lon > -140) %>%

  droplevels() %>%
  ggplot() +
  geom_point(aes(x = lon, y = lat, colour = pfma)) +
  geom_sf(data = coast)



