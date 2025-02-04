# Wrangle tow inf


library(tidyverse)
library(readxl)
library(janitor)
library(here)


# CHANGE TO LOCAL DB PULL
# tow_info <- read_excel(here("data-raw", "Seb_Trawl_Specs_202.xlsx"), 
#                        sheet = "Sheet1") %>%
#   clean_names() 

# CHANGE TO LOCAL DB PULL
tow_info <- read_csv(here("data-raw", "BCSI_TowInfo_20250121.csv")) %>%
  clean_names() 

tow_info$date %>% range()

filter(tow_info, unique_event == "IPES2024-026-126") %>% glimpse


# Fixing vessel names
tow_info$vessel <- forcats::fct_recode(
  tow_info$vessel_name,
  "Sea Crest" = "CFV SEACREST",
  "CCGS Sir John Franklin" = "CCGS SIR JOHN FRANKLIN",
  "Nordic Pearl" = "Nordic Pearl",
  "FV Anita J" = "FV ANITA J",
  "CCGS W.E. Ricker" = "CCGS W.E. RICKER",
  "CFV Ocean Selector" =   "CFV OCEAN SELECTOR",
  "CFV Frosti" = "CFV FROSTI",
  "CFV Viking Storm" = "CFV VIKING STORM",
  "CFV Viking Storm & CCGS W.E. Ricker" = "CFV VIKING STORM AND CCGS W.E. RICKER",
  "CFV Viking Storm & CCGS W.E. Ricker" = "CFV VIKING STORM & CCGS W.E. RICKER",
  "FV Columbia" =   "FV COLUMBIA",
  "Sea Crest" =   "SEA CREST",
  "Nordic Pearl" =   "NORDIC PEARL",
  "CCGS Sir John Franklin" = "Sir John Franklin"
)

saveRDS(tow_info, file = here("data", "tow_info.rds"))
      
dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))
gsidat <- readRDS(here::here("data", "chinook_gsi_counts_fitted_2025-01-27.rds"))
dat$year %>% range()
gsidat$year %>% range()

# old columns
# t_info <- select(tow_info,
#                  unique_event:trip_id, vessel, 
#                  net_desc:station_name,
#                  start_bottom_depth:end_bottom_depth,
#                  mouth_height_use:mouth_width_source)

# new columns
t_info <- select(tow_info,
                 unique_event, vessel,
                 net_desc,
                # start_bottom_depth:end_bottom_depth,
                 mouth_height_use:mouth_width_source)

tow_data <- left_join(dat, t_info, by = "unique_event")
gsitow_data <- left_join(gsidat, t_info, by = "unique_event")

# renaming single NA and also typo (1142 instead of 7742)
td <- tow_data %>%
  mutate(net_desc = if_else(is.na(net_desc), "LFS 7742", net_desc)) %>%
  mutate(net_desc = fct_recode(net_desc,
    "LFS 7742" = "LFS 1142"
    #"LFS 7742" = `NA`
  )) %>%
  group_by(year, vessel, net_desc) %>%
  summarise(n = n()) %>%
  arrange(year) 

# Checking tows with NAs in net_desc
filter(tow_data, is.na(net_desc))

 

print(td, n = 50)

tdw <- pivot_wider(td, names_from = net_desc, values_from = n) %>% print(n = 50)

td <- tow_data %>%
  group_by(year, vessel, net_desc) %>%
  summarise(n = n()) %>%
  arrange(year)


print(td, n = 45)

tdw <- pivot_wider(td, names_from = net_desc, values_from = n) %>% print(n = 50)


library(knitr)
library(kableExtra)
options(knitr.kable.NA = '')

tow_kb <- kable(tdw, format = "latex", align = "llccccccc",
                             caption = "Yearly vessel and net type  for all tows each year") %>%
  kable_styling(font_size = 6)


writeLines(tow_kb, here("tables","tows.tex"))


gsotow_kb <- kable(tdw, format = "latex", align = "llccccccc",
                caption = "Yearly vessel and net type  for all tows each year") %>%
  kable_styling(font_size = 6)


writeLines(tow_kb, here("tables","tows.tex"))

tow_data %>%
  filter(is.na(vessel))

tow_data %>%
  mutate(date = lubridate::as_date(date)) %>%
  arrange(desc(date))

gsitow_data %>%
  mutate(date = lubridate::as_date(date)) %>% 
  select(date, everything()) %>%
  arrange(desc(date))

