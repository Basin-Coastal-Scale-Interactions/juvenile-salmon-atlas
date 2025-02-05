# Wrangle table of yearly tows counts by ship and net

library(tidyverse)
library(readxl)
library(janitor)
library(here)

# All tows
tow_info <- read_csv(here("data-raw", "BCSI_TowInfo_20250121.csv")) %>%
  clean_names() 

tow_info$date %>% range()

# No day_night
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
gsidat <- readRDS(here::here("data", "chinook_gsi_counts_fitted.rds"))
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

# Checking tows with NAs in net_desc
filter(tow_data, is.na(net_desc))


# Fixing net values here
# jackie.king 7 November, 1:06 pm
# "The Ocean Selector is a groundfish boat, so we would not have used a bottom 
# trawl net.  I would expect the CanTrawl 250 to have been loaded."
td <- tow_data %>%
  mutate(net_desc = case_when(
    year >= 2021 & is.na(net_desc) ~ "LFS 7742",
    net_desc == "Ocean Selectors net" ~ "Ocean Selector's net",
    net_desc == "Rusty's backup net" ~ "CanTrawl 250",
    .default = net_desc
  )) %>%
  group_by(year, vessel, net_desc) %>%
  summarise(n = n()) %>%
  arrange(year)


print(td, n = 48)

tdw <- pivot_wider(td, names_from = net_desc, values_from = n) %>% print(n = 50)


library(knitr)
library(kableExtra)
options(knitr.kable.NA = '')

tow_kb <- kable(tdw, format = "latex", align = "llccccccc",
                  caption = "Yearly vessel and net type for all tows used in this 
                study. Note that CanTrawl 400/580 and LFS 1142 are much larger nets than 
                CanTrawl 250 and LFS 7742, with these latter two being comparable in size.") %>%
  kable_styling(font_size = 6) %>%
  kable_styling(latex_options = "hold_position")


writeLines(tow_kb, here("tables","tows.tex"))

# gsotow_kb <- kable(tdw, format = "latex", align = "llccccccc",
#                 caption = "Yearly vessel and net type  for all tows each year") %>%
#   kable_styling(font_size = 6)

print(tdw, n = 44)

tow_data %>%
  filter(is.na(vessel))

tow_data %>%
  mutate(date = lubridate::as_date(date)) %>%
  arrange(desc(date))

# gsitow_data %>%
#   mutate(date = lubridate::as_date(date)) %>% 
#   select(date, everything()) %>%
#   arrange(desc(date))

