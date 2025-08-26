# Wrangle table of yearly tows counts by ship and net

library(tidyverse)
library(readxl)
library(janitor)
library(here)
library(knitr)
library(kableExtra)
options(knitr.kable.NA = '')


# All tows
tow_info <- read_csv(here("data-raw", "BCSI_TowInfo_20250402.csv")) %>%
  clean_names() 

tow_info$date %>% range()

# No day_night
filter(tow_info, unique_event == "IPES2024-026-126") %>% glimpse


# Fixing vessel names
tow_info$vessel <- forcats::fct_recode(
  tow_info$vessel_name,
  "CCGS Neocaligus" = "CCGS NEOCALIGUS",
  "CFV Sea Crest" = "CFV SEACREST",
  "CFV Sea Crest" = "CFV SEA CREST",
  "CCGS Sir John Franklin" = "CCGS SIR JOHN FRANKLIN",
  "CFV Nordic Pearl" = "Nordic Pearl",
  "FV Anita J" = "FV ANITA J",
  "CCGS W.E. Ricker" = "CCGS W.E. RICKER",
  "CFV Ocean Selector" =   "CFV OCEAN SELECTOR",
  "CFV Frosti" = "CFV FROSTI",
  "CFV Viking Storm" = "CFV VIKING STORM",
  "CFV Viking Storm & CCGS W.E. Ricker" = "CFV VIKING STORM AND CCGS W.E. RICKER",
  "CFV Viking Storm & CCGS W.E. Ricker" = "CFV VIKING STORM & CCGS W.E. RICKER",
  "FV Columbia" =   "FV COLUMBIA",
  "CFV Sea Crest" =   "SEA CREST",
  "CFV Nordic Pearl" =   "NORDIC PEARL",
  "CCGS Sir John Franklin" = "Sir John Franklin"
)

unique(tow_info$vessel) |> as.character()


saveRDS(tow_info, file = here("data", "tow_info.rds"))
      
dat <- readRDS(here::here("data", "chinook_dat_allcoast.rds"))
gsidat <- readRDS(here::here("data", "chinook_gsi_counts_fitted.rds"))
dat$year %>% range()
gsidat$year %>% range()

dat$date %>% range()
gsidat$date %>% range()

unique(dat$unique_event) %>% length
unique(gsidat$unique_event) %>% length


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
gsitow_data <- gsidat |>
  select(-net_desc) |> 
  left_join(t_info, by = "unique_event")

# Checking tows with NAs in net_desc
filter(tow_data, is.na(net_desc))

filter(gsitow_data, grepl("Neocaligus", vessel)) %>%
  pull(year) |> table()

# Fixing net values here
# JK 7 November, 1:06 pm
# "The Ocean Selector is a groundfish boat, so we would not have used a bottom 
# trawl net.  I would expect the CanTrawl 250 to have been loaded."
# 
# Amy: LFS 1142 is the big offshore net they were testing.  they did a couple tows 
# with the codend open, just to test the gear etc, (would have been marked unusable), 
# then did some full tows with catch processed. 
# 2020 would have been the testing, 2022 would have been for the actual IYS survey 
# in which the whole survey used that net
# 
#   the 4x4 is even smaller than LFS 7742. 
# the net desc for the zub mini was updated to the CanTrawl 75.  
# 
#  Ocean Selector 2002:"Ocean Selector's net".
# "On this survey, the F.V. Ocean Selector towed a commercial mid-water trawl 
# which was modified for juvenile salmon fishing by the attachment of an 
# intermediate section of 3 in. (7.6 cm) polypropylene and a codend of 1.5 in. 
# (3.8 cm) knotted nylon lined with 0.25 in. mesh (64 mm) that were cut from 
# the Highseas Salmon program’s 2 mid-water trawl. The F.V. Ocean Selector 
# towed this modified mid-water trawl at 4-5 knots with Thyboron 92 type 3 
# trawl doors, and a 150 m pay-out of 1” Briden trawl warp per side. The F.V. 
# Ocean Selector achieved a mouth opening of approximately 28 m wide by 16 m 
# deep with this gear configuration and towing speed, which was comparable to 
# that achieved by the CCGS W.E. Ricker."


td <- tow_data %>%
  mutate(net_desc = case_when(
    year >= 2021 & is.na(net_desc) ~ "LFS 7742",
    net_desc == "Ocean Selectors net" ~ "Modified mid-water trawl",
    net_desc == "Rusty's backup net" ~ "CanTrawl 250",
    .default = net_desc
  )) %>%
  group_by(year, vessel, net_desc) %>%
  summarise(n = n()) %>%
  arrange(year)


print(td, n = 48)

gtd <- gsitow_data %>%
  mutate(net_desc = case_when(
    year >= 2021 & is.na(net_desc) ~ "LFS 7742",
    net_desc == "Ocean Selectors net" ~ "Modified mid-water trawl",
    net_desc == "Rusty's backup net" ~ "CanTrawl 250",
    net_desc == "Zubkowski mini 250" ~ "CanTrawl 75",
    year == 2017 & vessel == "CFV Sea Crest" ~  "CanTrawl 250",
    .default = net_desc
  )) %>%
  group_by(year, vessel, net_desc) %>%
  summarise(n = n()) %>%
  mutate(n = n/11) |>
  arrange(year)

atd <-full_join(td, gtd, by = join_by(year, vessel, net_desc)) |>
  arrange(year,vessel) |> 
  mutate(nx2 = ifelse(is.na(n.x), "-", as.character(n.x)),
         ny2 = ifelse(is.na(n.y), "-", as.character(n.y)),
         n = paste0(nx2, " (", ny2, ")")) 

print(select(atd, year, vessel, net_desc, n), n=50)


tdw <- atd |>
  select(year, vessel, net_desc, n) |>
  pivot_wider(names_from = net_desc, values_from = n) %>% print(n = 50) |>
  rename(Year = year,
         Vessel = vessel)


tow_kb <- kable(tdw, format = "latex", align = "llccccccc",
                  caption = "\\label{tab:tows}Yearly vessel and net type for all tows used in this 
                study, with numbers on the left being the count for tows used in the abundance model, 
                while numbers in parentheses are the count for tows used in the stock composition model. 
                Note that CanTrawl 400/580 and LFS 1142 are much larger nets than 
                CanTrawl 250 and LFS 7742, with these latter two being comparable in size, 
                while the 4 x 4 Trawl and CanTrawl 75 are much smaller and not comparable 
                in terms of fishing effort, thus samples from these are only used for 
                the stock composition model") %>%
  kable_styling(font_size = 8) %>%
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

### Size thresholds table

size_th <- read_csv(here("data-raw", "size_thresholds_chinook.csv")) %>%
  mutate(size_threshold = paste("\\textless", size_threshold))

size_kb <- kable(list(size_th[1:6,], size_th[7:12,]), format = "latex", align = "lr",
                 label = "size_th", escape = FALSE,  booktabs = TRUE,
                 col.names = c("Month", "Length (mm FL)"),
                caption = "Monthly size thresholds (fork length; mm) used to 
                classify juvenile chinook. These thresholds include both 
                subyearling and yearling fish.") %>%
  #kable_styling(font_size = 6) %>%
  kable_styling(latex_options = "hold_position")

writeLines(size_kb, here("tables","size_thresholds.tex"))


