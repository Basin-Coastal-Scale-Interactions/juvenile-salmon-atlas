# Comparison of sdmTMB and RTMB stock composition models

library(tidyverse)
library(sdmTMB)
library(RTMB)
library(ggsidekick)
library(here)

source("R/plot_map.R")

m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_sdmTMB.rds"))
pred <- readRDS("data/fits/gsi-prediction-normalized.rds")
prop_grid_all <- readRDS(here::here("data", "pred_grid_stock_prop_RTMB.rds"))

# downsizing data frames
pred_sdmTMB <- select(pred, utm_x_1000, utm_y_1000, month_adj, salmon_region = region, pred_st = pred_i, prob_st = prob_i)
pred_RTMB <- select(prop_grid_all, utm_x_1000, utm_y_1000, month_adj, salmon_region = region, pred_rt = pred,  prob_rt = prob)

glimpse(pred)
glimpse(prop_grid_all)

unique(pred_sdmTMB$month_adj)
unique(pred_RTMB$month_adj)

# removing months that are not in sdmTMB pred
pred_RTMB <- filter(pred_RTMB, month_adj %in% 1:9)

dim(pred_RTMB)
dim(pred_sdmTMB)

pj <- full_join(pred_sdmTMB, pred_RTMB, by = c("utm_x_1000", "utm_y_1000", "month_adj", "salmon_region")) %>%
  #mutate(X = utm_x_1000 * 1000, Y = utm_y_1000 * 1000) %>%
  mutate(diff_smr = prob_st - prob_rt)

ggplot(data = pj) +
  geom_point(aes(pred_st, pred_rt, color = salmon_region), size = 0.5) +
  facet_wrap(~month_adj, scales = "free") +
  ggsidekick::theme_sleek()

ggsave(here("figs", "stock-prop-pred-comparison.png"), width = 8, height = 7, units = "in")

# Map of difference in predicted proportion
pdiff <- pj %>%
  mutate(month_f = month(month_adj + 2 , label = TRUE, abbr = TRUE)) %>%
  ggplot(data = .) + # (dat, aes(X, Y, color = {{ column }})) +
  geom_raster( aes(x=utm_x_1000, y=utm_y_1000, fill = diff_smr)) +
  scale_fill_gradient2(name = "Difference in proportion of\nGSI assignment between\nsdmTMB and RTMB (RW by\ngroup) models",#"Predicted\nnumber of\njuveniles\nin GSI samples",
                       labels = scales::comma) +
  ggsidekick::theme_sleek() +
  geom_sf(data = all_coast_km, color = "gray80", fill = "gray90") +
  scale_x_continuous(name = NULL, limits = range(pj$utm_x_1000), expand = c(0, 0)) +
  scale_y_continuous(name = NULL, limits = range(pj$utm_y_1000), expand = c(0, 0)) +
  facet_grid(salmon_region ~ month_f)  +
  ggtitle("Juvenile chinoook difference in GSI proportion")

pdiff

ggsave(pdiff, filename = here::here("figs", "chinook_by_gsi_diff_sdmTMB_RTMB.png"), 
       width = 14, height = 12)

###

sims <- simulate(m_svc, n = 10, type = "mle-mvn")

glimpse(sims)
head(sims)
dim(sims)
