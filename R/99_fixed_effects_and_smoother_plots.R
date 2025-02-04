library(sdmTMB)
library(here)
library(tidyverse)

# loading sdmTMB model fits
m_chinook <- readRDS(here::here("data", "fits", "fits_list_mdmesh.rds"))$chinook
m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))



fix_pars <- tidy(m_chinook, effects = "fixed", conf.int = T) %>% 
  filter(term %in% c("day_nightDAY", "day_nightNIGHT", "survey_fipes", "scale_depth")) %>% 
  mutate(
    term = fct_recode(as.factor(term), "Headrope Depth" = "scale_depth", 
                      "IPES Survey" = "survey_fipes", 
                      "Diurnal Sampling" = "day_nightDAY",
                      "Nocturnal Sampling" = "day_nightNIGHT")
  )

fix_svc_pars <- tidy(m_svc, effects = "fixed", conf.int = T) %>%
  mutate(term = as.factor(sub("region", "", term)))

fix_pars

# need to address rows with empty values in `day_night`
table(m_chinook$data$day_night)

fp1 <- ggplot(fix_pars, aes(term, estimate, ymin = conf.low, ymax = conf.high)) +
  ggsidekick::theme_sleek() +
  geom_pointrange() +
  guides(fill = "none") +
  #theme(axis.title = element_blank()) +
  labs(y = "Parameter Estimate", x = "") +
  geom_hline(yintercept = 0, lty = 2, colour = "red") + 
  coord_flip()


fp2 <- ggplot(fix_svc_pars, aes(term, estimate, ymin = conf.low, ymax = conf.high)) +
  ggsidekick::theme_sleek() +
  geom_pointrange() +
  guides(fill = "none") +
  labs(y = "Parameter Estimate", x = "") +
  geom_hline(yintercept = 0, lty = 2, colour = "red") + coord_flip()


cowplot::plot_grid(fp1, fp2, labels = "AUTO")

ggsave(here("figs", "fixed_effects.png"), width = 8, height = 3, units = "in")


###############


psmooth_species <- plot_smoother(m_chinook)

psmooth_species2 <- psmooth_species+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

psmooth_gsi <- plot_smoother(m_svc, regions = TRUE)

psmooth_gsi2 <- psmooth_gsi + facet_wrap(~region, scales = "free_y", nrow = 5) +
  theme(plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm"))

ggsave(here("figs", "smoother_species.png"), psmooth_species2, width = 6, height = 3, units = "in")

ggsave(here("figs", "smoother_gsi.png"), psmooth_gsi2, width = 7, height = 9, units = "in")
