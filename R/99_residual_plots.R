library(sdmTMB)
library(here)

# loading sdmTMB model fits
m_chinook <- readRDS(here::here("data", "fits", "fits_list_mdmesh.rds"))$chinook
m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))

s_ch <- simulate(m_chinook, nsim = 200, type = "mle-mvn")
dr_ch <- dharma_residuals(s_ch, m_chinook, return_DHARMa = TRUE)

plot(dr_ch, title = "DHARMA residuals for species-wide model")

s_svc <- simulate(m_svc, nsim = 200, type = "mle-mvn")
dr_svc <- dharma_residuals(s_svc, m_svc, return_DHARMa = TRUE)

plot(dr_svc, title = "DHARMA residuals for stock-specific model")



png(here::here("figs", "dharma-species-wide.png"), width = 10, height = 5,
    res= 300, units = "in")
plot(dr_ch, title = "DHARMA residuals for species-wide model")
dev.off()

png(here::here("figs", "dharma-stock-specific.png"), width = 10, height = 5,
    res= 300, units = "in")
plot(dr_svc, title = "DHARMA residuals for stock-specific model")
dev.off()

