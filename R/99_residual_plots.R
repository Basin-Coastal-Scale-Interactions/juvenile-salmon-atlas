library(sdmTMB)
library(here)
library(tidyverse)

# loading sdmTMB model fits
m_chinook <- readRDS(here::here("data", "fits", "fits_list_mdmesh.rds"))$chinook
m_svc <- readRDS(here::here("data", "fits", "chinook_gsi_prop_svc_sdmTMB.rds"))

# checking that object is correct version and data
m_svc$version
table(m_svc$data$region)

s_ch <- simulate(m_chinook, nsim = 200, type = "mle-mvn")
dr_ch <- dharma_residuals(s_ch, m_chinook, return_DHARMa = TRUE,
                          test_uniformity = FALSE,
                          test_outliers = FALSE)

DHARMa::plotQQunif(simulationOutput = dr_ch, 
           testDispersion = FALSE,
           testUniformity = FALSE,
           testOutliers = FALSE,
           main = "DHARMA QQ residuals for abundance model")

plot(dr_ch, title = "DHARMA residuals for abundance model")

s_svc <- simulate(m_svc, nsim = 200, type = "mle-mvn")
dr_svc <- dharma_residuals(s_svc, m_svc, return_DHARMa = TRUE)

plot(dr_svc, title = "DHARMA residuals for stock-specific model")

DHARMa::plotQQunif(simulationOutput = dr_svc, 
                   testDispersion = FALSE,
                   testUniformity = FALSE,
                   testOutliers = FALSE,
                   main = "DHARMA QQ residuals for stock composition model")

png(here::here("figs", "dharma-species-wide.png"), width = 10, height = 5,
    res= 300, units = "in")
plot(dr_ch, title = "DHARMA residuals for species-wide model")
dev.off()

png(here::here("figs", "dharma-stock-specific.png"), width = 10, height = 5,
    res= 300, units = "in")
plot(dr_svc, title = "DHARMA residuals for stock-specific model")
dev.off()

png(here::here("figs", "dharma-stock-specific2.png"), width = 10, height = 5,
    res= 300, units = "in")
par(mfrow=c(1,2))

DHARMa::plotQQunif(simulationOutput = dr_ch, 
                   testDispersion = FALSE,
                   testUniformity = FALSE,
                   testOutliers = FALSE,
                   main = "DHARMa QQ residuals for abundance model")

DHARMa::plotQQunif(simulationOutput = dr_svc, 
                   testDispersion = FALSE,
                   testUniformity = FALSE,
                   testOutliers = FALSE,
                   main = "DHARMa QQ residuals for stock composition model")

par(mfrow = c(1,1))
dev.off()
