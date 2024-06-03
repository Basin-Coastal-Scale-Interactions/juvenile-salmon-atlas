### Predict monthly distributions in select inlets/regions in more detail

# ideally we want a finer scale grid than the one used for the whole BC coast,
# but for now we can just subset that grid

library(tidyverse)
library(ggplot2)
library(here)
library(sdmTMB)

source(here("R", "plot_map.R"))

dat_tbl_mout3 <- readRDS(here("data", "fits", "all_mout3_allcoast.rds"))
dat_tbl_dmesh <- readRDS(here("data", "fits", "all_dmesh_allcoast.rds"))

index_grid <- readRDS(here("data", "pred_grid_all_coast.rds")) %>%
  mutate(year = 2019L, 
         year_f = as.factor(year),
         survey_f = as.factor("hss"),
         target_depth = 0,
         day_night = "DAY",
         scale_depth = -0.5427145)


# Barkley Sound
barkley_grid <- filter(index_grid,
                       lon <= -124.5 & lon >= -126 & 
                         lat >= 48.6 & lat <= 49.5)


# northern VI inlets
northvi_grid <- filter(index_grid,
                       lon <= -126.5 & lon >= -129 & 
                         lat >= 49.5 & lat <= 50.7)


# northern VI inlets
qcss_grid <- filter(index_grid, 
                    lon <= -125.5 & lon >= -132 & 
                      lat >= 50 & lat <= 53)




# looking at min and max lat/lon to see equivalent UTM coordinates
barkley_grid %>%
  slice(c(which.min(lon), which.max(lon)))
barkley_grid %>%
  slice(c(which.min(lat), which.max(lat)))

barkley_grid_months <- barkley_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)    

northvi_grid_months <- northvi_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)    


qcss_grid_months <- qcss_grid %>%
  replicate_df(., "month", 1:12) %>%
  mutate(month = as.numeric(month), 
         month_f = month(month, label = TRUE),
         month_adj = ifelse(month > 2, month - 2, month + 10), # adjusting to start in March
         NULL)    



barkley_pred_months <- purrr::map(
  dat_tbl_dmesh$fit,
  ~ predict(.x, newdata = barkley_grid_months, 
            offset = rep.int(-6.287114, nrow(barkley_grid_months)), re_form_iid = NA) 
)

northvi_pred_months <- purrr::map(
  dat_tbl_dmesh$fit,
  ~ predict(.x, newdata = northvi_grid_months, 
            offset = rep.int(-6.287114, nrow(northvi_grid_months)), re_form_iid = NA) 
)

qcss_pred_months <- purrr::map(
  dat_tbl_dmesh$fit,
  ~ predict(.x, newdata = qcss_grid_months, 
            offset = rep.int(-6.287114, nrow(qcss_grid_months)), re_form_iid = NA) 
)

### Barkley
for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
  p <- plot_map(barkley_pred_months[[i]], exp(est), area = "Barkley",
                show_raw_data = TRUE, dat_tbl_dmesh$data[[i]]) + 
    annotation_scale(location = "tr", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month around Barkley Sound (dense mesh)")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = c(718406, 812406), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5388664, 5477664), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
    
  png(here::here("figs", paste0("Barkley_atlas_dmesh_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}


northvi_grid %>%
  slice(c(which.min(lon), which.max(lon)))
northvi_grid %>%
  slice(c(which.min(lat), which.max(lat))) 

### North VI
for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
 p <- plot_map(northvi_pred_months[[i]], exp(est), area = "BC",
                show_raw_data = TRUE, dat_tbl_dmesh$data[[i]]) + 
    annotation_scale(location = "tr", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month around northern WCVI (dense mesh)")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = c(500406, 680406), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5483664, 5618664), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
  png(here::here("figs", paste0("northvi_atlas_dmesh_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}

# Queen Charlotte Strait & Sound

qcss_grid %>%
  slice(c(which.min(lon), which.max(lon)))
qcss_grid %>%
  slice(c(which.min(lat), which.max(lat))) 

for (i in 1:nrow(dat_tbl_dmesh)) {
  species <- dat_tbl_dmesh$species[i]
  print(species)
  
  p <- plot_map(qcss_pred_months[[i]], exp(est), area = "BC",
                show_raw_data = TRUE, dat_tbl_dmesh$data[[i]]) + 
    annotation_scale(location = "bl", height = unit(0.1, "cm")) +
    ggtitle(paste0("Predicted distribution of juvenile ", species," by month around QC Strait and Sound (fine mesh)")) +
    # added new plot ranges based on UTM
    scale_x_continuous(name = NULL, limits = c(292406, 724406), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(5538664, 5872664), expand = c(0, 0)) +
    theme(legend.key.height = unit(0.06, 'npc'))
  
  png(here::here("figs", paste0("qcss_atlas_dmesh_", species, "_juv_months.png")),
      height = 7, width = 10, units = "in", res = 300)
  print(p) # need to print() plot inside for loop!
  dev.off() 
  print("Done!")
}

