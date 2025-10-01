# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(tidyverse)
library(mgcv)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

#-------------------------------------------------------------------------------
# Load monthly data already matched with covariates
#-------------------------------------------------------------------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_wnnd <- 
  read.csv("....csv")

# only look at locations-years that have at least one case
df_wnnd_filtered <- df_wnnd %>%
  mutate(year = year(Date)) %>%
  group_by(NUTS_ID, year) %>%
  filter(any(wnnd_sum > 0)) %>%
  ungroup() 

df_wnnd_filtered <- df_wnnd_filtered %>%
  mutate(wnnd_pa = ifelse(wnnd_sum>0,1,0))

#-------------------------------------------------------------------------------
# Load GAM analysis and plotting function 
#-------------------------------------------------------------------------------

source(".../helper scripts/gam_analysis.R")
source(".../helper scripts/plot_gam_f.R")

#-------------------------------------------------------------------------------
# Apply GAM analysis and plotting function to temperature/R0 covariates
#-------------------------------------------------------------------------------

### 1. Temperature as predictor, incidence as outcome

temp_analysis <- gam_analysis(df_wnnd_filtered, 
                              "temp_3mo_avg")
temp_plots <- plot_gam_f(df_wnnd_filtered, temp_analysis$model, "temp_3mo_avg")

temp_plots$model_plot <- temp_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) 

temp_plots$data_plot <- temp_plots$data_plot + 
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  annotate("text", x = 10, 
           y = max(df_wnnd_filtered$incidence_per_100K, na.rm = TRUE),
           label = "tau[B] == 0.413", hjust = 0, size = 4, parse = T)
temp_plots[[3]]

### 2. Temperature as predictor, presence/absence as outcome

temp_binary_analysis <- gam_analysis(df_wnnd_filtered, 
                                     "temp_3mo_avg", binary = T)
temp_binary_plots <- plot_gam_f(df_wnnd_filtered,
                                temp_binary_analysis$model, "temp_3mo_avg",
                                binary = T)

temp_binary_plots$model_plot <- temp_binary_plots$model_plot +  
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) 
temp_binary_plots$model_plot

### 3. R0 from daily NUTS3 temperature as predictor, incidence as outcome

R0_from_daily_analysis <- gam_analysis(df_wnnd_filtered, "r0_from_daily_3mo_avg")
R0_from_daily_plots <- plot_gam_f(df_wnnd_filtered,
                                  R0_from_daily_analysis$model,
                                  "r0_from_daily_3mo_avg", r0_predictor = T)

R0_from_daily_plots$model_plot <- R0_from_daily_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, d"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_daily_plots$data_plot <- R0_from_daily_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, d"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

### 4. R0 from monthly NUTS3 temperature as predictor, incidence as outcome

R0_from_monthly_analysis <- gam_analysis(df_wnnd_filtered, 
                                         "r0_from_monthly_3mo_avg")
R0_from_monthly_plots <- plot_gam_f(df_wnnd_filtered,
                                    R0_from_monthly_analysis$model,
                                  "r0_from_monthly_3mo_avg", r0_predictor = T)

R0_from_monthly_plots$model_plot <- R0_from_monthly_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, m"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_monthly_plots$data_plot <- R0_from_monthly_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, m"}))
  )  +  
  scale_x_continuous(limits = c(0,1))

### 5. R0 from "seasonal" NUTS3 temperature as predictor, incidence as outcome

R0_from_seasonal_analysis <- gam_analysis(df_wnnd_filtered,
                                          "r0_from_seasonal_3mo_avg")
R0_from_seasonal_plots <- plot_gam_f(df_wnnd_filtered,
                                     R0_from_seasonal_analysis$model,
                                    "r0_from_seasonal_3mo_avg", r0_predictor = T)

# With shorter label notation

R0_from_seasonal_plots$model_plot_short <- R0_from_seasonal_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_seasonal_plots$data_plot_short <- R0_from_seasonal_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) +  
  annotate("text", x = 0, 
           y = max(df_wnnd_filtered$incidence_per_100K, na.rm = TRUE),
           label = "tau[B] == 0.44", hjust = 0, size = 4, parse = T)

# With long label notation

R0_from_seasonal_plots$model_plot <- R0_from_seasonal_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_seasonal_plots$data_plot <- R0_from_seasonal_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

### 6. R0 from daily gridded temperature as predictor, incidence as outcome

R0_from_daily_gridded_analysis <- gam_analysis(df_wnnd_filtered,
                                               "r0_from_daily_native_3mo_avg")
R0_from_daily_gridded_plots <- plot_gam_f(df_wnnd_filtered,
                                          R0_from_daily_gridded_analysis$model,
                                          "r0_from_daily_native_3mo_avg",
                                          r0_predictor = T)

R0_from_daily_gridded_plots$model_plot <- R0_from_daily_gridded_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R["0, grid"]^{"rel, d"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_daily_gridded_plots$data_plot <- R0_from_daily_gridded_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R["0, grid"]^{"rel, d"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

### 7. R0 from monthly gridded temperature as predictor, incidence as outcome

R0_from_monthly_gridded_analysis <- gam_analysis(df_wnnd_filtered,
                                                 "r0_from_monthly_native_3mo_avg")
R0_from_monthly_gridded_plots <- plot_gam_f(df_wnnd_filtered,
                                            R0_from_monthly_gridded_analysis$model,
                                          "r0_from_monthly_native_3mo_avg",
                                          r0_predictor = T)

R0_from_monthly_gridded_plots$model_plot <- R0_from_monthly_gridded_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R["0, grid"]^{"rel, m"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_monthly_gridded_plots$data_plot <- R0_from_monthly_gridded_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R["0, grid"]^{"rel, m"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

### 8. R0 from seasonal gridded temperature as predictor, incidence as outcome

R0_from_seasonal_gridded_analysis <- gam_analysis(df_wnnd_filtered,
                                                  "r0_from_seasonal_native_3mo_avg")
R0_from_seasonal_gridded_plots <- plot_gam_f(df_wnnd_filtered,
                                             R0_from_seasonal_gridded_analysis$model,
                                            "r0_from_seasonal_native_3mo_avg",
                                            r0_predictor = T)

R0_from_seasonal_gridded_plots$model_plot <- R0_from_seasonal_gridded_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R["0, grid"]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_seasonal_gridded_plots$data_plot <- R0_from_seasonal_gridded_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R["0, grid"]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

#-------------------------------------------------------------------------------
# Save ggplots
#-------------------------------------------------------------------------------

plot_list1 = list(temp_plots$model_plot, 
                  temp_plots$data_plot, 
                  temp_binary_plots$model_plot)

plot_grid1 = cowplot::plot_grid(plotlist = plot_list1, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B', 'C'))

plot_grid1

ggsave(".../main_gam_monthly_plot.tiff", 
       plot = plot_grid1,
       width = 8, height = 6, 
       dpi = 600, units = "in", compression = "lzw")

plot_list2 = list(R0_from_daily_plots$model_plot,
                  R0_from_daily_plots$data_plot,
                  R0_from_monthly_plots$model_plot,
                  R0_from_monthly_plots$data_plot,
                  R0_from_seasonal_plots$model_plot,
                  R0_from_seasonal_plots$data_plot)

plot_grid2 = cowplot::plot_grid(plotlist = plot_list2, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B',
                                                                    'C', 'D',
                                                                    'E', 'F'))

plot_grid2

ggsave(".../r0_gam_monthly_plot.tiff", 
       plot = plot_grid2,
       width = 8, height = 9, 
       dpi = 600, units = "in", compression = "lzw")

plot_list3 = list(R0_from_daily_gridded_plots$model_plot,
                  R0_from_daily_gridded_plots$data_plot,
                  R0_from_monthly_gridded_plots$model_plot,
                  R0_from_monthly_gridded_plots$data_plot,
                  R0_from_seasonal_gridded_plots$model_plot,
                  R0_from_seasonal_gridded_plots$data_plot)

plot_grid3 = cowplot::plot_grid(plotlist = plot_list3, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B',
                                                                    'C', 'D',
                                                                    'E', 'F'))

plot_grid3

ggsave(".../r0_gam_monthly_plot_native.tiff", 
       plot = plot_grid3,
       width = 8, height = 9, 
       dpi = 600, units = "in", compression = "lzw")

plot_list4 = list(temp_plots$model_plot, 
                  temp_plots$data_plot,
                  R0_from_seasonal_plots$model_plot_short,
                  R0_from_seasonal_plots$data_plot_short)

plot_grid4 = cowplot::plot_grid(plotlist = plot_list4, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B',
                                                                    'C', 'D'))

plot_grid4

ggsave(".../main_gam_monthly_plot_2.tiff", 
       plot = plot_grid4,
       width = 8, height = 6, 
       dpi = 600, units = "in", compression = "lzw")

ggsave(".../binary_gam_monthly.tiff", 
       plot = temp_binary_plots$model_plot,
       width = 4, height = 3, 
       dpi = 600, units = "in", compression = "lzw")

