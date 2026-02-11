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

df_wnnd <- 
  read.csv("data/WNND_data_with_covariates.csv")

# only look at locations-years that have at least one case
df_wnnd_filtered <- df_wnnd %>%
  mutate(year = year(Date)) %>%
  group_by(NUTS_ID, year) %>%
  filter(any(wnnd_sum > 0)) %>%
  ungroup() 

#-------------------------------------------------------------------------------
# Load GAM analysis and plotting function 
#-------------------------------------------------------------------------------

source("code/helper scripts/gam_analysis.R")
source("code/helper scripts/plot_gam_f.R")

#-------------------------------------------------------------------------------
# Apply GAM analysis and plotting function to temperature/R0 covariates
#-------------------------------------------------------------------------------

### 1. Temperature as predictor, incidence as outcome
gam_analysis(df = df_wnnd_filtered, predictor = "temp")
gam_analysis(df = df_wnnd_filtered, predictor = "temp_2mo_avg")
gam_analysis(df = df_wnnd_filtered, predictor = "temp_4mo_avg")
gam_analysis(df = df_wnnd_filtered, predictor = "temp_3mo_avg")
gam_analysis(df = df_wnnd_filtered, predictor = "temp_pop_weight")
gam_analysis(df = df_wnnd_filtered, predictor = "temp_2mo_avg_pop_weight")
gam_analysis(df = df_wnnd_filtered, predictor = "temp_4mo_avg_pop_weight")

temp_analysis <- gam_analysis(df = df_wnnd_filtered, 
                              predictor = "temp_3mo_avg_pop_weight")

temp_plots <- plot_gam_f(df = df_wnnd_filtered, model = temp_analysis$model, 
                         predictor = "temp_3mo_avg_pop_weight")

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
  ) 

temp_plots

### 2. Temperature as predictor adjusted for additional variables, incidence as outcome

temp_analysis_adjusted <- gam_analysis(df = df_wnnd_filtered, 
                                       predictor = "temp_3mo_avg_pop_weight", 
                                       adjustments = c("urban_fabric", 
                                         "rain_3mo_cum_pop_weight", 
                                         "percent_65plus", 
                                         "arable_land",
                                         "gdp_per_capita"
                                       ))

temp_plots_adjusted <- plot_gam_f(df = df_wnnd_filtered, 
                                  model = temp_analysis_adjusted$model, 
                                  predictor="temp_3mo_avg_pop_weight", 
                                  adjustments = c("urban_fabric", 
                                    "rain_3mo_cum_pop_weight", 
                                    "arable_land",
                                    "percent_65plus",
                                    "gdp_per_capita"
                                  ))

temp_plots_adjusted$model_plot <- temp_plots_adjusted$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) 

temp_plots_adjusted$data_plot <- temp_plots_adjusted$data_plot + 
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  ) 
temp_plots_adjusted

### 3. R0 from "seasonal" NUTS3 (window/nuts3) temperature as predictor, 
# incidence as outcome

gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_monthly")
gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_seasonal_2mo_avg")
gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_seasonal_4mo_avg")
gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_seasonal_3mo_avg")
gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_monthly_pop_weight")
gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_seasonal_2mo_avg_pop_weight")
gam_analysis(df = df_wnnd_filtered, predictor = "r0_from_seasonal_4mo_avg_pop_weight")

R0_from_seasonal_analysis <- gam_analysis(df = df_wnnd_filtered,
                                          predictor = "r0_from_seasonal_3mo_avg_pop_weight")

R0_from_seasonal_plots <- plot_gam_f(df = df_wnnd_filtered,
                                     model = R0_from_seasonal_analysis$model,
                                     predictor = "r0_from_seasonal_3mo_avg_pop_weight", 
                                     r0_predictor = T)

R0_from_seasonal_plots$model_plot <- R0_from_seasonal_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0, 1.2))

R0_from_seasonal_plots$data_plot <- R0_from_seasonal_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)) 

R0_from_seasonal_plots

### 4. R0 from "seasonal" NUTS3 (window/nuts3) temperature as predictor 
# adjusted for additional variables, incidence as outcome

R0_from_seasonal_analysis_adjusted <- gam_analysis(df = df_wnnd_filtered,
                                                   predictor = "r0_from_seasonal_3mo_avg_pop_weight",
                                                   adjustments = c("urban_fabric", 
                                                                   "rain_3mo_cum_pop_weight", 
                                                                   "percent_65plus", 
                                                                   "arable_land", 
                                                                   "gdp_per_capita"))

R0_from_seasonal_plots_adjusted <- plot_gam_f(df = df_wnnd_filtered,
                                              model = R0_from_seasonal_analysis_adjusted$model,
                                              predictor = "r0_from_seasonal_3mo_avg_pop_weight", 
                                              adjustments = c("urban_fabric", 
                                                              "rain_3mo_cum_pop_weight", 
                                                              "percent_65plus", 
                                                              "arable_land", 
                                                              "gdp_per_capita"),
                                              r0_predictor = T)

R0_from_seasonal_plots_adjusted$model_plot <- R0_from_seasonal_plots_adjusted$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0, 1.2))

R0_from_seasonal_plots_adjusted$data_plot <- R0_from_seasonal_plots_adjusted$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)) 

R0_from_seasonal_plots_adjusted

#-------------------------------------------------------------------------------
# Save plots
#-------------------------------------------------------------------------------

plot_list = list(temp_plots$model_plot + ggtitle("Unadjusted") +
                    theme(plot.title = element_text(size = 12)), 
                  temp_plots_adjusted$model_plot + ggtitle("Adjusted") +
                    theme(plot.title = element_text(size = 12)),
                  temp_plots$data_plot + ggtitle("Raw data") +
                    theme(plot.title = element_text(size = 12)),
                  R0_from_seasonal_plots$model_plot + ggtitle("Unadjusted") +
                    theme(plot.title = element_text(size = 12)),
                  R0_from_seasonal_plots_adjusted$model_plot + ggtitle("Adjusted") +
                    theme(plot.title = element_text(size = 12)),
                  R0_from_seasonal_plots$data_plot + ggtitle("Raw data") +
                    theme(plot.title = element_text(size = 12)))

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=3, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B',
                                                                    'C', 'D',
                                                                    'E', 'F'))

plot_grid

ggsave("figures/monthly/Figure_4.tiff",
       plot = plot_grid4,
       width = 11.5, height = 7,
       dpi = 300, units = "in", compression = "lzw")
ggsave("figures/monthly/Figure_4.svg",
       plot = plot_grid4,
       width = 11.5, height = 7,
       units = "in")

