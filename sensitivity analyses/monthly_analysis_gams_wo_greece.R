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

# calculate some stats:
# highest incidence value after excluding Greece
highest_incidence_wo_Greece <- 
  max(df_wnnd_filtered$incidence_per_100K[!(grepl("^EL", df_wnnd_filtered$NUTS_ID))])
highest_incidence_wo_Greece

# % positive incidence observations above that threshold
sum(df_wnnd_filtered$incidence_per_100K > highest_incidence_wo_Greece) /
  sum(df_wnnd_filtered$incidence_per_100K > 0)

# % of observations above 24.4°C that are from Greece
sum(grepl("^EL", 
          df_wnnd_filtered$NUTS_ID[df_wnnd_filtered$temp_3mo_avg_pop_weight > 
                                     24.4])) / 
  sum(df_wnnd_filtered$temp_3mo_avg_pop_weight > 24.4)

# filter out greece
df_wnnd_filtered <- df_wnnd_filtered %>%
  filter(!grepl("^EL", NUTS_ID))

df_wnnd_filtered$country <- factor(substr(df_wnnd_filtered$NUTS_ID, 1, 2))

#-------------------------------------------------------------------------------
# Load GAM analysis and plotting function 
#-------------------------------------------------------------------------------

source("code/helper scripts/gam_analysis.R")
source("code/helper scripts/plot_gam_f.R")

#-------------------------------------------------------------------------------
# Apply GAM analysis and plotting function to temperature/R0 covariates
#-------------------------------------------------------------------------------

### 1. Temperature as predictor, incidence as outcome

temp_analysis <- gam_analysis(df = df_wnnd_filtered, 
                              predictor = "temp_3mo_avg_pop_weight")

# generate plots
temp_plots <- plot_gam_f(df = df_wnnd_filtered, 
                         model = temp_analysis$model, 
                         predictor = "temp_3mo_avg_pop_weight")

temp_plots$shared_model_plot <- temp_plots$shared_model_plot +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2))

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

temp_plots$shared_data_plot <- temp_plots$shared_data_plot + 
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

# generate plots
temp_plots_adjusted <- plot_gam_f(df = df_wnnd_filtered, 
                                  model = temp_analysis_adjusted$model, 
                                  predictor="temp_3mo_avg_pop_weight", 
                                  adjustments = c("urban_fabric", 
                                                  "rain_3mo_cum_pop_weight", 
                                                  "arable_land",
                                                  "percent_65plus",
                                                  "gdp_per_capita"
                                  ))

temp_plots_adjusted$shared_model_plot <- temp_plots_adjusted$shared_model_plot +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) 

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

temp_plots_adjusted$shared_data_plot <- temp_plots_adjusted$shared_data_plot + 
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  ) 

temp_plots_adjusted

### 3. R0 from "seasonal" NUTS3 temperature as predictor, incidence as outcome

R0_from_seasonal_analysis <- gam_analysis(df = df_wnnd_filtered,
                                          predictor = "r0_from_seasonal_3mo_avg_pop_weight",
                                          k_shared = 10, k_country = 3)

# generate plots
R0_from_seasonal_plots <- plot_gam_f(df = df_wnnd_filtered,
                                     model = R0_from_seasonal_analysis$model,
                                     predictor = "r0_from_seasonal_3mo_avg_pop_weight", 
                                     r0_predictor = T)

R0_from_seasonal_plots$shared_model_plot <- R0_from_seasonal_plots$shared_model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)
  )

R0_from_seasonal_plots$model_plot <- R0_from_seasonal_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)
  )

R0_from_seasonal_plots$data_plot <- R0_from_seasonal_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)) 

R0_from_seasonal_plots$shared_data_plot <- R0_from_seasonal_plots$shared_data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)) 

R0_from_seasonal_plots

### 4. R0 from "seasonal" NUTS3 temperature as predictor adjusted for additional variables, incidence as outcome

R0_from_seasonal_analysis_adjusted <- gam_analysis(df = df_wnnd_filtered,
                                                   predictor = "r0_from_seasonal_3mo_avg_pop_weight",
                                                   adjustments = c("urban_fabric", 
                                                                   "rain_3mo_cum_pop_weight", 
                                                                   "percent_65plus", 
                                                                   "arable_land", 
                                                                   "gdp_per_capita"),
                                                   k_shared = 10, k_country = 3)

# generate plots
R0_from_seasonal_plots_adjusted <- plot_gam_f(df = df_wnnd_filtered,
                                              model = R0_from_seasonal_analysis_adjusted$model,
                                              predictor = "r0_from_seasonal_3mo_avg_pop_weight", 
                                              adjustments = c("urban_fabric", 
                                                              "rain_3mo_cum_pop_weight", 
                                                              "percent_65plus", 
                                                              "arable_land", 
                                                              "gdp_per_capita"),
                                              r0_predictor = T)

R0_from_seasonal_plots_adjusted$shared_model_plot <- R0_from_seasonal_plots_adjusted$shared_model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)
  )

R0_from_seasonal_plots_adjusted$model_plot <- R0_from_seasonal_plots_adjusted$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)
  )

R0_from_seasonal_plots_adjusted$data_plot <- R0_from_seasonal_plots_adjusted$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)) 

R0_from_seasonal_plots_adjusted$shared_data_plot <- R0_from_seasonal_plots_adjusted$shared_data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel"}))
  )  +  
  coord_cartesian(xlim = c(0,1)) 

R0_from_seasonal_plots_adjusted

#-------------------------------------------------------------------------------
# Save plots
#-------------------------------------------------------------------------------

plot_list = list(temp_plots$shared_model_plot + ggtitle("Unadjusted") +
                   coord_cartesian(ylim = c(0, 1.1)) +
                   theme(plot.title = element_text(size = 8)), 
                 temp_plots_adjusted$shared_model_plot + ggtitle("Adjusted") +
                   coord_cartesian(ylim = c(0, 1.1)) +
                   theme(plot.title = element_text(size = 8)),
                 temp_plots$shared_data_plot + ggtitle("Data") +
                   theme(plot.title = element_text(size = 8)),
                 R0_from_seasonal_plots$shared_model_plot + ggtitle("Unadjusted") +
                   coord_cartesian(ylim = c(0, 0.3)) +
                   theme(plot.title = element_text(size = 8)),
                 R0_from_seasonal_plots_adjusted$shared_model_plot + ggtitle("Adjusted") +
                   coord_cartesian(ylim = c(0, 0.3)) +
                   theme(plot.title = element_text(size = 8)),
                 R0_from_seasonal_plots$shared_data_plot + ggtitle("Data") +
                   theme(plot.title = element_text(size = 8)))

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol = 3, label_size = 10,
                               align = "h", axis = "b", labels = c('A', 'B',
                                                                   'C', 'D',
                                                                   'E', 'F'))

plot_grid

ggsave("figures/sensitivity/Figure_S6.tiff",
       plot = plot_grid,
       bg = "white",
       width = 6.3, height = 3.83,
       dpi = 600, units = "in", compression = "lzw")
ggsave("figures/sensitivity/Figure_S6.svg",
       plot = plot_grid,
       bg = "white",
       width = 6.3, height = 3.83,
       units = "in")

