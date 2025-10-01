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

#-------------------------------------------------------------------------------
# Load GAM analysis and plotting function 
#-------------------------------------------------------------------------------

source(".../helper scripts/gam_analysis_conf_urban.R")
source(".../helper scripts/plot_gam_conf_urban_f.R")
source(".../helper scripts/gam_analysis_conf_all.R")
source(".../helper scripts/plot_gam_conf_all_f.R")

#-------------------------------------------------------------------------------
# GAM analysis with temperature and urban fabric as predictors
#-------------------------------------------------------------------------------

temp_urban_analysis <- gam_analysis_conf_urban(df_wnnd_filtered,
                                               "temp_3mo_avg")
temp_urban_plots <- plot_gam_conf_urban_f(df_wnnd_filtered, 
                                          temp_urban_analysis$model,
                                          "temp_3mo_avg")

temp_urban_plots$model_plot <- temp_urban_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) 

temp_urban_plots$data_plot <- temp_urban_plots$data_plot + 
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  ) 
temp_urban_plots[[3]]

#-------------------------------------------------------------------------------
# GAM analysis with R0 from seasonal and urban fabric as predictors
#-------------------------------------------------------------------------------

R0_from_seasonal_urban_analysis <- gam_analysis_conf_urban(df_wnnd_filtered,
                                               "r0_from_seasonal_3mo_avg")
R0_from_seasonal_urban_plots <- plot_gam_conf_urban_f(df_wnnd_filtered,
                                                      R0_from_seasonal_urban_analysis$model,
                                                      "r0_from_seasonal_3mo_avg",
                                                      r0_predictor = T)

R0_from_seasonal_urban_plots$model_plot <- R0_from_seasonal_urban_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_seasonal_urban_plots$data_plot <- R0_from_seasonal_urban_plots$data_plot + 
  scale_x_continuous(limits = c(0,1)) +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  ) 
R0_from_seasonal_urban_plots[[3]]

#-------------------------------------------------------------------------------
# GAM analysis with R0 from seasonal and "all" confounders as predictors
#-------------------------------------------------------------------------------

temp_all_analysis <- gam_analysis_conf_all(df_wnnd_filtered,
                                           "temp_3mo_avg")
temp_all_plots <- plot_gam_conf_all_f(df_wnnd_filtered, 
                                          temp_all_analysis$model,
                                          "temp_3mo_avg")

temp_all_plots$model_plot <- temp_all_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) 

temp_all_plots$data_plot <- temp_all_plots$data_plot + 
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  ) 
temp_all_plots[[3]]

#-------------------------------------------------------------------------------
# GAM analysis with R0 from seasonal and "all" confounders as predictors
#-------------------------------------------------------------------------------

R0_from_seasonal_all_analysis <- gam_analysis_conf_all(df_wnnd_filtered,
                                                           "r0_from_seasonal_3mo_avg")
R0_from_seasonal_all_plots <- plot_gam_conf_all_f(df_wnnd_filtered,
                                                      R0_from_seasonal_all_analysis$model,
                                                      "r0_from_seasonal_3mo_avg",
                                                      r0_predictor = T)

R0_from_seasonal_all_plots$model_plot <- R0_from_seasonal_all_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

R0_from_seasonal_all_plots$data_plot <- R0_from_seasonal_all_plots$data_plot + 
  scale_x_continuous(limits = c(0,1)) +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  ) 
R0_from_seasonal_all_plots[[3]]

#-------------------------------------------------------------------------------
# Save ggplots
#-------------------------------------------------------------------------------

plot_list1 = list(temp_urban_plots$model_plot, 
                  temp_urban_plots$data_plot, 
                  R0_from_seasonal_urban_plots$model_plot,
                  R0_from_seasonal_urban_plots$data_plot)

plot_grid1 = cowplot::plot_grid(plotlist = plot_list1, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B', 'C'))

title1 <- cowplot::ggdraw() +
  cowplot::draw_label("GAMs adj. for urban fabric %", fontface = 'bold', x = 0.5, hjust = 0.5, size = 12) 

plot_grid1 <- cowplot::plot_grid(title1, plot_grid1, ncol = 1, rel_heights = c(0.05, 1))

plot_grid1

ggsave(".../confounder_gam_urban.tiff", 
       plot = plot_grid1,
       bg = "white",
       width = 8, height = 6, 
       dpi = 300, units = "in", compression = "lzw")

plot_list2 = list(temp_all_plots$model_plot, 
                  temp_all_plots$data_plot, 
                  R0_from_seasonal_all_plots$model_plot,
                  R0_from_seasonal_all_plots$data_plot)

plot_grid2 = cowplot::plot_grid(plotlist = plot_list2, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B', 'C'))

title2 <- cowplot::ggdraw() +
  cowplot::draw_label("GAMs adjusted for all hypothesized confounding variables", fontface = 'bold', x = 0.5, hjust = 0.5, size = 12) 

plot_grid2 <- cowplot::plot_grid(title2, plot_grid2, ncol = 1, rel_heights = c(0.05, 1))

plot_grid2

ggsave(".../confounder_gam_all.tiff", 
       plot = plot_grid2,
       bg = "white",
       width = 8, height = 6, 
       dpi = 300, units = "in", compression = "lzw")
