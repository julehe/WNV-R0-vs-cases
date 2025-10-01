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

positive_obs <- df_wnnd_filtered %>% filter(incidence_per_100K>0)

thresh <- quantile(positive_obs$incidence_per_100K, 0.95)
thresh

df_wnnd_filtered <- filter(df_wnnd_filtered,incidence_per_100K <= thresh)

target_col <- which(colnames(df_wnnd_filtered) == "incidence_per_100K")
predictor_col <- which(colnames(df_wnnd_filtered) %in% 
                         c("temp_3mo_avg", "r0_from_seasonal_3mo_avg"))

#-------------------------------------------------------------------------------
# Load GAM analysis, plotting, and permutation test function 
#-------------------------------------------------------------------------------

source(".../helper scripts/gam_analysis.R")
source(".../helper scripts/plot_gam_f.R")
source(".../helper scripts/kendall_pairwise_permutation_test.R")

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

kendall_tau_temp <- round(cor.test(df_wnnd_filtered$temp_3mo_avg,
                                   df_wnnd_filtered$incidence_per_100K,
                                   method = "kendall")$estimate,3)
print(kendall_tau_temp)

label = paste0("tau[B] == ",kendall_tau_temp)

temp_plots$data_plot <- temp_plots$data_plot + 
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
  )  +  
  annotate("text", x = 10, 
           y = max(df_wnnd_filtered$incidence_per_100K, na.rm = TRUE),
           label = label, hjust = 0, size = 4, parse = T)
temp_plots[[3]]

### 2. R0 from "seasonal" NUTS3 temperature as predictor, incidence as outcome

R0_from_seasonal_analysis <- gam_analysis(df_wnnd_filtered,
                                          "r0_from_seasonal_3mo_avg")
R0_from_seasonal_plots <- plot_gam_f(df_wnnd_filtered,
                                     R0_from_seasonal_analysis$model,
                                     "r0_from_seasonal_3mo_avg", r0_predictor = T)

R0_from_seasonal_plots$model_plot <- R0_from_seasonal_plots$model_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) 

kendall_tau_R0 <- round(cor.test(df_wnnd_filtered$r0_from_seasonal_3mo_avg,
                                 df_wnnd_filtered$incidence_per_100K,
                                 method = "kendall")$estimate,3)
print(kendall_tau_R0)

label = paste0("tau[B] == ",kendall_tau_R0)

R0_from_seasonal_plots$data_plot <- R0_from_seasonal_plots$data_plot +
  labs(  
    x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
  )  +  
  scale_x_continuous(limits = c(0,1)) +  
  annotate("text", x = 0, 
           y = max(df_wnnd_filtered$incidence_per_100K, na.rm = TRUE),
           label = label, hjust = 0, size = 4, parse = T)

#-------------------------------------------------------------------------------
# Test significance of difference between rank correlations
#-------------------------------------------------------------------------------

res <- kendall_pairwise_permutation_test(
  df = df_wnnd_filtered,
  target_col = target_col,
  predictor_cols = predictor_col,
  nsim = 10000,
  two_sided = T,
  seed = 123
)

print(res)

#-------------------------------------------------------------------------------
# Save ggplots
#-------------------------------------------------------------------------------

plot_list1 = list(temp_plots$model_plot, temp_plots$data_plot, 
                  R0_from_seasonal_plots$model_plot, 
                  R0_from_seasonal_plots$data_plot)

plot_grid1 = cowplot::plot_grid(plotlist = plot_list1, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B', 
                                                                    'C', 'D'))

title1 <- cowplot::ggdraw() +
  cowplot::draw_label("GAMs fitted to data without upper 5% incidence observations", fontface = 'bold', x = 0.5, hjust = 0.5, size = 12) 

plot_grid1 <- cowplot::plot_grid(title1, plot_grid1, ncol = 1, rel_heights = c(0.05, 1))

plot_grid1

ggsave(".../main_gam_monthly_plot_wo_outliers.tiff", 
       plot = plot_grid1,
       bg = "white",
       width = 8, height = 6, 
       dpi = 300, units = "in", compression = "lzw")
