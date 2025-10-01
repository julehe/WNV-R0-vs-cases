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

df_wnnd_filtered$country <- substr(df_wnnd_filtered$NUTS_ID, 1, 2)
df_wnnd_filtered$country <- factor(df_wnnd_filtered$country)

# Exploratory Visual
plot_gam_incidence_data <- ggplot(data = df_wnnd_filtered,
                                  aes(x = temp_3mo_avg, y = incidence_per_100K, 
                                      color=country, group=country)) +
  geom_point(
    shape = 19,
    size = 1,
    alpha = 0.2,
    col="black"
  ) +
  geom_smooth(
  ) +
  geom_vline(xintercept = 24.4, linetype = "dashed", color = "red") +
  facet_wrap(~ country) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)")),  
    y = "Monthly WNND per 100K"  
  ) +  
  coord_cartesian(xlim = c(10, 29), ylim = c(0, 
                                             10)) +
  theme_bw(base_size=12)+
  theme(
    legend.position =  "none",
    legend.position.inside = c(0.2, 0.75),
    panel.grid = element_blank()
  ) 
plot_gam_incidence_data

# check data above 24.5
df_high_temps <- df_wnnd_filtered %>% filter(temp_3mo_avg > 24.5)
df_high_temps_greece <- df_high_temps %>% filter(country == "EL")
nrow(df_high_temps_greece)/nrow(df_high_temps)

target_col <- which(colnames(df_wnnd_filtered)== "incidence_per_100K")
predictor_col <- which(colnames(df_wnnd_filtered) %in% 
                         c("temp_3mo_avg", "r0_from_seasonal_3mo_avg"))

countries <- c("EL", "IT", "RO", "HU")

#-------------------------------------------------------------------------------
# Load GAM analysis, plotting, and permutation test function 
#-------------------------------------------------------------------------------

source(".../helper scripts/gam_analysis.R")
source(".../helper scripts/plot_gam_f.R")
source(".../helper scripts/kendall_pairwise_permutation_test.R")

#-------------------------------------------------------------------------------
# Loop countries and 
# Apply GAM analysis and plotting function to temperature/R0 covariates
#-------------------------------------------------------------------------------

plot_list_temp <- list()
plot_list_temp_data <- list()
plot_list_r0 <- list()
plot_list_r0_data <- list()
kendall_tau_list <- list()
counter = 1

for(country_label in countries){
  print(country_label)
  
  df_wnnd_filtered_temp <- df_wnnd_filtered %>% filter(country == country_label)
  
  ### 1. Temperature as predictor, incidence as outcome
  
  temp_analysis <- gam_analysis(df_wnnd_filtered_temp, 
                                "temp_3mo_avg")
  temp_plots <- plot_gam_f(df_wnnd_filtered_temp, temp_analysis$model, 
                           "temp_3mo_avg", n_bins = 25)
  
  temp_plots$model_plot <- temp_plots$model_plot +
    labs(  
      x = expression(paste("3-month moving avg. temperature (", degree, "C)"))
    )  +  
    scale_x_continuous(limits = c(10,NA),
                       breaks = seq(10, 28, by=2)) 
  
  kendall_tau_temp <- round(cor.test(df_wnnd_filtered_temp$temp_3mo_avg,
                                     df_wnnd_filtered_temp$incidence_per_100K,
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
             y = max(df_wnnd_filtered_temp$incidence_per_100K, na.rm = TRUE),
             label = label, hjust = 0, size = 4, parse = T)
  
  print(temp_plots$model_plot)
  print(temp_plots$data_plot)
  print(temp_plots$optimum)
  plot_list_temp[[counter]] <- temp_plots$model_plot
  plot_list_temp_data[[counter]] <- temp_plots$data_plot
  
  ### 2. R0 from "seasonal" NUTS3 temperature as predictor, incidence as outcome
  
  R0_from_seasonal_analysis <- gam_analysis(df_wnnd_filtered_temp,
                                            "r0_from_seasonal_3mo_avg")
  R0_from_seasonal_plots <- plot_gam_f(df_wnnd_filtered_temp,
                                       R0_from_seasonal_analysis$model,
                                       "r0_from_seasonal_3mo_avg", 
                                       r0_predictor = T, n_bins = 25)
  
  R0_from_seasonal_plots$model_plot <- R0_from_seasonal_plots$model_plot +
    labs(  
      x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
    )  +  
    scale_x_continuous(limits = c(0,1)) 
  
  kendall_tau_R0 <- round(cor.test(df_wnnd_filtered_temp$r0_from_seasonal_3mo_avg,
                                   df_wnnd_filtered_temp$incidence_per_100K,
                                   method = "kendall")$estimate,3)
  print(kendall_tau_R0)
  
  label = paste0("tau[B] == ",kendall_tau_R0)
  
  R0_from_seasonal_plots$data_plot <- R0_from_seasonal_plots$data_plot +
    labs(  
      x = expression(paste("3-month moving avg. " * R[0]^{"rel, s"}))
    )  +  
    scale_x_continuous(limits = c(0,1)) +  
    annotate("text", x = 0, 
             y = max(df_wnnd_filtered_temp$incidence_per_100K, na.rm = TRUE),
             label = label, hjust = 0, size = 4, parse = T)
  
  print(R0_from_seasonal_plots$model_plot)
  print(R0_from_seasonal_plots$data_plot)
  plot_list_r0[[counter]] <- R0_from_seasonal_plots$model_plot
  plot_list_r0_data[[counter]] <- R0_from_seasonal_plots$data_plot
  
  # Test significance of difference between rank correlations
  kendall_taus <- kendall_pairwise_permutation_test(
    df = df_wnnd_filtered_temp,
    target_col = target_col,
    predictor_cols = predictor_col,
    nsim = 10000,
    two_sided = T,
    seed = 123
  )
  
  print(kendall_taus)
  
  kendall_tau_list[counter] <- kendall_taus
  
  counter = counter + 1 
}

#-------------------------------------------------------------------------------
# Save ggplots
#-------------------------------------------------------------------------------

grid_list <- list()
country_names <- c("Greece", "Italy", "Romania", "Hungary")

for(i in 1:4){
  plots = list(plot_list_temp[[i]], 
               plot_list_temp_data[[i]], 
               plot_list_r0[[i]],
               plot_list_r0_data[[i]])
  
  plot_grid = cowplot::plot_grid(plotlist = plots, ncol=2, label_size = 12,
                                  align = "h", axis = "b", labels = c('A', 'B', 
                                                                      'C', 'D'))
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(country_names[i], fontface = 'bold', x = 0.5, hjust = 0.5, size = 12) 
  
  plot_grid <- cowplot::plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05, 1))
  
  print(plot_grid)
  
  grid_list[[i]] <- plot_grid
}

ggsave(".../gams_greece.tiff",
       plot = grid_list[[1]],
       bg = "white",
       width = 8, height = 6,
       dpi = 300, units = "in", compression = "lzw")

ggsave(".../gams_italy.tiff",
       plot = grid_list[[2]],
       bg = "white",
       width = 8, height = 6,
       dpi = 300, units = "in", compression = "lzw")

ggsave(".../gams_romania.tiff",
       plot = grid_list[[3]],
       bg = "white",
       width = 8, height = 6,
       dpi = 300, units = "in", compression = "lzw")

ggsave(".../gams_hungary.tiff",
       plot = grid_list[[4]],
       bg = "white",
       width = 8, height = 6,
       dpi = 300, units = "in", compression = "lzw")

