# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(tidyverse)
library(brms)
library(ggplot2)
library(cowplot)

#-------------------------------------------------------------------------------
# Load monthly aggregated data already matched with covariates
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
# BRMS model
#-------------------------------------------------------------------------------

# fit_zinb1 <- brm(wnnd_sum ~ s(temp_3mo_avg, bs="cr", k=15) + offset(log(population)),
#                  data = df_wnnd_filtered,
#                  family = zero_inflated_negbinomial(link = "log"),
#                  chains = 4, cores = 4,
#                  seed = 123,
#                  control = list(adapt_delta = 0.99, max_treedepth = 12))

#saveRDS(fit_zinb1, ".../zero_inflated_fit.rds")
fit_zinb1 <- readRDS(".../zero_inflated_fit.rds")

summary(fit_zinb1)
#post <- posterior_samples(fit_zinb1, pars = "b_zi_Intercept")
#prob_zi <- plogis(post$b_zi_Intercept)
#mean(prob_zi)
plot(fit_zinb1)
zi_est <- posterior_summary(fit_zinb1, variable = "zi")
zi_est
pp_check(fit_zinb1, type = "bars", ndraws = 100) +
  coord_cartesian(xlim = c(0, 5))
bayes_R2(fit_zinb1)

ce_df <- conditional_effects(fit_zinb1, effects = "temp_3mo_avg")$temp_3mo_avg

n_bins = 30

# plot model prediction with binned data
plot_gam_incidence_brms <- ggplot() + 
  stat_summary_bin(data = df_wnnd_filtered, 
                   aes(x = temp_3mo_avg, y = incidence_per_100K, color='1Data'), 
                   fun.data = function(y) {
                     mean_y <- mean(y)
                     sd_y <- sd(y)/sqrt(length(y))
                     return(data.frame(y = mean_y, ymin = mean_y - sd_y, ymax = mean_y + sd_y))
                   }, 
                   bins = n_bins, 
                   linewidth = 0.5,
                   geom = 'errorbar') +
  stat_summary_bin(data = df_wnnd_filtered,
                   aes(x = temp_3mo_avg, y = incidence_per_100K, color='1Data'),
                   fun = 'mean', bins = n_bins,
                   size = 1, geom='point') +
  geom_ribbon(
    data = ce_df,
    aes(x = temp_3mo_avg, ymin = 100000*(lower__/population), 
        ymax = 100000*(upper__/population), color = "2GAM Model"),
    fill = "red",
    alpha = 0.3,
    linetype = 0
  ) +
  geom_line(
    data = ce_df,
    aes(x = temp_3mo_avg, y = 100000*(estimate__/population), color = "2GAM Model"),
    linewidth = 0.7,
    linetype = "solid"
  ) +
  geom_rect(data = data.frame(xmin = 23.0, xmax = 25.8, ymin = 0, ymax = Inf, label = "R0 Optimum"),  
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = label),  
            fill = "#56B4E9",
            alpha = 0.3,
            linetype = 0) + 
  geom_segment(data = data.frame(x = 24.4, label = "R0 Optimum"),  
               aes(x = x, xend = x, y = 0, yend = Inf, color = label),  
               linewidth = 0.7, 
               linetype = "dashed") + 
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)")),  
    y = "Monthly WNND per 100K"  
  ) +  
  scale_color_manual(
    values = c("1Data" = "black","2GAM Model" = "red", "R0 Optimum" = "#56B4E9"),  
    labels = c("Data mean ± SE" ,"GAM fit", expression(R[0]^{"rel"} * " optimum"))
  ) +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +
  theme_bw(base_size=12)+
  guides() + 
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.75),
    panel.grid = element_blank()
  )

# Plot data
plot_gam_incidence_data <- ggplot() +
  geom_point(
    data = df_wnnd_filtered,
    aes(x = temp_3mo_avg, y = incidence_per_100K, color='Data'),
    shape = 19,
    size = 1,
    alpha = 0.2
  ) +
  labs(  
    x = expression(paste("3-month moving avg. temperature (", degree, "C)")),  
    y = "Monthly WNND per 100K"  
  ) +  
  scale_x_continuous(limits = c(10,NA),
                     breaks = seq(10, 28, by=2)) +  
  scale_color_manual(
    values = c("Data" = "black"),  
    labels = c("Data")
  ) + 
  theme_bw(base_size=12)+
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.75),
    panel.grid = element_blank()
  ) 

plot_list = list(plot_gam_incidence_brms, 
                 plot_gam_incidence_data)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=2, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B'))

title <- cowplot::ggdraw() +
  cowplot::draw_label("GAM fitted with zero-inflated distribution", fontface = 'bold', x = 0.5, hjust = 0.5, size = 12) 

plot_grid <- cowplot::plot_grid(title, plot_grid, ncol = 1, rel_heights = c(0.05, 1))

plot_grid

ggsave(".../brms_gam_monthly.tiff", 
       plot = plot_grid,
       width = 8, height = 3, 
       bg = "white",
       dpi = 300, units = "in", compression = "lzw")
