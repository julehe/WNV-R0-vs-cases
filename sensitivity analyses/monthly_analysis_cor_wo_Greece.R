# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(tidyverse)
library(multcompView)
library(cowplot)
library(grid)
library(gridExtra)
source("code/helper scripts/kendall_pairwise_permutation_test.R")

#-------------------------------------------------------------------------------
# Load monthly data already matched with covariates
#-------------------------------------------------------------------------------

df_wnnd <- 
  read.csv("data/WNND_data_with_covariates.csv")

# only look at locations-years that have at least one case
df_wnnd <- df_wnnd %>%
  mutate(year = year(Date)) %>%
  group_by(NUTS_ID, year) %>%
  filter(any(wnnd_sum > 0)) %>%
  ungroup()

# filter out greece
df_wnnd <- df_wnnd %>%
  filter(!grepl("^EL", NUTS_ID))

# ppcor::pcor.test(df_wnnd$temp_3mo_avg_pop_weight, df_wnnd$incidence_per_100K, 
#                  df_wnnd$perm_irrigated, method = "kendall")
# 
# ppcor::pcor.test(df_wnnd$r0_from_seasonal_3mo_avg_pop_weight, df_wnnd$incidence_per_100K, 
#                  df_wnnd$perm_irrigated, method = "kendall")

# df_wnnd <- df_wnnd %>%
#   mutate(year = year(Date)) %>%
#   group_by(NUTS_ID) %>%
#   filter(any(wnnd_sum > 0)) %>%
#   ungroup()
# 
# 
# df_wnnd <- df_wnnd %>%
#   mutate(year = year(Date)) %>%
#   group_by(NUTS_ID) %>%
#   filter(any(wnnd_sum > 0)) %>%
#   mutate(first_year = min(year[wnnd_sum > 0], na.rm = TRUE)) %>%
#   filter(year >= first_year) %>%
#   ungroup()

# df_wnnd <- filter(df_wnnd, temp_3mo_avg>15)

# seasonal
# df_wnnd <- df_wnnd %>%
#   #filter(temp_3mo_avg>15) %>%
#   dplyr::select(NUTS_ID, Date, incidence_per_100K, wnnd_sum, temp_3mo_avg,
#                 r0_from_seasonal_3mo_avg, r0_from_daily_native_3mo_avg) %>%
#   mutate(year = year(Date)) %>%
#   group_by(NUTS_ID, year) %>%
#   filter(any(wnnd_sum > 0)) %>%
#   mutate(incidence_per_100K = incidence_per_100K-mean(incidence_per_100K),
#          temp_3mo_avg = temp_3mo_avg-mean(temp_3mo_avg),
#          r0_from_daily_native_3mo_avg = r0_from_daily_native_3mo_avg-mean(r0_from_daily_native_3mo_avg),
#          r0_from_seasonal_3mo_avg = r0_from_seasonal_3mo_avg-mean(r0_from_seasonal_3mo_avg)) %>%
#   ungroup() 
  

# spatial
# df_wnnd_x <- df_wnnd %>%
#   mutate(year = year(Date),
#          month = month(Date)) %>%
#   group_by(NUTS_ID, year) %>%
#   summarise(incidence_per_100K = sum(wnnd_sum)/mean(population),
#          temp_3mo_avg = mean(temp_3mo_avg[month %in% 7:9], na.rm=T),
#          r0_from_seasonal_native_3mo_avg = mean(r0_from_seasonal_native_3mo_avg[month %in% 7:9], na.rm=T),
#          r0_from_seasonal_3mo_avg = mean(r0_from_seasonal_3mo_avg[month %in% 7:9], na.rm=T),
#          .groups = "drop") %>%
#   group_by(NUTS_ID) %>%
#   summarise(incidence_per_100K = mean(incidence_per_100K),
#             temp_3mo_avg = mean(temp_3mo_avg),
#             r0_from_seasonal_native_3mo_avg = mean(r0_from_seasonal_native_3mo_avg),
#             r0_from_seasonal_3mo_avg = mean(r0_from_seasonal_3mo_avg),
#             .groups = "drop") %>%
#   mutate(incidence_per_100K = incidence_per_100K-mean(incidence_per_100K, na.rm = T),
#          temp_3mo_avg = temp_3mo_avg-mean(temp_3mo_avg, na.rm = T),
#          r0_from_seasonal_native_3mo_avg = r0_from_seasonal_native_3mo_avg-mean(r0_from_seasonal_native_3mo_avg, na.rm = T),
#          r0_from_seasonal_3mo_avg = r0_from_seasonal_3mo_avg-mean(r0_from_seasonal_3mo_avg, na.rm = T))

# yearly --> weird looks like when looking at spatial and interannual signal 
# together there is higher rank correlation with R0 but not or extremely marginal when separate...?
# what happens if I leave spatial signal in the seasonal? this is weird...
# df_wnnd_x <- df_wnnd %>%
#   mutate(year = year(Date),
#          month = month(Date)) %>%
#   group_by(NUTS_ID, year) %>%
#   summarise(incidence_per_100K = sum(wnnd_sum)/mean(population),
#             temp_3mo_avg = mean(temp_3mo_avg[month %in% 7:9], na.rm=T),
#             r0_from_seasonal_native_3mo_avg = mean(r0_from_seasonal_native_3mo_avg[month %in% 7:9], na.rm=T),
#             r0_from_seasonal_3mo_avg = mean(r0_from_seasonal_3mo_avg[month %in% 7:9], na.rm=T),
#             .groups = "drop") %>%
#   group_by(NUTS_ID) %>%
#   filter(any(incidence_per_100K>0)) %>%
#   mutate(incidence_per_100K = incidence_per_100K-mean(incidence_per_100K),
#          temp_3mo_avg = temp_3mo_avg-mean(temp_3mo_avg),
#          r0_from_seasonal_native_3mo_avg = r0_from_seasonal_native_3mo_avg-mean(r0_from_seasonal_native_3mo_avg),
#          r0_from_seasonal_3mo_avg = r0_from_seasonal_3mo_avg-mean(r0_from_seasonal_3mo_avg)) %>%
#   ungroup()

colnames(df_wnnd)

# Keep only relevant columns
df_wnnd <- df_wnnd[,c(which(grepl("pop_weight", colnames(df_wnnd)) & 
                              (
                                grepl("r0", colnames(df_wnnd)) |
                                  grepl("temp", colnames(df_wnnd))
                              )
                            ),
                      ncol(df_wnnd))]

# We will also apply the analysis on a filtered dataset with temps>15°C
df_wnnd_upper_temps <- df_wnnd %>% 
  filter(temp_3mo_avg_pop_weight>15)

#-------------------------------------------------------------------------------
# Calculate rank correlations
#-------------------------------------------------------------------------------

predictors <- names(df_wnnd)[c(1:26)]

# this order is just for better overview when printing results
order_results <- c("temp_3mo_avg_pop_weight", "r0_from_daily_3mo_avg_pop_weight",
                   "r0_from_monthly_3mo_avg_pop_weight", "r0_from_seasonal_3mo_avg_pop_weight",
                   "r0_from_daily_native_3mo_avg_pop_weight", "r0_from_monthly_native_3mo_avg_pop_weight",
                   "r0_from_seasonal_native_3mo_avg_pop_weight",
                   "temp_2mo_avg_pop_weight", "r0_from_daily_2mo_avg_pop_weight",
                   "r0_from_monthly_2mo_avg_pop_weight", "r0_from_seasonal_2mo_avg_pop_weight",
                   "r0_from_daily_native_2mo_avg_pop_weight", "r0_from_monthly_native_2mo_avg_pop_weight",
                   "r0_from_seasonal_native_2mo_avg_pop_weight",
                   "temp_4mo_avg_pop_weight", "r0_from_daily_4mo_avg_pop_weight",
                   "r0_from_monthly_4mo_avg_pop_weight", "r0_from_seasonal_4mo_avg_pop_weight",
                   "r0_from_daily_native_4mo_avg_pop_weight", "r0_from_monthly_native_4mo_avg_pop_weight",
                   "r0_from_seasonal_native_4mo_avg_pop_weight",
                   "temp_pop_weight", "r0_from_daily_pop_weight",
                   "r0_from_monthly_pop_weight", 
                   "r0_from_daily_native_pop_weight", "r0_from_monthly_native_pop_weight")

# calculate kendall tau's on "full" dataset
kendalls_tauB <- sapply(c(1:26), 
                        function(x) cor.test(pull(df_wnnd, 27),
                                             pull(df_wnnd, x),
                                             method = "kendall")$estimate)

names(kendalls_tauB) <- predictors
kendalls_tauB[order_results]

# calculate kendall tau's on restricted dataset
kendalls_tauB_restricted <- sapply(c(1:26), 
                                   function(x) cor.test(pull(df_wnnd_upper_temps, 27), 
                                                        pull(df_wnnd_upper_temps, x),
                                                        method = "kendall")$estimate)

names(kendalls_tauB_restricted) <- predictors
kendalls_tauB_restricted[order_results]

#-------------------------------------------------------------------------------
# Test significance of difference between rank correlations
#-------------------------------------------------------------------------------

# keep only 3-month moving average values of predictors
df_wnnd <- df_wnnd[,c(
  which(grepl("3mo", colnames(df_wnnd))), 
                      27)
  ]
#df_wnnd <- df_wnnd[,c(2,6,10,14,19,22,24,27)]
df_wnnd_upper_temps <- df_wnnd_upper_temps[,c(
  which(grepl("3mo", colnames(df_wnnd_upper_temps))), 
  27)
]
#df_wnnd_upper_temps <- df_wnnd_upper_temps[,c(2,6,10,14,19,22,24,27)]
predictors <- names(df_wnnd)[c(1:7)]

# loop through the full and restricted dataset and store results in lists
dfs <- list(df_wnnd, df_wnnd_upper_temps)
res <- list()

for(i in 1:2){
  # apply pairwise permutation test
  res[[i]] <- kendall_pairwise_permutation_test(
    df = dfs[[i]],
    target_col = 8,                    
    predictor_cols = c(1:7),     
    nsim = 10000,
    two_sided = T,
    seed = 123
  )
  
  # add tau's and their p-values 
  pred_df <- dfs[[i]][, c(1:7)]
  res[[i]]$taus <- lapply(predictors, \(p) #\(p) is equivalent to function(p){}
                          tibble(
                            predictor = p,
                            tau = cor.test(dfs[[i]][[8]], 
                                           rank(pred_df[[p]]), 
                                           method = "kendall")$estimate,
                            p = cor.test(dfs[[i]][[8]], 
                                         rank(pred_df[[p]]), 
                                         method = "kendall")$p.value
                          )) |> bind_rows()
  
  # build matrix of p-values of pairwise difference of tau's
  p_mat <- matrix(1, nrow = length(predictors), ncol = length(predictors),
                  dimnames = list(predictors, predictors))
  for (j in seq_len(nrow(res[[i]]$comparisons))) {
    row <- res[[i]]$comparisons$predictor1[j]
    col <- res[[i]]$comparisons$predictor2[j]
    p_val <- res[[i]]$comparisons$p_value[j]
    p_mat[row, col] <- p_val
    p_mat[col, row] <- p_val
  }

  # generate group letters based on statistically non-significant differences
  group_letters <- multcompLetters(p_mat, threshold = 0.05)$Letters

  res[[i]]$taus <- res[[i]]$taus %>%
    mutate(group = group_letters[predictor]) %>%
    arrange(-desc(tau))
  
  res[[i]]$taus$predictor <- factor(res[[i]]$taus$predictor,
                                    levels = res[[i]]$taus$predictor)
  
}

# some ordering for plotting
predictor_order <- c(
  "temp_3mo_avg_pop_weight",
  "r0_from_daily_native_3mo_avg_pop_weight",
  "r0_from_monthly_native_3mo_avg_pop_weight",
  "r0_from_seasonal_native_3mo_avg_pop_weight",
  "r0_from_daily_3mo_avg_pop_weight",
  "r0_from_monthly_3mo_avg_pop_weight",
  "r0_from_seasonal_3mo_avg_pop_weight"
)

# Convert predictor to factor with specified order
res[[1]]$taus$predictor <- factor(res[[1]]$taus$predictor,
                                  levels = predictor_order)
res[[2]]$taus$predictor <- factor(res[[2]]$taus$predictor,
                                  levels = predictor_order)

# plot for results on "full" dataset
plot1 <- ggplot(res[[1]]$taus, aes(x = predictor, y = tau)) + 
  geom_col(fill = "#56B4E9", width = 0.7) + 
  geom_text(aes(label = round(tau, 3)), vjust = -0.5, size = 4) +
  geom_text(aes(label = group), vjust = -2, size = 4) +
  labs(
    title = "\"Full\" dataset"
  ) +
  scale_x_discrete(labels = c(
    "temp_3mo_avg_pop_weight" = expression(paste("Temperature")),
    "r0_from_daily_native_3mo_avg_pop_weight" = expression(R["0, grid"]^{"rel, d"}),
    "r0_from_monthly_native_3mo_avg_pop_weight" = expression(R["0, grid"]^{"rel, m"}),
    "r0_from_seasonal_native_3mo_avg_pop_weight" = expression(R["0, grid"]^{"rel, s"}),
    "r0_from_daily_3mo_avg_pop_weight" = expression(R[0]^{"rel, d"}),
    "r0_from_monthly_3mo_avg_pop_weight" = expression(R[0]^{"rel, m"}),
    "r0_from_seasonal_3mo_avg_pop_weight" = expression(R[0]^{"rel, s"})
  )) + 
  scale_y_continuous(limits = c(0, 0.52), expand = c(0, 0)) + 
  theme_bw(base_size = 14) +  
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) 

plot1

# plot for results on "restricted" dataset
plot2 <- ggplot(res[[2]]$taus, aes(x = predictor, y = tau)) + 
  geom_col(fill = "#56B4E9", width = 0.7) + 
  geom_text(aes(label = round(tau, 3)), vjust = -0.5, size = 4) +
  geom_text(aes(label = group), vjust = -2, size = 4) +
  labs(
    title = "Limited to temperatures >15°C"
  ) +
  scale_x_discrete(labels = c(
    "temp_3mo_avg" = expression(paste("Temperature")),
    "r0_from_daily_native_3mo_avg" = expression(R["0, grid"]^{"rel, d"}),
    "r0_from_monthly_native_3mo_avg" = expression(R["0, grid"]^{"rel, m"}),
    "r0_from_seasonal_native_3mo_avg" = expression(R["0, grid"]^{"rel, s"}),
    "r0_from_daily_3mo_avg" = expression(R[0]^{"rel, d"}),
    "r0_from_monthly_3mo_avg" = expression(R[0]^{"rel, m"}),
    "r0_from_seasonal_3mo_avg" = expression(R[0]^{"rel, s"})
  )) + 
  scale_y_continuous(limits = c(0, 0.47), expand = c(0, 0)) + 
  theme_bw(base_size = 14) +  
  theme(
    plot.title = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) 

plot2

# Combine plots in grid
plot_list = list(plot1,
                 plot2)

plot_grid = cowplot::plot_grid(plotlist = plot_list, ncol=1, label_size = 12,
                                align = "h", axis = "b", labels = c('A', 'B'))

y.grob <- textGrob(expression("Kendall's " * tau[B] * " with monthly WNND per 100K"), 
                   gp=gpar(col="black", fontsize=12), rot=90)

x.grob <- textGrob("3-month moving average", 
                   gp=gpar(col="black", fontsize=12))

grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob))

ggsave("figures/monthly/kendalls_monthly.tiff", 
       plot = grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob)),
       width = 8, height = 7, 
       dpi = 300, units = "in", compression = "lzw")
ggsave("figures/monthly/kendalls_monthly.svg", 
       plot = grid.arrange(arrangeGrob(plot_grid, left = y.grob, bottom = x.grob)),
       width = 8, height = 7, 
       units = "in")
