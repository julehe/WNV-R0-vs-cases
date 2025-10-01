# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(tidyverse)
library(multcompView)
library(cowplot)
library(grid)
library(gridExtra)
source(".../helper scripts/kendall_pairwise_permutation_test.R")

#-------------------------------------------------------------------------------
# Load yearly data already matched with covariates
#-------------------------------------------------------------------------------

# The covariate data is openly available from:
# - Copernicus (https://cds.climate.copernicus.eu/datasets/derived-era5-land-daily-statistics?tab=overview)
# - Eurostat (https://ec.europa.eu/eurostat/de/)
# - CORINE (https://land.copernicus.eu/en/products/corine-land-cover)
# Restrictions apply to the Human WNND data which is only available upon request
# at TESSy (https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
df_wnnd <- read.csv("....csv")

# aggregate over years to "spatial-only" dataset 
df_wnnd <- df_wnnd %>%
  group_by(NUTS_ID) %>%
  summarise(wnnd_sum = sum(wnnd_sum),
            temp_summer = mean(temp_summer),
            r0_from_daily_summer_native = mean(r0_from_daily_summer_native),
            r0_from_monthly_summer_native = mean(r0_from_monthly_summer_native),
            r0_from_seasonal_summer_native = mean(r0_from_seasonal_summer_native),
            r0_from_daily_summer = mean(r0_from_daily_summer),
            r0_from_monthly_summer = mean(r0_from_monthly_summer),
            population = mean(population),
            r0_from_seasonal_summer = mean(r0_from_seasonal_summer),
            mean_incidence_per_100K = mean(incidence_per_100K))

df_wnnd$wnnd_sum <- round(df_wnnd$wnnd_sum)
df_wnnd$population <- round(df_wnnd$population)
# total incidence over the period
df_wnnd$incidence_per_100K <- 100000*df_wnnd$wnnd_sum/(df_wnnd$population)

colnames(df_wnnd)

# Keep only relevant columns
df_wnnd <- df_wnnd[,c(3:8,10,12)]

#-------------------------------------------------------------------------------
# Calculate rank correlations
#-------------------------------------------------------------------------------

predictors <- names(df_wnnd)[c(1:7)]

# calculate kendall tau's 
kendalls_tauB <- sapply(c(1:7), function(x) cor.test(pull(df_wnnd, 8), 
                                                     pull(df_wnnd, x), 
                                                     method = "kendall")$estimate)

names(kendalls_tauB) <- predictors
kendalls_tauB

#-------------------------------------------------------------------------------
# Test significance of difference between rank correlations
#-------------------------------------------------------------------------------

# apply pairwise permutation test
res <- kendall_pairwise_permutation_test(
  df = df_wnnd,
  target_col = 8,                    
  predictor_cols = c(1:7),     
  nsim = 10000,
  two_sided = T,
  seed = 123
)

# add tau's and their p-values 
pred_df <- df_wnnd[, c(1:7)]

res$taus <- lapply(predictors, \(p) #\(p) equivalent to function(p){}
  tibble(
    predictor = p,
    tau = cor.test(df_wnnd[[8]], rank(pred_df[[p]]), method = "kendall")$estimate,
    p = cor.test(df_wnnd[[8]], rank(pred_df[[p]]), method = "kendall")$p.value
  )) |> bind_rows()

# build matrix of p-values of pairwise difference of tau's
p_mat <- matrix(1, nrow = length(predictors), ncol = length(predictors),
                dimnames = list(predictors, predictors))
for (j in seq_len(nrow(res$comparisons))) {
  row <- res$comparisons$predictor1[j]
  col <- res$comparisons$predictor2[j]
  p_val <- res$comparisons$p_value[j]
  p_mat[row, col] <- p_val
  p_mat[col, row] <- p_val  
}

# generate group letters based on statistically non-significant differences
group_letters <- multcompLetters(p_mat, threshold = 0.05)$Letters

res$taus <- res$taus %>%
  mutate(group = group_letters[predictor]) %>%
  arrange(-desc(tau))

res$taus$predictor <- factor(res$taus$predictor,
                             levels = res$taus$predictor)

# print result
print(res)
