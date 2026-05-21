# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
#
# Description: This provides a helper function for computing posterior statistics
# of the optimum of the shared smooth and country-specific smooths
#
# ==============================================================================

library(mgcv)
library(tidyverse)
library(MASS)

#-------------------------------------------------------------------------------
# Define GAM optimum calculator
#-------------------------------------------------------------------------------

compute_optimum <- function(df, model, predictor, 
                            country_id = NA) {

  set.seed(123)
  # optimum of shared smooth
  pd_shared <- data.frame(
    seq(10, 40, length = 1000),
    country = "EL",
    # values of adjustment variables don't matter since they are specified as simple linear fixed effects
    urban_fabric = median(df$urban_fabric),
    rain_3mo_cum_pop_weight = median(df$rain_3mo_cum_pop_weight),
    percent_65plus = median(df$percent_65plus),
    arable_land = median(df$arable_land),
    gdp_per_capita = median(df$gdp_per_capita)
  )
  names(pd_shared)[1] <- c(predictor)
  
  
  Xp_shared <- predict(model, pd_shared, type = "lpmatrix",
                       exclude = c("s(country,temp_3mo_avg_pop_weight)"))
  
  beta <- coef(model)
  Vb <- vcov(model)
  
  n <- 5000
  br <- mvrnorm(n, beta, Vb)
  x_opt_shared <- numeric(n)
  
  for (i in seq_len(n)) {
    f_i <- Xp_shared %*% br[i, ]
    x_opt_shared[i] <- pd_shared[,1][which.max(f_i)]
  }
  
  # sanity check
  if(any(x_opt_shared == 40)){
    print(paste(sum(x_opt_shared == 40), "Samples of shared temperature optimum > 40°C"))
  }
  
  CI_optimal <- quantile(x_opt_shared, c(0.025, 0.5, 0.975))
  CI_optimal["mean"] <- mean(x_opt_shared)
  
  if(is.na(country_id)){
    CI_optimal_country <- NA
  }else{
    
    set.seed(123)
    # optimum of country smooth
    pd_country <- data.frame(
      seq(10, 40, length = 1000),
      country = country_id,
      # values of adjustment variables don't matter since they are specified as simple linear fixed effects
      urban_fabric = median(filter(df, country == country_id)$urban_fabric),
      rain_3mo_cum_pop_weight = median(filter(df, country == country_id)$rain_3mo_cum_pop_weight),
      percent_65plus = median(filter(df, country == country_id)$percent_65plus),
      arable_land = median(filter(df, country == country_id)$arable_land),
      gdp_per_capita = median(filter(df, country == country_id)$gdp_per_capita)
    )
    names(pd_country)[1] <- predictor
    
    Xp_country <- predict(model, pd_country, type = "lpmatrix")
    
    x_opt_country <- numeric(5000)
    
    for (i in seq_len(5000)) {
      f_i <- Xp_country %*% br[i, ]
      x_opt_country[i] <- pd_country[[predictor]][which.max(f_i)]
    }
    
    # sanity check 
    if(any(x_opt_country == 40)){
      print(paste(sum(x_opt_country == 40), "Samples of country temperature optimum > 40°C"))
    }
    
    CI_optimal_country <- quantile(x_opt_country, c(0.025, 0.5, 0.975))
    CI_optimal_country["mean"] <- mean(x_opt_country)
  }
  
  return(list(shared_optimum = CI_optimal,
              country_optimum = CI_optimal_country
              ))
  
}
