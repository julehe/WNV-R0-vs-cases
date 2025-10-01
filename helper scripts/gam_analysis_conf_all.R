# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(mgcv)
library(tidyverse)

#-------------------------------------------------------------------------------
# Define GAM analysis function 
#-------------------------------------------------------------------------------

gam_analysis_conf_all <- function(df, predictor, binary = F, upper_k = 20){
  # try different k
  gams_list <- list()
  for(k in c(3:upper_k)){
    
    if(binary == T){
      
      formula <- as.formula(
        paste(
          "wnnd_pa ~ s(",
          predictor,
          ", bs='cr', k=k) + urban_fabric + rain_3mo_cum", 
          "+ percent_65plus + perm_irrigated + gdp_per_capita", sep=""
        )
      )
      model <- gam(formula,
                   family = binomial(), 
                   data = df,
                   method = "REML")
      gams_list[[as.character(k)]] <- model 
      print(k)
      
    }else{
      
      formula <- as.formula(
        paste(
          "wnnd_sum ~ s(",
          predictor,
          ", bs='cr', k=k) + urban_fabric + rain_3mo_cum", 
          "+ percent_65plus + perm_irrigated + gdp_per_capita", sep=""
        )
      )
      
      model <- gam(formula, 
                   family = nb(), 
                   offset = log(population), 
                   data = df,
                   method="REML")
      gams_list[[as.character(k)]] <- model 
      print(k)
    }
  }
  # select model based on AIC
  aics <- sapply(gams_list, AIC)
  print(aics)
  # Get local minimum, i.e., k where k+1 will not decrease AIC
  k_selected <- names(aics)[which(diff(aics) > 0)[1]] 
  print(k_selected)
  model_selected <- gams_list[[k_selected]] 
  
  # some model checks
  par(mfrow = c(2, 2))
  gam.check(model_selected, rep=500)
  print(k.check(model_selected))
  print(summary(model_selected))
  plot(model_selected)
  print(model_selected$outer.info)
  print(sum(influence(model_selected)))
  
  return(list(model = model_selected,
              k = k_selected,
              aics = aics,
              summary = summary(model_selected)))
}