# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
#
# Description: This provides a helper function for running fitting the GAMs and
# retrieving results
#
# ==============================================================================

library(mgcv)
library(tidyverse)

#-------------------------------------------------------------------------------
# Define GAM analysis function 
#-------------------------------------------------------------------------------

gam_analysis <- function(df, predictor, adjustments = c(NA), 
                         binary = F, k_shared = 20, k_country = 5,
                         country_effects = T){
  if(any(is.na(adjustments))){
    formula_adjustment <- ""
    print("No adjustments used")
  }else{
    formula_adjustment <- paste0("+",paste0(adjustments, collapse = "+")) 
  }
  
  if(country_effects == T){
    formula_country <- paste0("+ s(country,", 
                              predictor, 
                              ", bs='sz', k=k_country, xt = list(bs='cr'))")
  }else{
    formula_country <- ""
  }
    
  if(binary == T){
    formula <- as.formula(
      paste(
        "wnnd_pa ~ s(",
        predictor,
        ", bs='cr', k=k_shared)",
        formula_country,
        formula_adjustment,
        sep=""
      )
    )
    
    model <- gam(formula,
                 family = binomial(), 
                 data = df,
                 method = "REML")
    
  }else{
    
    formula <- as.formula(
      paste(
        "wnnd_sum ~ s(",
        predictor,
        ", bs='cr', k=k_shared)",
        formula_country,
        formula_adjustment,
        sep=""
      )
    )
    
    model <- gam(formula, 
                 family = nb(), 
                 offset = log(population), 
                 data = df,
                 select = T,
                 method="REML")
  }
  
  print(AIC(model))
  
  # some model checks
  par(mfrow = c(2, 2))
  gam.check(model, rep=500)
  print(k.check(model))
  print(summary(model))
  plot(model)
  print(model$outer.info)
  print(sum(influence(model)))
  
  return(list(model = model,
              summary = summary(model)))
}