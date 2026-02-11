# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(mgcv)
library(tidyverse)

#-------------------------------------------------------------------------------
# Define GAM analysis function 
#-------------------------------------------------------------------------------

gam_analysis <- function(df, predictor, adjustments = c(NA), 
                         binary = F, k = 20){
  if(any(is.na(adjustments))){
    formula_adjustment <- ""
    print("No adjustments used")
  }else{
    formula_adjustment <- paste0("+",paste0(adjustments, collapse = "+")) 
  }
    
  if(binary == T){
    formula <- as.formula(
      paste(
        "wnnd_pa ~ s(",
        predictor,
        ", bs='cr', k=k)",
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
        ", bs='cr', k=k)",
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
  model_selected <- model
  
  # some model checks
  par(mfrow = c(2, 2))
  gam.check(model_selected, rep=500)
  print(k.check(model_selected))
  print(summary(model_selected))
  plot(model_selected)
  print(model_selected$outer.info)
  print(sum(influence(model_selected)))
  
  return(list(model = model_selected,
              summary = summary(model_selected)))
}