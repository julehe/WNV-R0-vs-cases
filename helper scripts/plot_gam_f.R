# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(mgcv)
library(tidyverse)
library(ggplot2)
library(MASS)

#-------------------------------------------------------------------------------
# Define GAM plotting function 
#-------------------------------------------------------------------------------

plot_gam_f <- function(df, model, predictor, adjustments = c(NA), 
                       r0_predictor = F, type = "monthly", binary = F, 
                       n_bins = 30){
  
  if(r0_predictor){
    lower_bound = 0
  }else{
    lower_bound = 10
  }
  
  new_data <- data.frame(
    seq(lower_bound,
        max(df[[predictor]], na.rm = TRUE),
        length.out = 1000)
  )

  names(new_data) <- predictor 
  
  if(all(!is.na(adjustments))){
    new_data[adjustments] <- lapply(model$model[adjustments], median)
  }
  
  # predict with GAM on new data
  prediction <- predict(
    model, se.fit = TRUE, newdata = new_data,
    unconditional = TRUE, type = "link"
  )
  
  # extract mean prediction and 95CI of that mean prediction
  new_data$mean_prediction <- prediction$fit
  new_data$prediction_upper_95CI <- new_data$mean_prediction + 1.96*prediction$se.fit
  new_data$prediction_lower_95CI <- new_data$mean_prediction - 1.96*prediction$se.fit 
  
  # apply inverse of link function to get to response scale
  new_data$mean_prediction  <- model$family$linkinv(new_data$mean_prediction)
  new_data$prediction_upper_95CI <- model$family$linkinv(new_data$prediction_upper_95CI)
  new_data$prediction_lower_95CI <- model$family$linkinv(new_data$prediction_lower_95CI)
  
  if(r0_predictor){
    optimal_pred <- NA
    CI_optimal <- NA
  }else{
    set.seed(123)
    #extract optimum predictor value of GAM mean
    optimal_pred <- new_data[[predictor]][which.max(new_data$mean_prediction)]
    
    # generate 95% CI of optimum predictor value 
    # (following https://webhomes.maths.ed.ac.uk/~swood34/mgcv/tampere/mgcv-advanced.pdf)
    pd <- data.frame(
        seq(10, 40, length = 1000)
    )
    names(pd) <- predictor 
    
    if(all(!is.na(adjustments))){
      pd[adjustments] <- lapply(model$model[adjustments], median, na.rm = TRUE)
    }
    
    Xp <- predict(model, pd, type = "lpmatrix")
    
    beta <- coef(model)
    Vb <- vcov(model)
    
    n <- 5000
    br <- mvrnorm(n, beta, Vb)
    x_opt <- numeric(n)
    
    for (i in seq_len(n)) {
      f_i <- Xp %*% br[i, ]
      x_opt[i] <- pd[,1][which.max(f_i)]
    }
    
    # sanity check 
    if(any(x_opt == 40)){
      print(paste(sum(x_opt == 40), "Samples of temperature optimum > 40°C"))
    }
    
    CI_optimal <- quantile(x_opt, c(0.025, 0.5, 0.975))
    CI_optimal["mean"] <- mean(x_opt)
  }
  
  
  if(binary == T){
    y <- "wnnd_pa"
    if(type == "monthly"){
      y_label <- "WNND presence probability" 
    }else if(type == "yearly"){
      y_label <- "Yearly WNND probability"
    }else{
      y_label <- "WNND 2010-2023"
    }
  }else{
    # scale predictions to incidence per 100K
    new_data$mean_prediction <- new_data$mean_prediction*100000
    new_data$prediction_upper_95CI <- new_data$prediction_upper_95CI*100000
    new_data$prediction_lower_95CI <- new_data$prediction_lower_95CI*100000
    y <- "incidence_per_100K"
    y_label <- 
    if(type == "monthly"){
      y_label <- "Monthly WNND per 100K"  
    }else if(type == "yearly"){
      y_label <- "Yearly WNND per 100K" 
    }else{
      y_label <- "Avg. yearly WNND per 100K" 
    }
  }
  
  # Plot GAM prediction with (binned) data
  if(binary == T){
    plot_gam <- ggplot() + 
      # data
      geom_point(
        data = df,
        aes(x = .data[[predictor]], y = .data[[y]], color = "1Data"),
        shape = 19,
        size = 1,
        alpha = 0.1
      ) 
  }else{
    plot_gam <- ggplot() +
      # binned data
      stat_summary_bin(data = df, 
                       aes(x = .data[[predictor]], y = .data[[y]], color='1Data'), 
                       fun.data = function(y) {
                         mean_y <- mean(y)
                         sd_y <- sd(y)/sqrt(length(y))
                         return(data.frame(y = mean_y, ymin = mean_y - sd_y, ymax = mean_y + sd_y))
                       }, 
                       bins = n_bins, 
                       linewidth = 0.5,
                       geom = 'errorbar') +
      stat_summary_bin(data = df,
                       aes(x = .data[[predictor]], y = .data[[y]], color='1Data'),
                       fun = 'mean', bins = n_bins,
                       size = 1, geom='point')
  }
  plot_gam <- plot_gam +
    # GAM prediction
    geom_ribbon(
      data = new_data,
      aes(x = .data[[predictor]], ymin = prediction_lower_95CI, 
          ymax = prediction_upper_95CI, color = "2GAM Model"),
      fill = "red",
      alpha = 0.3,
      linetype = 0
    ) +
    geom_line(
      data = new_data,
      aes(x = .data[[predictor]], y = mean_prediction, color = "2GAM Model"),
      linewidth = 0.7,
      linetype = "solid"
    ) +
    labs(  
      y = y_label
    ) +
    theme_bw(base_size=12)+
    guides() + 
    theme(
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.75),
      panel.grid = element_blank()
    )
  
  if(r0_predictor == T){
    plot_gam <- plot_gam +  
      scale_color_manual(
        values = c("1Data" = "black","2GAM Model" = "red"),  
        labels = c("Data mean ± SE" ,"GAM fit")
      )
  }else if(binary == T){
    plot_gam <- plot_gam +
      # Theoretical R0 optimum
      geom_rect(data = data.frame(xmin = 23.0, xmax = 25.8, ymin = 0, ymax = Inf, 
                                  label = "R0 Optimum"),  
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                    color = label),  
                fill = "#56B4E9",
                alpha = 0.3,
                linetype = 0) +  
      geom_segment(data = data.frame(x = 24.4, label = "R0 Optimum"),  
                   aes(x = x, xend = x, y = 0, yend = Inf, color = label),  
                   linewidth = 0.7, 
                   linetype = "dashed") +
          scale_color_manual(
            values = c("1Data" = "black","2GAM Model" = "red", 
                       "R0 Optimum" = "#56B4E9"),  
            labels = c("Data mean" ,"GAM fit", 
                       "GAM Optimum", expression(R[0]^{"rel"} * " optimum"))
      )
  }else{
    plot_gam <- plot_gam +
      # Theoretical R0 optimum
      geom_rect(data = data.frame(xmin = 23.0, xmax = 25.8, ymin = 0, ymax = Inf, 
                                  label = "R0 Optimum"),  
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                    color = label),  
                fill = "#56B4E9",
                alpha = 0.3,
                linetype = 0) +  
      geom_segment(data = data.frame(x = 24.4, label = "R0 Optimum"),  
                   aes(x = x, xend = x, y = 0, yend = Inf, color = label),  
                   linewidth = 0.7, 
                   linetype = "dashed") +  
      scale_color_manual(
        values = c("1Data" = "black","2GAM Model" = "red", 
                   "R0 Optimum" = "#56B4E9"),  
        labels = c("Data mean ± SE" ,"GAM fit", 
                   "GAM Optimum", expression(R[0]^{"rel"} * " optimum"))
      )
  }
  if(binary == T){
    plot_data <- NULL
  }else{
    plot_data <- ggplot() +
      geom_point(
        data = df,
        aes(x = .data[[predictor]], y = .data[[y]]),
        shape = 19,
        size = 1,
        alpha = 0.2
      ) +
      labs(  
        y = y_label
      ) +
      theme_bw(base_size=12)+
      theme(
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.75),
        panel.grid = element_blank()
      ) 
    
    if(r0_predictor == F){
      plot_data <- plot_data +  
        geom_rect(data = data.frame(xmin = 23.0, xmax = 25.8, ymin = 0, ymax = Inf, 
                                    label = "R0 Optimum"),  
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                      color = label),  
                  fill = "#56B4E9",
                  alpha = 0.3,
                  linetype = 0) +  
        geom_segment(data = data.frame(x = 24.4, label = "R0 Optimum"),  
                     aes(x = x, xend = x, y = 0, yend = Inf, color = label),  
                     linewidth = 0.7, 
                     linetype = "dashed") +
        scale_color_manual(
          values = c("R0 Optimum" = "#56B4E9"),  
          labels = c(expression(R[0]^{"rel"} * " optimum"))
        ) 
    }
  }
  
  return(list(model_plot = plot_gam, data_plot = plot_data, 
              optimum = optimal_pred, CI_optimum = CI_optimal))
}