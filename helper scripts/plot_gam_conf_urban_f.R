# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
# ==============================================================================

library(mgcv)
library(tidyverse)
library(ggplot2)

#-------------------------------------------------------------------------------
# Define GAM plotting function 
#-------------------------------------------------------------------------------

plot_gam_conf_urban_f <- function(df, model, predictor, r0_predictor = F, type = "monthly",
                       binary = F, n_bins = 30){
  
  if(r0_predictor){
    new_data <- data.frame(
      seq(0,
          #min(df[[predictor]], na.rm = TRUE),
          max(df[[predictor]], na.rm = TRUE),
          length.out = 1000),
      median(df$urban_fabric)
    )
  }else{
    new_data <- data.frame(
      seq(10,
          #min(df[[predictor]], na.rm = TRUE),
          max(df[[predictor]], na.rm = TRUE),
          length.out = 1000),
      median(df$urban_fabric)
    )
  }
  
  names(new_data) <- c(predictor, "urban_fabric")
  # predict with GAM on new data
  prediction <- predict(
    model, se.fit = TRUE, newdata = new_data,
    unconditional = TRUE, type = "link",
    exclude = c("s(urban_fabric)")
  )
  
  # extract mean prediction and 95CI of that mean prediction
  new_data$mean_prediction <- prediction$fit
  new_data$prediction_upper_95CI <- new_data$mean_prediction + 1.96*prediction$se.fit
  new_data$prediction_lower_95CI <- new_data$mean_prediction - 1.96*prediction$se.fit 
  
  # apply inverse of link function to get to response scale
  new_data$mean_prediction  <- model$family$linkinv(new_data$mean_prediction)
  new_data$prediction_upper_95CI <- model$family$linkinv(new_data$prediction_upper_95CI)
  new_data$prediction_lower_95CI <- model$family$linkinv(new_data$prediction_lower_95CI)
  
  #extract optimum predictor value
  optimal_pred <- new_data[[predictor]][which.max(new_data$mean_prediction)]
  
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
      legend.position.inside = c(0.25, 0.75),
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
        values = c("1Data" = "black","2GAM Model" = "red", "R0 Optimum" = "#56B4E9"),  
        labels = c("Data" ,"GAM fit", expression(R[0]^{"rel"} * " optimum"))
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
        values = c("1Data" = "black","2GAM Model" = "red", "R0 Optimum" = "#56B4E9"),  
        labels = c("Data mean ± SE" ,"GAM fit", expression(R[0]^{"rel"} * " optimum"))
      )
  }
  if(binary == T){
    plot_data <- NULL
  }else{
    plot_data <- ggplot() +
      geom_point(
        data = df,
        aes(x = .data[[predictor]], y = .data[[y]], color='Data'),
        shape = 19,
        size = 1,
        alpha = 0.2
      ) +
      labs(  
        y = y_label
      ) +
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
  }
  
  return(list(model_plot = plot_gam, data_plot = plot_data, 
              optimum= optimal_pred))
}