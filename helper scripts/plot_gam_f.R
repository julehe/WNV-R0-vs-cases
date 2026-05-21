# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
#
# Description: This provides a helper function for plotting of GAM results
#
# ==============================================================================

library(mgcv)
library(tidyverse)
library(ggplot2)
library(MASS)
library(rlang)

#-------------------------------------------------------------------------------
# Define GAM plotting function 
#-------------------------------------------------------------------------------

plot_gam_f <- function(df, model, predictor, adjustments = c(NA), 
                       r0_predictor = F, type = "monthly", binary = F, 
                       n_bins = 15, country_effects = T){
  
  if(r0_predictor){
    lower_bound = 0
  }else{
    lower_bound = 10
  }
  if(country_effects == T){
    country_max <- df %>%
      group_by(country) %>%
      summarise(
        max_x = max(.data[[predictor]], na.rm = TRUE),
        .groups = "drop"
      )
    
    new_data <- country_max %>%
      group_by(country) %>%
      reframe(
        seq(lower_bound, 
            #max(df[[predictor]], na.rm = TRUE),
            max_x, 
            length.out = 1000)
      )
    
    names(new_data)[2] <- predictor
    
    if(length(adjustments) > 0 && all(!is.na(adjustments))){
      
      country_medians <- df %>%
        group_by(country) %>%
        summarise(
          across(all_of(adjustments), \(x) median(x, na.rm = TRUE)),
          .groups = "drop"
        )
      
      new_data <- dplyr::left_join(new_data, country_medians, by = "country")
    }
  }else{
    new_data <- data.frame(
      seq(lower_bound,
          max(df[[predictor]], na.rm = TRUE),
          length.out = 1000)
    )
    
    names(new_data) <- predictor 
    
    if(length(adjustments) > 0 && all(!is.na(adjustments))){
      new_data[adjustments] <- lapply(model$model[adjustments], median)
    }
  }
  
  # predict with GAM on new data
  prediction <- predict(
    model, se.fit = TRUE, newdata = new_data,
    unconditional = TRUE, type = "link"
  )
  
  exclude_term = paste0("s(country,",predictor,")")
  # shared smooth only
  prediction_shared <- predict(
    model, se.fit = TRUE, newdata = new_data,
    unconditional = TRUE, type = "link",
    exclude = c(exclude_term,
                "urban_fabric",
                "rain_3mo_cum_pop_weight",
                "percent_65plus",
                "arable_land",
                "gdp_per_capita")
  )
  
  # extract mean prediction and 95CI of that mean prediction
  new_data$mean_prediction <- prediction$fit
  new_data$prediction_upper_95CI <- new_data$mean_prediction + 1.96*prediction$se.fit
  new_data$prediction_lower_95CI <- new_data$mean_prediction - 1.96*prediction$se.fit 
  
  new_data$mean_prediction_shared <- prediction_shared$fit
  new_data$prediction_upper_95CI_shared <- new_data$mean_prediction_shared + 1.96*prediction_shared$se.fit
  new_data$prediction_lower_95CI_shared <- new_data$mean_prediction_shared - 1.96*prediction_shared$se.fit 
  
  # apply inverse of link function to get to response scale
  new_data$mean_prediction  <- model$family$linkinv(new_data$mean_prediction)
  new_data$prediction_upper_95CI <- model$family$linkinv(new_data$prediction_upper_95CI)
  new_data$prediction_lower_95CI <- model$family$linkinv(new_data$prediction_lower_95CI)
  
  new_data$mean_prediction_shared  <- model$family$linkinv(new_data$mean_prediction_shared)
  new_data$prediction_upper_95CI_shared <- model$family$linkinv(new_data$prediction_upper_95CI_shared)
  new_data$prediction_lower_95CI_shared <- model$family$linkinv(new_data$prediction_lower_95CI_shared)
  
  if(binary == F){
    # scale predictions to incidence per 100K
    new_data$mean_prediction <- new_data$mean_prediction*100000
    new_data$prediction_upper_95CI <- new_data$prediction_upper_95CI*100000
    new_data$prediction_lower_95CI <- new_data$prediction_lower_95CI*100000
    
    new_data$mean_prediction_shared <- new_data$mean_prediction_shared*100000
    new_data$prediction_upper_95CI_shared <- new_data$prediction_upper_95CI_shared*100000
    new_data$prediction_lower_95CI_shared <- new_data$prediction_lower_95CI_shared*100000
  }
  
  country_labels = c(
    AT = "Austria",
    BG = "Bulgaria",
    CY = "Cyprus",
    CZ = "Czech Republic",
    DE = "Germany",
    EL = "Greece",
    ES = "Spain",
    FR = "France",
    HR = "Croatia",
    HU = "Hungary",
    IT = "Italy",
    NL = "Netherlands",
    PT = "Portugal",
    RO = "Romania",
    SI = "Slovenia",
    SK = "Slovakia"
  )
  
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
        size = 0.5,
        alpha = 0.1
      ) 
  }else{
    plot_gam <- ggplot() +
      # binned data
      stat_summary_bin(data = df, 
                       aes(x = .data[[predictor]], y = .data[[y]], color='1Data'), 
                       fun.data = function(y) {
                         mean_y <- mean(y)
                         se_y <- sd(y)/sqrt(length(y))
                         return(data.frame(y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y))
                       }, 
                       bins = n_bins, 
                       linewidth = 0.2,
                       geom = 'errorbar') +
      stat_summary_bin(data = df,
                       aes(x = .data[[predictor]], y = .data[[y]], color='1Data'),
                       fun = 'mean', bins = n_bins,
                       size = 0.5, geom='point')
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
      aes(x = .data[[predictor]], y = mean_prediction, color = "2GAM Model",
          group = if(country_effects == T){
            country
          }else{
            1
          }),
      linewidth = 0.4,
      linetype = "solid"
    ) +
    labs(  
      y = y_label
    ) +
    theme_bw(base_size = 7)+
    guides() + 
    theme(
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.3, 0.75),
      legend.key.size = unit(0.75,"line"),
      panel.grid = element_blank()
    )
  
  if(r0_predictor == T){
    plot_gam <- plot_gam +  
        geom_rug(
          data = df,
          aes(x = .data[[predictor]]),
          inherit.aes = FALSE,
          alpha = 0.1,
          sides = "b"
        ) +
      scale_color_manual(
        values = c("1Data" = "black","2GAM Model" = "red"),  
        labels = c("Data mean ± SE" ,"GAM fit")
      )
  }else if(binary == T){
    plot_gam <- plot_gam +
      geom_rug(
        data = df,
        aes(x = .data[[predictor]]),
        inherit.aes = FALSE,
        alpha = 0.1,
        sides = "b"
      ) +
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
                   linewidth = 0.4, 
                   linetype = "dashed") +
          scale_color_manual(
            values = c("1Data" = "black","2GAM Model" = "red", 
                       "R0 Optimum" = "#56B4E9"),  
            labels = c("Data" ,"GAM fit", 
                       expression(R[0]^{"rel"} * " optimum"))
      )
  }else{
    plot_gam <- plot_gam +
      geom_rug(
        data = df,
        aes(x = .data[[predictor]]),
        inherit.aes = FALSE,
        alpha = 0.1,
        sides = "b"
      ) +
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
                   linewidth = 0.4, 
                   linetype = "dashed") +  
      scale_color_manual(
        values = c("1Data" = "black","2GAM Model" = "red", 
                   "R0 Optimum" = "#56B4E9"),  
        labels = c("Data mean ± SE" ,"GAM fit",
                   expression(R[0]^{"rel"} * " optimum"))
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
        size = 0.5,
        alpha = 0.2
      ) +
      labs(  
        y = y_label
      ) +
      theme_bw(base_size = 7)+
      theme(
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.75),
        legend.key.size = unit(0.75,"line"),
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
                     linewidth = 0.4, 
                     linetype = "dashed") +
        scale_color_manual(
          values = c("R0 Optimum" = "#56B4E9"),  
          labels = c(expression(R[0]^{"rel"} * " optimum"))
        ) 
    }
  }
  
  if (country_effects == T){
    plot_gam <- plot_gam + facet_wrap(~country, #scales = "free_y", 
                                      ncol = 4, 
                                      labeller = labeller(country = country_labels)) +
      theme(legend.position = "bottom",
            strip.background = element_blank())
    
    plot_data_shared <- plot_data
    
    plot_data <- plot_data + facet_wrap(~country, #scales = "free_y", 
                                        ncol = 4, 
                                        labeller = labeller(country = country_labels)) +
      theme(legend.position = "bottom",
            strip.background = element_blank())
    
    plot_gam_shared <- ggplot() +
      geom_ribbon(
        data = new_data,
        aes(x = .data[[predictor]], ymin = prediction_lower_95CI_shared, 
            ymax = prediction_upper_95CI_shared, color = "2GAM Model"),
        fill = "red",
        alpha = 0.3,
        linetype = 0
      ) +
      geom_line(
        data = new_data,
        aes(x = .data[[predictor]], y = mean_prediction_shared, 
            group = ,
            group = if(country_effects == T){
              country
            }else{
              1
            },
            color = "2GAM Model"),
        linewidth = 0.4,
        linetype = "solid"
      ) +
      labs(  
        y = y_label
      ) +
      theme_bw(base_size = 7)+
      guides() + 
      theme(
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.75),
        legend.key.size = unit(0.75,"line"),
        panel.grid = element_blank()
      )
    
    if(r0_predictor == T){
      plot_gam_shared <- plot_gam_shared +  
        geom_rug(
          data = df,
          aes(x = .data[[predictor]]),
          inherit.aes = FALSE,
          alpha = 0.03,
          sides = "b"
        ) +
        scale_color_manual(
          values = c("2GAM Model" = "red"),  
          labels = c("GAM fit")
        )
    }else{
      plot_gam_shared <- plot_gam_shared +
        geom_rug(
          data = df,
          aes(x = .data[[predictor]]),
          inherit.aes = FALSE,
          alpha = 0.03,
          sides = "b"
        ) +
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
                     linewidth = 0.4, 
                     linetype = "dashed") +  
        scale_color_manual(
          values = c("2GAM Model" = "red", 
                     "R0 Optimum" = "#56B4E9"),  
          labels = c("GAM fit",
                     expression(R[0]^{"rel"} * " optimum"))
        )
    }
  }
  
  return(list(shared_model_plot = plot_gam_shared, model_plot = plot_gam, 
              shared_data_plot = plot_data_shared,
              data_plot = plot_data))
}