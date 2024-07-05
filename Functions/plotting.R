#' Plot inferred contact rate.
#'
#' This function generates a plot of the estimated contact rate along with its 95%-confidence interval.
#'
#' @param U_plot Data frame containing the estimated contact rate values.
#' @param real_beta_df Data frame containing the real contact rate values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate (optional).
#' @param lengthscale Numeric value indicating the length scale.
#' @return A ggplot object displaying the contact rate plot.
#' @export
plot_contact_rate_with_CI <- function(U_plot, 
                                      df_beta = NULL, 
                                      latency_rate, recovery_rate, 
                                      fatality_rate = NULL, 
                                      lengthscale) {
    
  if (!is.null(df_beta)) { # df_beta exists
    
    ggplot() +
      geom_ribbon(data = U_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = df_beta, aes(x = t, y = beta, color = "Simulated contact rate"), linetype = "dashed") +
      geom_line(data = U_plot, aes(x = t, y = U_scaled, color = "Inferred contact rate"), linewidth = 1) +
      labs(x = "Time", y = "Contact rate", title = "Contact rate with 95%-confidence interval",
           color = "Legend") +  
      scale_color_manual(values = c("Inferred contact rate" = "#009E73", "Simulated contact rate" = "#E69F00"),
                         labels = c("Inferred contact rate", "Simulated contact rate"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#66FFCC"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Length scale:", lengthscale), hjust = 0, vjust = 2, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Latency rate:", round(latency_rate,4)), hjust = 0, vjust = 4, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Recovery rate:", round(recovery_rate,4)), hjust = 0, vjust = 6, size = 3) +
      if (fatality_rate > 0) annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Fatality rate:", round(fatality_rate,4)), hjust = 0, vjust = 8, size = 3) else NULL 

  } else { # df_beta does not exist
    
    ggplot() +
      geom_ribbon(data = U_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = U_plot, aes(x = t, y = U_scaled, color = "Inferred contact rate"), linewidth = 1) +
      labs(x = "Time", y = "Contact rate", title = "Contact rate with 95%-confidence interval",
           color = "Legend") +  
      scale_color_manual(values = c("Inferred contact rate" = "#009E73"),
                         labels = c("Inferred contact rate"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#66FFCC"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Length scale:", lengthscale), hjust = 0, vjust = 2, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Latency rate:", round(latency_rate,4)), hjust = 0, vjust = 4, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Recovery rate:", round(recovery_rate,4)), hjust = 0, vjust = 6, size = 3) +
      if (fatality_rate > 0) annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Fatality rate:", round(fatality_rate,4)), hjust = 0, vjust = 8, size = 3) else NULL 
  }
}

#' Plot inferred compartment counts from simulated data.
#'
#' This function generates a plot of the observed and inferred compartment counts.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param obs Data frame containing the simulated compartment counts.
#' @param obs_with_noise Data frame containing the noisy compartment counts (optional).
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the simulated, noisy and inferred compartment counts.
#' @export
plot_compartments <- function(model,
                              obs, 
                              obs_with_noise = NULL, 
                              X_plot) {
  
  plot_title <- "Compartments counts"
  
  cols_SEIRD <- c("S" = "#E69F00", 
                  "E" = "#56B4E9",
                  "I" = "#009E73", 
                  "R" = "#F0E442",
                  "D" = "#0072B2")
  
  cols_SEIR <- c("S" = "#E69F00", 
                 "E" = "#56B4E9",
                 "I" = "#009E73", 
                 "R" = "#F0E442")
  
  labels_SEIRD <- c("S","E","I","R","D")
  
  labels_SEIR <- c("S","E","I","R")

  obs <- if (model == 'SEIRD') {
    obs %>%
      mutate(
        S = ifelse(S <= 0, 0.01, S),
        E = ifelse(E <= 0, 0.01, E),
        I = ifelse(I <= 0, 0.01, I),
        R = ifelse(R <= 0, 0.01, R),
        D = ifelse(D <= 0, 0.01, D)
      )
  } else {
    obs %>%
      mutate(
        S = ifelse(S <= 0, 0.01, S),
        E = ifelse(E <= 0, 0.01, E),
        I = ifelse(I <= 0, 0.01, I),
        R = ifelse(R <= 0, 0.01, R)
      )
  }
  
  if (!is.null(obs_with_noise)) {
    
    obs_with_noise <- if (model == 'SEIRD') {
      obs_with_noise %>%
        mutate(
          S = ifelse(S <= 0, 0.01, S),
          E = ifelse(E <= 0, 0.01, E),
          I = ifelse(I <= 0, 0.01, I),
          R = ifelse(R <= 0, 0.01, R),
          D = ifelse(D <= 0, 0.01, D)
        )
    } else {
      obs_with_noise %>%
        mutate(
          S = ifelse(S <= 0, 0.01, S),
          E = ifelse(E <= 0, 0.01, E),
          I = ifelse(I <= 0, 0.01, I),
          R = ifelse(R <= 0, 0.01, R)
        )
    }
  }
  
  X_plot <- if (model == 'SEIRD') {
    X_plot %>%
      mutate(
        S = ifelse(S <= 0, 0.01, S),
        E = ifelse(E <= 0, 0.01, E),
        I = ifelse(I <= 0, 0.01, I),
        R = ifelse(R <= 0, 0.01, R),
        D = ifelse(D <= 0, 0.01, D)
      )
  } else {
    X_plot %>%
      mutate(
        S = ifelse(S <= 0, 0.01, S),
        E = ifelse(E <= 0, 0.01, E),
        I = ifelse(I <= 0, 0.01, I),
        R = ifelse(R <= 0, 0.01, R)
      )
  }
  
  if (model == 'SEIRD') {
    
    if (!is.null(obs_with_noise)) { # obs_with_noise exists
      
      p <- ggplot() +
        geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
        geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
        geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "log(compartment counts)", title = plot_title) +
        theme_minimal() +
        scale_y_continuous(trans = 'log10') + 
        scale_color_manual(values = cols_SEIRD,
                           labels = labels_SEIRD,
                           name = "Compartments") +
        scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data", "Inferred data"),
                              guide = guide_legend(override.aes = list(color = "grey"))) +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey")))
      
    } else { # obs_with_noise does not exist
      
      p <- ggplot() +
        geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "log(compartment counts)", title = plot_title) +
        theme_minimal() +
        scale_y_continuous(trans = 'log10') + 
        scale_color_manual(values = cols_SEIRD,
                           labels = labels_SEIRD,
                           name = "Compartments") +
        scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data", "Inferred data"),
                              guide = guide_legend(override.aes = list(color = "grey")))
      
    }
    
  } else if (model == 'SEIR') {
    
    if (!is.null(obs_with_noise)) { # obs_with_noise exists
      
      p <- ggplot() +
        geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "log(compartment counts)", title = plot_title) +
        theme_minimal() +
        scale_y_continuous(trans = 'log10') + 
        scale_color_manual(values = cols_SEIR,
                           labels = labels_SEIR,
                           name = "Compartments") +
        scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data", "Inferred data"),
                              guide = guide_legend(override.aes = list(color = "grey"))) +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey")))
      
    } else { # obs_with_noise does not exist
      
      p <- ggplot() +
        geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "log(compartment counts)", title = plot_title) +
        theme_minimal() +
        scale_y_continuous(trans = 'log10') + 
        scale_color_manual(values = cols_SEIR,
                           labels = labels_SEIR,
                           name = "Compartments") +
        scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data", "Inferred data"),
                              guide = guide_legend(override.aes = list(color = "grey")))
      
    }
    
  }

  return(p)  
}

#' Plot inferred compartment counts separately.
#'
#' This function generates plots of the different SEIRD compartments containing
#' each simulated data, noisy data, and inferred data.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param obs Data frame containing the simulated compartment counts.
#' @param obs_with_noise Data frame containing the noisy compartment counts (optional).
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A list containing five ggplot objects comparing the SEIRD compartments.
#' @export
plot_compartments_separately <- function(model, obs, obs_with_noise = NULL, X_plot) {
  
  cols_SEIRD <- c("S" = "#E69F00", 
                  "E" = "#56B4E9",
                  "I" = "#009E73", 
                  "R" = "#F0E442",
                  "D" = "#0072B2")
  
  cols_SEIR <- c("S" = "#E69F00", 
                 "E" = "#56B4E9",
                 "I" = "#009E73", 
                 "R" = "#F0E442")
  
  if (!is.null(obs_with_noise)) { # obs_with_noise available
    
    plot_S <- ggplot() +
      geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
      geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
      geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minS, ymax = maxS), fill = "#E69F00", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of S") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey"))) +
      scale_shape_manual(values = c("Noisy data" = 1),
                         name = "Points",
                         labels = c("Noisy data"),
                         guide = guide_legend(override.aes = list(color = "grey")))
    
    plot_E <- ggplot() +
      geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
      geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
      geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minE, ymax = maxE), fill = "#56B4E9", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of E") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey"))) +
      scale_shape_manual(values = c("Noisy data" = 1),
                         name = "Points",
                         labels = c("Noisy data"),
                         guide = guide_legend(override.aes = list(color = "grey")))
    
    plot_I <- ggplot() +
      geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
      geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
      geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minI, ymax = maxI), fill = "#009E73", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of I") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey"))) +
      scale_shape_manual(values = c("Noisy data" = 1),
                         name = "Points",
                         labels = c("Noisy data"),
                         guide = guide_legend(override.aes = list(color = "grey")))
    
    plot_R <- ggplot() +
      geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
      geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
      geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minR, ymax = maxR), fill = "#F0E442", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of R") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey"))) +
      scale_shape_manual(values = c("Noisy data" = 1),
                         name = "Points",
                         labels = c("Noisy data"),
                         guide = guide_legend(override.aes = list(color = "grey")))
    
    if (model == 'SEIRD') {
      
      plot_D <- ggplot() +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
        geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
        geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
        geom_ribbon(data = X_plot, aes(x = t, ymin = minD, ymax = maxD), fill = "#0072B2", alpha = 0.3) +
        theme_minimal() +
        ylab("Counts") +
        xlab("Time") + 
        ggtitle("Detail of D") +
        scale_color_manual(values = cols_SEIRD,
                           name = "Compartment") +
        scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data", "Inferred data"),
                              guide = guide_legend(override.aes = list(color = "grey"))) +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey")))
      
    }
    
    if (model == 'SEIRD') {
      res <- list(plot_S, plot_E, plot_I, plot_R, plot_D)
    } else {
      res <- list(plot_S, plot_E, plot_I, plot_R)
    }
    
  } else { # no obs_with_noise available
    
    plot_S <- ggplot() +
      geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
      geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minS, ymax = maxS), fill = "#E69F00", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of S") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey"))) 
    
    plot_E <- ggplot() +
      geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
      geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minE, ymax = maxE), fill = "#56B4E9", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of E") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey")))
    
    plot_I <- ggplot() +
      geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
      geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minI, ymax = maxI), fill = "#009E73", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of I") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey")))
    
    plot_R <- ggplot() +
      geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
      geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minR, ymax = maxR), fill = "#F0E442", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of R") +
      scale_color_manual(values = cols_SEIRD,
                         name = "Compartment") +
      scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                            name = "Lines",
                            labels = c("Simulated data", "Inferred data"),
                            guide = guide_legend(override.aes = list(color = "grey")))
    
    if (model == 'SEIRD') {
      
      plot_D <- ggplot() +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
        geom_ribbon(data = X_plot, aes(x = t, ymin = minD, ymax = maxD), fill = "#0072B2", alpha = 0.3) +
        theme_minimal() +
        ylab("Counts") +
        xlab("Time") + 
        ggtitle("Detail of D") +
        scale_color_manual(values = cols_SEIRD,
                           name = "Compartment") +
        scale_linetype_manual(values = c("Data" = "dashed", "Inferred" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data", "Inferred data"),
                              guide = guide_legend(override.aes = list(color = "grey")))
      
    }
    
    if (model == 'SEIRD') {
      res <- list(plot_S, plot_E, plot_I, plot_R, plot_D)
    } else {
      res <- list(plot_S, plot_E, plot_I, plot_R)
    }
    
  }
  
  return(res)
}

################################################################################
################################################################################

#' Plot inferred compartment counts.
#'
#' This function generates a plot of the observed and inferred compartment counts.
#'
#' @param obs Data frame containing the observed compartment counts.
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the observed and inferred compartment counts.
#' @export
plot_data_real <- function(obs, X_plot) {
  plot_title <- "Compartments counts"
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442",
            "D" = "#0072B2")
  
  obs <- obs %>%
    mutate(
      #S = ifelse(S <= 0, 0.001, S),
      #E = ifelse(E <= 0, 0.001, E),
      #I = ifelse(I <= 0, 0.001, I),
      R = ifelse(R <= 0, 0.001, R),
      D = ifelse(D <= 0, 0.001, D)
    )

  X_plot <- X_plot %>%
    mutate(
      S = ifelse(S <= 0, 0.001, S),
      E = ifelse(E <= 0, 0.001, E),
      I = ifelse(I <= 0, 0.001, I),
      R = ifelse(R <= 0, 0.001, R),
      D = ifelse(D <= 0, 0.001, D)
    )
  
  ggplot() +
    #geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
    #geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
    #geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
    labs(x = "Time", y = "Count (log-scale)", title = plot_title) +
    theme_minimal() +
    scale_y_continuous(trans = 'log10') + 
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Observed data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) 
}

#' Plot inferred compartment counts separately.
#'
#' This function generates plots of the different SEIRD compartments containing
#' each simulated data, noisy data, and inferred data.
#'
#' @param obs Data frame containing the observed compartment counts.
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A list containing five ggplot objects comparing the SEIRD compartments.
#' @export
plot_compartment_real <- function(obs,X_plot) {
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442",
            "D" = "#0072B2")
  
  plot_S <- ggplot() +
    #geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minS, ymax = maxS), fill = "#E69F00", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of S") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    # scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
    #                       name = "Lines",
    #                       labels = c("Observed data", "Inferred data"),
    #                       guide = guide_legend(override.aes = list(color = "grey")))
    scale_linetype_manual(values = c("Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_E <- ggplot() +
    #geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minE, ymax = maxE), fill = "#56B4E9", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of E") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    # scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
    #                       name = "Lines",
    #                       labels = c("Observed data", "Inferred data"),
    #                       guide = guide_legend(override.aes = list(color = "grey")))
    scale_linetype_manual(values = c("Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_I <- ggplot() +
    #geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minI, ymax = maxI), fill = "#009E73", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of I") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    # scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
    #                       name = "Lines",
    #                       labels = c("Observed data", "Inferred data"),
    #                       guide = guide_legend(override.aes = list(color = "grey")))
    scale_linetype_manual(values = c("Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_R <- ggplot() +
    geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minR, ymax = maxR), fill = "#F0E442", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of R") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Observed data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) 
  
  plot_D <- ggplot() +
    geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minD, ymax = maxD), fill = "#0072B2", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of D") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Observed data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) 
  
  return(list(plot_S, plot_E, plot_I, plot_R, plot_D))
}
