#' Plot inferred transmission rate.
#'
#' This function generates a plot of the estimated transmission rate along with its 95%-confidence interval.
#'
#' @param U_plot Data frame containing the estimated transmission rate values.
#' @param df_beta Data frame containing the real transmission rate values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate (optional).
#' @param lengthscale Numeric value indicating the length scale.
#' @return A ggplot object displaying the transmission rate plot.
#' @export
plot_transmission_rate_with_CI <- function(U_plot, 
                                      df_beta = NULL, 
                                      latency_rate, recovery_rate, 
                                      fatality_rate = NULL, 
                                      lengthscale) {
    
  if (!is.null(df_beta)) { # df_beta exists
    
    ggplot() +
      geom_ribbon(data = U_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = df_beta, aes(x = t, y = beta, color = "Simulated transmission rate"), linetype = "dashed") +
      geom_line(data = U_plot, aes(x = t, y = U_scaled, color = "Inferred transmission rate"), linewidth = 1) +
      labs(x = "Time", 
           y = TeX("Transmission rate $\\beta(t)$"), 
           title = TeX("Transmission rate $\\beta(t)$ with 95\\%-confidence interval"),
           color = "Legend") +  
      scale_color_manual(values = c("Inferred transmission rate" = "#009E73", "Simulated transmission rate" = "#E69F00"),
                         labels = c("Inferred transmission rate", "Simulated transmission rate"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#66FFCC"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Length scale:", lengthscale), hjust = 0, vjust = -10, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Latency rate:", round(latency_rate,4)), hjust = 0, vjust = -8, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Recovery rate:", round(recovery_rate,4)), hjust = 0, vjust = -6, size = 3) +
      if (fatality_rate > 0) annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Fatality rate:", round(fatality_rate,4)), hjust = 0, vjust = -4, size = 3) else NULL 

  } else { # df_beta does not exist
    
    ggplot() +
      geom_ribbon(data = U_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = U_plot, aes(x = t, y = U_scaled, color = "Inferred transmission rate"), linewidth = 1) +
      labs(x = "Time", 
           y = TeX("Transmission rate $\\beta(t)$"), 
           title = TeX("Transmission rate $\\beta(t)$ with 95\\%-confidence interval"),
           color = "Legend") +  
      scale_color_manual(values = c("Inferred transmission rate" = "#009E73"),
                         labels = c("Inferred transmission rate"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#66FFCC"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Length scale:", lengthscale), hjust = 0, vjust = -10, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Latency rate:", round(latency_rate,4)), hjust = 0, vjust = -8, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Recovery rate:", round(recovery_rate,4)), hjust = 0, vjust = -6, size = 3) +
      if (fatality_rate > 0) annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Fatality rate:", round(fatality_rate,4)), hjust = 0, vjust = -4, size = 3) else NULL 
  }
}

#' Plot inferred reproduction number.
#'
#' This function generates a plot of the estimated reproduction number along with its 95%-confidence interval.
#'
#' @param R_plot Data frame containing the estimated reproduction number values.
#' @param df_R Data frame containing the real reproduction number values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate (optional).
#' @param lengthscale Numeric value indicating the length scale.
#' @return A ggplot object displaying the transmission rate plot.
#' @export
plot_reproduction_number_with_CI <- function(R_plot, 
                                           df_R = NULL, 
                                           latency_rate, recovery_rate, 
                                           fatality_rate = NULL, 
                                           lengthscale) {
  
  if (!is.null(df_R)) { # df_R exists
    
    ggplot() +
      geom_ribbon(data = R_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = df_R, aes(x = t, y = R, color = "Simulated reproduction number"), linetype = "dashed") +
      geom_line(data = R_plot, aes(x = t, y = R, color = "Inferred reproduction number"), linewidth = 1) +
      labs(x = "Time", 
           y = TeX("Reproduction number $R(t)$"), 
           title = TeX("Reproduction number $R(t)$ with 95\\%-confidence interval"),
           color = "Legend") +  
      scale_color_manual(values = c("Inferred reproduction number" = "#0072B2", "Simulated reproduction number" = "#E69F00"),
                         labels = c("Inferred reproduction number", "Simulated reproduction number"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#56B4E9"),
                        labels = "Approx. 95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Length scale:", lengthscale), hjust = 0, vjust = -10, size = 3) +
      annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Latency rate:", round(latency_rate,4)), hjust = 0, vjust = -8, size = 3) +
      annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Recovery rate:", round(recovery_rate,4)), hjust = 0, vjust = -6, size = 3) +
      if (fatality_rate > 0) annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Fatality rate:", round(fatality_rate,4)), hjust = 0, vjust = -4, size = 3) else NULL 
    
  } else { # df_R does not exist
    
    ggplot() +
      geom_ribbon(data = R_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = R_plot, aes(x = t, y = R, color = "Inferred reproduction number"), linewidth = 1) +
      labs(x = "Time", 
           y = TeX("Reproduction number $R(t)$"), 
           title = TeX("Reproduction number $R(t)$ with 95\\%-confidence interval"),
           color = "Legend") +  
      scale_color_manual(values = c("Inferred reproduction number" = "#0072B2"),
                         labels = c("Inferred reproduction number"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#56B4E9"),
                        labels = "Approx.  95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Length scale:", lengthscale), hjust = 0, vjust = -10, size = 3) +
      annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Latency rate:", round(latency_rate,4)), hjust = 0, vjust = -8, size = 3) +
      annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Recovery rate:", round(recovery_rate,4)), hjust = 0, vjust = -6, size = 3) +
      if (fatality_rate > 0) annotate("text", x = min(R_plot$t), y = max(R_plot$R), label = paste("Fatality rate:", round(fatality_rate,4)), hjust = 0, vjust = -4, size = 3) else NULL 
    
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
  
  cols_SEIRD <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442", D = "#0072B2")
  labels_SEIRD <- c(S = "Susceptible", E = "Exposed", I = "Infected", R = "Recovered", D = "Deceased")
  
  cols_SEIR <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442")
  labels_SEIR <- c(S = "Susceptible", E = "Exposed", I = "Infected", R = "Recovered")

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

#' Plot inferred compartment counts from simulated data.
#'
#' This function generates a plot of the observed and inferred compartment counts,
#' excluding the S compartment and without log scaling.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param obs Data frame containing the simulated compartment counts.
#' @param obs_with_noise Data frame containing the noisy compartment counts (optional).
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the simulated, noisy, and inferred compartment counts.
#' @export
plot_compartments_except_S <- function(model,
                              obs, 
                              obs_with_noise = NULL, 
                              X_plot) {
  
  plot_title <- "Compartments counts (excluding S)"
  
  cols_SEIRD <- c(E = "#56B4E9", I = "#009E73", R = "#F0E442", D = "#0072B2")
  labels_SEIRD <- c(E = "Exposed", I = "Infected", R = "Recovered", D = "Deceased")
  
  cols_SEIR <- c(E = "#56B4E9", I = "#009E73", R = "#F0E442")
  labels_SEIR <- c(E = "Exposed", I = "Infected", R = "Recovered")
  
  if (model == 'SEIRD') {
    
    if (!is.null(obs_with_noise)) { # obs_with_noise exists
      
      p <- ggplot() +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "Compartment counts", title = plot_title) +
        theme_minimal() +
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
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "Compartment counts", title = plot_title) +
        theme_minimal() +
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
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "Compartment counts", title = plot_title) +
        theme_minimal() +
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
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
        geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
        geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
        labs(x = "Time", y = "Compartment counts", title = plot_title) +
        theme_minimal() +
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
  
  cols_SEIRD <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442", D = "#0072B2")

  cols_SEIR <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442")
  
  if (!is.null(obs_with_noise)) { # obs_with_noise available
    
    plot_S <- ggplot() +
      geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
      geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
      geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
      geom_ribbon(data = X_plot, aes(x = t, ymin = minS, ymax = maxS), fill = "#E69F00", alpha = 0.3) +
      theme_minimal() +
      ylab("Counts") +
      xlab("Time") + 
      ggtitle("Detail of Susceptible") +
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
      ggtitle("Detail of Exposed") +
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
      ggtitle("Detail of Infected") +
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
      ggtitle("Detail of Recovered") +
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
        ggtitle("Detail of Deceased") +
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
      ggtitle("Detail of Susceptible") +
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
      ggtitle("Detail of Exposed") +
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
      ggtitle("Detail of Infected") +
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
      ggtitle("Detail of Recovered") +
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
        ggtitle("Detail of Deceased") +
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

#' Plot Observed data
#' 
#' Plot the Observed compartment counts.
#'
#' @param obs Data frame, observed data.
#' @param log Logical indicating whether to plot the counts on a log scale.
#' 
#' @return ggplot object: the counts plot.
#' @export
plotting_real_data <- function(obs, log = TRUE) {

  cols_SEIR <- c(R = "#F0E442")
  labels_SEIR <- c(R = "Recovered")
  

  obs <- obs %>%
    mutate(
      R = ifelse(R <= 0, 0.001, R)
      )
    
  # Create the counts plot
  real_compartments <- ggplot(obs, aes(x = t)) +
    geom_line(aes(y = R, color = 'R')) +
    theme_minimal() +
    ylab("Compartment counts") +
    xlab("Time") + 
    ggtitle("Real observations") +
    scale_color_manual(values = cols_SEIR,
                       labels = labels_SEIR, 
                       name = "Compartments") +
    guides(color = guide_legend(order = 1)) +
    if (log) scale_y_continuous(trans = 'log10', name = "log(compartment counts)") else NULL
  
  # Return the plot
  return(real_compartments)
}

#' Plot inferred compartment counts.
#'
#' This function generates a plot of the observed and inferred compartment counts.
#'
#' @param obs Data frame containing the observed compartment counts.
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the observed and inferred compartment counts.
#' @export
plot_compartments_real <- function(obs, X_plot) {
  plot_title <- "Compartments counts"
  
  cols_SEIR <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442")
  labels_SEIR <- c(S = "Susceptible", E = "Exposed", I = "Infected", R = "Recovered")
  
  obs <- obs %>%
    mutate(
      R = ifelse(R <= 0, 0.001, R)
    )

  X_plot <- X_plot %>%
    mutate(
      S = ifelse(S <= 0, 0.001, S),
      E = ifelse(E <= 0, 0.001, E),
      I = ifelse(I <= 0, 0.001, I),
      R = ifelse(R <= 0, 0.001, R)
    )
  
  ggplot() +
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
plot_compartments_separately_real <- function(obs,X_plot) {
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442")
  
  plot_S <- ggplot() +
    geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minS, ymax = maxS), fill = "#E69F00", alpha = 0.3) +
    theme_minimal() +
    ylab("Counts") +
    xlab("Time") +
    ggtitle("Detail of Susceptible") +
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
    geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minE, ymax = maxE), fill = "#56B4E9", alpha = 0.3) +
    theme_minimal() +
    ylab("Counts") +
    xlab("Time") +
    ggtitle("Detail of Exposed") +
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
    geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minI, ymax = maxI), fill = "#009E73", alpha = 0.3) +
    theme_minimal() +
    ylab("Counts") +
    xlab("Time") + 
    ggtitle("Detail of Infected") +
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
    ylab("Counts") +
    xlab("Time") + 
    ggtitle("Detail of Recovered") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Observed data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) 
  
  
  return(list(plot_S, plot_E, plot_I, plot_R))
}

################################################################################
################################################################################

#' Plot condition numbers over time
#'
#' This function generates and displays four separate line plots for the condition numbers 
#' of four matrices (S_obs, S_ode, P_X, P_U) over time. The plots are arranged in a 2x2 grid.
#'
#' @param data A data frame containing the condition numbers of four matrices over time.
#' The data frame should have the following columns:
#' - `time_grid`: The time points.
#' - `cond_S_obs`: Condition number of the S_obs matrix.
#' - `cond_S_ode`: Condition number of the S_ode matrix.
#' - `cond_P_X`: Condition number of the P_X matrix.
#' - `cond_P_U`: Condition number of the P_U matrix.
#'
#' @return A grid of four plots, each representing the condition number of a different matrix over time.
#' The plots are displayed in a 2x2 arrangement.
#'
#' @import ggplot2
#' @import gridExtra
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with a data frame `cond_plot`:
#' plot_condition_numbers(cond_plot)
#' }

plot_condition_numbers <- function(data) {
  
  # Plot for cond_S_obs
  p1 <- ggplot(data, aes(x = time_grid, y = cond_S_obs)) +
    geom_line(color = "#E69F00") +
    labs(title = TeX("Condition number: $S_{obs}$"),
         x = TeX("Time"),
         y = TeX("Condition number ($\\kappa$)")) +
    theme_minimal()
  
  # Plot for cond_S_ode
  p2 <- ggplot(data, aes(x = time_grid, y = cond_S_ode)) +
    geom_line(color = "#56B4E9") +
    labs(title = TeX("Condition number: $S_{ODE}$"),
         x = TeX("Time"),
         y = TeX("Condition number ($\\kappa$)")) +
    theme_minimal()
  
  # Plot for cond_P_X
  p3 <- ggplot(data, aes(x = time_grid, y = cond_P_X)) +
    geom_line(color = "#009E73") +
    labs(title = TeX("Condition number: $P_X$"),
         x = TeX("Time"),
         y = TeX("Condition number ($\\kappa$)")) +
    theme_minimal()
  
  # Plot for cond_P_U
  p4 <- ggplot(data, aes(x = time_grid, y = cond_P_U)) +
    geom_line(color = "#F0E442") +
    labs(title = TeX("Condition number: $P_U$"),
         x = TeX("Time"),
         y = TeX("Condition number ($\\kappa$)")) +
    theme_minimal()
  
  res <- list(p1, p2, p3, p4)
  
  return(res)
}
