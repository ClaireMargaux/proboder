#' Plot simulated contact rate
#'
#' This function generates a plot of the simulated contact rate
#'
#' @param real_beta_df Data frame containing the real contact rate values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate.
#' @return A ggplot object displaying the contact rate plot.
#' @export
plot_sim_contact_rate <- function(real_beta_df, latency_rate, recovery_rate, fatality_rate) {
  ggplot() +
    geom_line(data = real_beta_df, aes(x = t, y = beta, color = "Simulated contact rate"), linetype = "dashed") +
    labs(x = "Time", y = "Simulated contact rate", title = "Simulated contact rate",
         color = "Legend") +  
    scale_color_manual(values = c("Simulated contact rate" = "#E69F00"),
                       labels = c("Simulated contact rate"), name = "Line") + 
    theme_minimal() +
    ylim(c(0,1)) +
    theme(legend.position = "top") +
    guides(color = guide_legend(order = 1))
}

#' Plot contact rate inferred from simulated data
#'
#' This function generates a plot of the estimated contact rate along with its 95%-confidence interval.
#'
#' @param U_plot Data frame containing the estimated contact rate values.
#' @param real_beta_df Data frame containing the real contact rate values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate.
#' @param lengthscale Numeric value indicating the length scale.
#' @return A ggplot object displaying the contact rate plot.
#' @export
plot_contact_rate_with_CI <- function(U_plot, real_beta_df = NULL, latency_rate, recovery_rate, fatality_rate, lengthscale) {
    ggplot() +
      geom_ribbon(data = U_plot, aes(x = t, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = real_beta_df, aes(x = t, y = beta, color = "Simulated contact rate"), linetype = "dashed") +
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
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Latency Rate:", latency_rate), hjust = 0, vjust = -3, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Recovery Rate:", recovery_rate), hjust = 0, vjust = -1, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Fatality Rate:", fatality_rate), hjust = 0, vjust = 1, size = 3) +
      annotate("text", x = min(U_plot$t), y = max(U_plot$U_scaled), label = paste("Length Scale:", lengthscale), hjust = 0, vjust = 3, size = 3)
}

#' Plot inferred compartment counts.
#'
#' This function generates a plot of the observed and inferred compartment counts.
#'
#' @param obs Data frame containing the simulated compartment counts.
#' @param obs_with_noise Data frame containing the noisy compartment counts.
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the simulated, noisy and inferred compartment counts.
#' @export
plot_data_sim <- function(obs, obs_with_noise, X_plot) {
  plot_title <- "Compartments counts"
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442",
            "D" = "#0072B2")
  
  ggplot() +
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
    labs(x = "Time", y = "Count (log-scale)", title = plot_title) +
    theme_minimal() +
    scale_y_continuous(trans = 'log10') + 
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
}

#' Plot inferred compartment counts separately.
#'
#' This function generates plots of the different SEIRD compartments containing
#' each simulated data, noisy data, and inferred data.
#'
#' @param obs Data frame containing the observed compartment counts.
#' @param obs_with_noise Data frame containing the noisy compartment counts.
#' @param X_plot Data frame containing the inferred compartment counts.
#' @return A list containing five ggplot objects comparing the SEIRD compartments.
#' @export
plot_compartment <- function(obs,obs_with_noise,X_plot) {
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442",
            "D" = "#0072B2")
  
  plot_S <- ggplot() +
    geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
    geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
    geom_line(data = X_plot, aes(x = t, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minS, ymax = maxS), fill = "#E69F00", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of S") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_E <- ggplot() +
    geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
    geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
    geom_line(data = X_plot, aes(x = t, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minE, ymax = maxE), fill = "#56B4E9", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of E") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_I <- ggplot() +
    geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
    geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
    geom_line(data = X_plot, aes(x = t, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
    geom_ribbon(data = X_plot, aes(x = t, ymin = minI, ymax = maxI), fill = "#009E73", alpha = 0.3) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Detail of I") +
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_R <- ggplot() +
    geom_line(data = X_plot, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
    geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
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
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
  
  plot_D <- ggplot() +
    geom_line(data = X_plot, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
    geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
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
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
  
  return(list(plot_S, plot_E, plot_I, plot_R, plot_D))
}
