###########################################
####### SIMULATE DATA USING LSODA #########
###########################################

#' Simulate data using LSODA
#' 
#' Simulate compartment counts (with or without noise) using the function lsoda from the package deSolve.
#' 
#' @param model String, type of model to be used for simulation ('SEIRD' or 'SEIR' available).
#' @param noise Variance of the Gaussian white noise to be added to the data (if set to 0, no noise is added).
#' @param seed Seed to be used for the noise added to the data.
#' @param steps Size of the time steps of the generated data.
#' @param max_time Last time step of the generated data.
#' @param pop Population size of generated data.
#' @param beta Function; transmission rate to be used for simulation.
#' @examples
#' # Beta function: sine
#' # beta <- function(t){0.05*sin(t/10)+0.8}
#' # Beta function: log
#' # beta <- function(t){0.03*log(t+1)+0.8}
#' @param xstart Starting values for the compartments.
#' @param lambda Latency rate.
#' @param gamma Recovery rate.
#' @param eta Fatality rate.
#' 
#' @return List containing three data frames: 
#'   (1) (noisy) simulated counts
#'   (2) predefined transmission rate
#'   (3) corresponding reproduction number
#' @export
simulate_data_LSODA <- function(model, noise = 0, seed = 5, steps = 1, max_time = 30, 
                                lambda, gamma, eta = 0, pop, beta, xstart) {
  
  # Set seed for reproducibility
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # Create time grid
  grid <- seq(0, max_time, by = steps)
  n <- length(grid)
  
  # Define parameters
  if (model == 'SEIRD') {
    parms <- c(lambda = lambda, gamma = gamma, eta = eta, pop = pop)
  }
  if (model == 'SEIR') {
    parms <- c(lambda = lambda, gamma = gamma, pop = pop)
  }
  
  # ODE solver function
  library(deSolve)
  
  if (model == 'SEIRD') {
    SPCmod <- function(t, x, parms) {
      with(as.list(c(parms, x)), {
        beta_val <- beta(t)
        dS <- -beta_val * S * I / pop
        dE <- beta_val * S * I / pop - lambda * E
        dI <- lambda * E - gamma * I - eta * I
        dR <- gamma * I
        dD <- eta * I
        list(c(dS, dE, dI, dR, dD))
      })
    }
  }
  
  if (model == 'SEIR') {
    SPCmod <- function(t, x, parms) {
      with(as.list(c(parms, x)), {
        beta_val <- beta(t)
        dS <- -beta_val * S * I / pop
        dE <- beta_val * S * I / pop - lambda * E
        dI <- lambda * E - gamma * I - eta * I
        dR <- gamma * I
        list(c(dS, dE, dI, dR))
      })
    }
  }
  
  # Solve ODE
  out <- lsoda(xstart, grid, SPCmod, parms)
  
  # Create data frame for observations
  if (model == 'SEIRD') {
    obs <- data.frame(t = grid, S = out[,"S"], E = out[,"E"], I = out[,"I"], R = out[,"R"], D = out[,"D"])
  }
  if (model == 'SEIR') {
    obs <- data.frame(t = grid, S = out[,"S"], E = out[,"E"], I = out[,"I"], R = out[,"R"])
  }
  
  # Create data frame for beta values
  beta_val <- beta(grid)
  df_beta <- data.frame(t = grid, beta = beta_val)
  
  # Create data frame for reproduction number
  df_R <- data.frame(t = grid, R = NA)
  for (i in 1:length(grid)) {
    t <- grid[i]
    S_t <- obs$S[i]
    beta_t <- beta(t)
    if (model == 'SEIR') {
      df_R$R[i] <- (beta_t / (gamma + lambda)) * (S_t / pop)
    } else if (model == 'SEIRD') {
      df_R$R[i] <- (beta_t / (gamma + lambda + eta)) * (S_t / pop)
    }
  }
  
  # Add noise if specified
  if (noise > 0) {
    
    if (model == 'SEIRD') {
      obs_with_noise <- obs[,c("S","E","I","R","D")] + matrix(rnorm(n * 5, mean = 0, sd = sqrt(noise)), nrow = n)
      obs_with_noise <- data.frame(t = grid, obs_with_noise)
      colnames(obs_with_noise) <- c("t", "S", "E", "I", "R", "D")
    }
    
    if (model == 'SEIR') {
      obs_with_noise <- obs[,c("S","E","I","R")] + matrix(rnorm(n * 4, mean = 0, sd = sqrt(noise)), nrow = n)
      obs_with_noise <- data.frame(t = grid, obs_with_noise)
      colnames(obs_with_noise) <- c("t", "S", "E", "I", "R")
    }
    
    result <- list(obs = obs, 
                   obs_with_noise = obs_with_noise, 
                   df_beta = df_beta,
                   df_R = df_R)
  } else {
    result <- list(obs = obs, 
                   df_beta = df_beta,
                   df_R = df_R)
  }
  
  return(result)
}

#' Plot simulated data, transmission rate and reproduction number
#' 
#' Plot the simulated compartment counts, the predefined transmission rate and the corresponding reproduction rate.
#' 
#' @param model String, type of model to be used for simulation ('SEIRD' or 'SEIR' available).
#' @param sim List of data frames from simulation.
#' @param latency_rate Numeric, latency rate.
#' @param recovery_rate Numeric, recovery rate.
#' @param fatality_rate Numeric, fatality rate (if provided).
#' @param log Logical indicating whether to plot the counts on a log scale.
#' 
#' @return A list containing three ggplot objects: 
#'  - the counts plot
#'  - the transmission rate plot 
#'  - and the reproduction number plot.
#' @export
plotting_simulated_data_lsoda <- function(model, sim, 
                                          latency_rate, recovery_rate,
                                          fatality_rate = NULL, 
                                          log = TRUE) {
  
  obs <- sim$obs
  df_beta <- sim$df_beta
  df_R <- sim$df_R
  
  cols_SEIRD <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442", D = "#0072B2")
  labels_SEIRD <- c(S = "Susceptible", E = "Exposed", I = "Infected", R = "Recovered", D = "Deceased")
  
  cols_SEIR <- c(S = "#E69F00", E = "#56B4E9", I = "#009E73", R = "#F0E442")
  labels_SEIR <- c(S = "Susceptible", E = "Exposed", I = "Infected", R = "Recovered")
  
  if (!is.null(sim$obs_with_noise)) { 
    obs_with_noise <- sim$obs_with_noise
  }
  
  if (model == 'SEIRD') {
    
    if (is.null(sim$obs_with_noise)) { # No obs with noise
      
      # Plot without S compartment and without log scaling
      simulated_compartments_except_S <- ggplot(obs, aes(x = t)) +
        geom_line(aes(y = E, color = 'E')) +
        geom_line(aes(y = I, color = 'I')) +
        geom_line(aes(y = R, color = 'R')) +
        geom_line(aes(y = D, color = 'D')) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations (excluding S)") +
        scale_color_manual(values = cols_SEIRD[names(cols_SEIRD) != "S"],
                           labels = labels_SEIRD[names(labels_SEIRD) != "S"], 
                           name = "Compartments") +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Fatality Rate:", round(fatality_rate,4)), hjust = 0, vjust = 1, size = 3)
      
      # To avoid warnings
      if (log == TRUE) {
        obs <- obs %>%
          mutate(
            S = ifelse(S <= 0, 0.001, S),
            E = ifelse(E <= 0, 0.001, E),
            I = ifelse(I <= 0, 0.001, I),
            R = ifelse(R <= 0, 0.001, R),
            D = ifelse(D <= 0, 0.001, D)
          )
      }
      
      # Create the counts plot
      simulated_compartments <- ggplot(obs, aes(x = t)) +
        geom_line(aes(y = S, color = 'S')) +
        geom_line(aes(y = E, color = 'E')) +
        geom_line(aes(y = I, color = 'I')) +
        geom_line(aes(y = R, color = 'R')) +
        geom_line(aes(y = D, color = 'D')) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations") +
        scale_color_manual(values = cols_SEIRD,
                           labels = labels_SEIRD, 
                           name = "Compartments") +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Fatality Rate:", round(fatality_rate,4)), hjust = 0, vjust = 1, size = 3) +
        if (log) scale_y_continuous(trans = 'log10', name = "log(compartment counts)") else NULL
      
    } else { # Obs with noise exist
      
      # Plot without S compartment and without log scaling
      simulated_compartments_except_S <- ggplot(obs, aes(x = t)) +
        geom_line(aes(y = E, color = 'E')) +
        geom_line(aes(y = I, color = 'I')) +
        geom_line(aes(y = R, color = 'R')) +
        geom_line(aes(y = D, color = 'D')) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations (excluding S)") +
        scale_color_manual(values = cols_SEIRD[names(cols_SEIRD) != "S"],
                           labels = labels_SEIRD[names(labels_SEIRD) != "S"], 
                           name = "Compartments") +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey"))) +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Fatality Rate:", round(fatality_rate,4)), hjust = 0, vjust = 1, size = 3)
      
      # To avoid warnings
      if (log == TRUE) {
        obs_with_noise <- obs_with_noise %>%
          mutate(
            S = ifelse(S <= 0, 0.001, S),
            E = ifelse(E <= 0, 0.001, E),
            I = ifelse(I <= 0, 0.001, I),
            R = ifelse(R <= 0, 0.001, R),
            D = ifelse(D <= 0, 0.001, D)
          )
      }
      
      # Create the counts plot
      simulated_compartments <- ggplot(obs, aes(x = t)) +
        geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations") +
        scale_color_manual(values = cols_SEIRD,
                           labels = labels_SEIRD, 
                           name = "Compartments") +
        scale_linetype_manual(values = c("Data" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data"),
                              guide = guide_legend(override.aes = list(color = "grey"))) +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey"))) +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Fatality Rate:", round(fatality_rate,4)), hjust = 0, vjust = 1, size = 3) +
        if (log) scale_y_continuous(trans = 'log10', name = "log(compartment counts)") else NULL
      
    }
  }
  
  if (model == 'SEIR') {
    
    if (is.null(sim$obs_with_noise)) { # No obs with noise 
      
      # Plot without S compartment and without log scaling
      simulated_compartments_except_S <- ggplot(obs, aes(x = t)) +
        geom_line(aes(y = E, color = 'E')) +
        geom_line(aes(y = I, color = 'I')) +
        geom_line(aes(y = R, color = 'R')) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations (excluding S)") +
        scale_color_manual(values = cols_SEIR[names(cols_SEIR) != "S"],
                           labels = labels_SEIR[names(labels_SEIR) != "S"], 
                           name = "Compartments") +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3)
      
      # To avoid warnings
      if (log == TRUE) {
        obs <- obs %>%
          mutate(
            S = ifelse(S <= 0, 0.001, S),
            E = ifelse(E <= 0, 0.001, E),
            I = ifelse(I <= 0, 0.001, I),
            R = ifelse(R <= 0, 0.001, R),
          )
      }
      
      # Create the counts plot
      simulated_compartments <- ggplot(obs, aes(x = t)) +
        geom_line(aes(y = S, color = 'S')) +
        geom_line(aes(y = E, color = 'E')) +
        geom_line(aes(y = I, color = 'I')) +
        geom_line(aes(y = R, color = 'R')) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations") +
        scale_color_manual(values = cols_SEIR,
                           labels = labels_SEIR, 
                           name = "Compartments") +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3) +
        if (log) scale_y_continuous(trans = 'log10', name = "log(compartment counts)") else NULL
      
    } else { # Obs with noise exists
      
      # Plot without S compartment and without log scaling
      simulated_compartments_except_S <- ggplot(obs, aes(x = t)) +
        geom_line(aes(y = E, color = 'E')) +
        geom_line(aes(y = I, color = 'I')) +
        geom_line(aes(y = R, color = 'R')) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations (excluding S)") +
        scale_color_manual(values = cols_SEIR[names(cols_SEIR) != "S"],
                           labels = labels_SEIR[names(labels_SEIR) != "S"], 
                           name = "Compartments") +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey"))) +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3)
      
      # To avoid warnings
      if (log == TRUE) {
        obs_with_noise <- obs_with_noise %>%
          mutate(
            S = ifelse(S <= 0, 0.001, S),
            E = ifelse(E <= 0, 0.001, E),
            I = ifelse(I <= 0, 0.001, I),
            R = ifelse(R <= 0, 0.001, R)
          )
      }
      
      # Create the counts plot
      simulated_compartments <- ggplot(obs, aes(x = t)) +
        geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data")) +
        geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
        geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
        theme_minimal() +
        ylab("Compartment counts") +
        xlab("Time") + 
        ggtitle("Simulated observations") +
        scale_color_manual(values = cols_SEIR,
                           labels = labels_SEIR, 
                           name = "Compartments") +
        scale_linetype_manual(values = c("Data" = "solid"),
                              name = "Lines",
                              labels = c("Simulated data"),
                              guide = guide_legend(override.aes = list(color = "grey"))) +
        scale_shape_manual(values = c("Noisy data" = 1),
                           name = "Points",
                           labels = c("Noisy data"),
                           guide = guide_legend(override.aes = list(color = "grey"))) +
        guides(color = guide_legend(order = 1)) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Latency Rate:", round(latency_rate,4)), hjust = 0, vjust = -3, size = 3) +
        annotate("text", x = min(obs$t), y = max(obs$E), label = paste("Recovery Rate:", round(recovery_rate,4)), hjust = 0, vjust = -1, size = 3) +
        if (log) scale_y_continuous(trans = 'log10', name = "log(compartment counts)") else NULL
      
    }
    
  }
  
  # Create the transmission rate plot
  simulated_beta <- ggplot() +
    geom_line(data = df_beta, aes(x = t, y = beta, color = "Simulated transmission rate"), linetype = "dashed") +
    theme_minimal() +
    labs(x = "Time", 
         y = TeX("Transmission rate $\\beta(t)$"), 
         title = TeX("Simulated transmission rate $\\beta(t)$"),
         color = "Legend") + 
    scale_color_manual(values = c("Simulated transmission rate" = "#E69F00"),
                       labels = c("Simulated transmission rate"), name = "Line") +
    theme(legend.position = "top")
  
  # Create the reproduction number plot
  simulated_R <- ggplot() +
    geom_line(data = df_R, aes(x = t, y = R, color = "Simulated reproduction number"), linetype = "dashed") +
    theme_minimal() +
    geom_hline(yintercept = 1, color = "#0072B2", linetype = "solid", linewidth = 0.5, alpha = 0.5) +
    labs(x = "Time", 
         y = TeX("Reproduction number $R(t)$"), 
         title = TeX("Simulated reproduction number $R(t)$"),
         color = "Legend") +   
    scale_color_manual(values = c("Simulated reproduction number" = "#56B4E9"),
                       labels = c("Simulated reproduction number"), name = "Line") +
    theme(legend.position = "top")
  
  # Return the plots as a list
  list(simulated_compartments = simulated_compartments, 
       simulated_compartments_except_S = simulated_compartments_except_S, 
       simulated_beta = simulated_beta, 
       simulated_R = simulated_R)
}
