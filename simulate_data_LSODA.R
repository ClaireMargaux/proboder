###########################################
####### SIMULATE DATA USING LSODA #########
###########################################

#' Simulate data using LSODA
#' 
#' Simulate compartment counts (with or without noise) using the function lsoda from the package deSolve.
#' 
#' @param noise Variance of the Gaussian white noise to be added to the data (if set to 0, no noise is added).
#' @param seed Seed to be used for the noise added to the data.
#' @param steps Size of the time steps of the generated data.
#' @param max_time Last time step of the generated data.
#' @param pop Population size of generated data.
#' @param beta Function; contact rate to be used for simulation.
#' @param xstart Starting values for the compartments.
#' @param lambda Latency rate.
#' @param gamma Recovery rate.
#' @param eta Fatality rate.
#' 
#' @return List containing two data frames: (1) (noisy) simulated counts and (2) predefined contact rate.
#' @export
simulate_data_LSODA <- function(noise = 0, seed = 5, steps = 1, max_time = 30, 
                                lambda, gamma, eta, pop, beta, xstart) {
  
  # Set seed for reproducibility
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  # Create time grid
  grid <- seq(0, max_time, by = steps)
  n <- length(grid)
  
  # Define parameters
  parms <- c(lambda = lambda, gamma = gamma, eta = eta, pop = pop)
  
  # ODE solver function
  library(deSolve)
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
  
  # Solve ODE
  out <- lsoda(xstart, grid, SPCmod, parms)
  
  # Create data frame for observations
  obs <- data.frame(t = grid, S = out[,"S"], E = out[,"E"], I = out[,"I"], R = out[,"R"], D = out[,"D"])
  
  # Create data frame for beta values
  beta_val <- beta(grid)
  df_beta <- data.frame(t = grid, beta = beta_val)
  
  # Add noise if specified
  if (noise > 0) {
    obs_with_noise <- obs[,c("S","E","I","R","D")] + matrix(rnorm(n * 5, mean = 0, sd = sqrt(noise)), nrow = n)
    obs_with_noise <- data.frame(t = grid, obs_with_noise)
    colnames(obs_with_noise) <- c("t", "S", "E", "I", "R", "D")
    result <- list(obs_with_noise, df_beta)
  } else {
    result <- list(obs, df_beta)
  }
  
  return(result)
}
  
  # # Beta function: sine
  # beta <- function(t){0.05*sin(t/10)+0.8}

  # # Beta function: log
  # beta <- function(t){0.03*log(t+1)+0.8}

#' Plot Simulated Data and Beta Values
#' 
#' Plot the simulated compartment counts and the predefined contact rate.
#' 
#' @param obs Data frame containing the observed counts.
#' @param df_beta Data frame containing the predefined contact rate.
#' @param log Logical indicating whether to plot the counts on a log scale.
#' 
#' @return A list containing two ggplot objects: the counts plot and the beta plot.
#' @export
plotting_simulated_data_lsoda <- function(obs, df_beta, log = TRUE) {
  
  # Create the counts plot
  countsplot <- ggplot(obs, aes(x = t)) +
    geom_line(aes(y = S, color = 'S')) +
    geom_line(aes(y = E, color = 'E')) +
    geom_line(aes(y = I, color = 'I')) +
    geom_line(aes(y = R, color = 'R')) +
    geom_line(aes(y = D, color = 'D')) +
    theme_minimal() +
    ylab("count") +
    xlab("time") + 
    ggtitle("Observations without noise") +
    scale_color_manual(values = c('S' = '#E69F00', 'E' = '#56B4E9', 'I' = '#009E73', 'R' = '#F0E442', 'D' = '#0072B2')) +
    if (log) scale_y_continuous(trans = 'log10', name = "log(count)") else NULL
  
  # Create the beta plot
  betaplot <- ggplot(df_beta, aes(x = t, y = beta)) +
    geom_line(color = '#E69F00') +
    theme_minimal() +
    ylab("beta") +
    xlab("time") +
    ggtitle("Simulated beta") 
  
  # Return the plots as a list
  list(countsplot, betaplot)
}

################################
######## SAVE THE DATA #########
################################

# save(df, file = file.path(directory, "simulated_data_LSODA.Rdata"))
# save(df_with_noise, file = file.path(directory, "simulated_noisy_data_LSODA.Rdata"))
# save(df_params, file = file.path(directory, "simulated_params_LSODA.Rdata"))
# saveRDS(beta_val, file = file.path(directory, "simulated_beta_LSODA.Rds"))
# save(df_beta, file = file.path(directory, "simulated_beta_LSODA.Rdata"))