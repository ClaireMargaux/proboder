#################################### PROBODER ##################################
################################ Claire Descombes ##############################
########################## Saving, loading, plotting ########################### 

load_data <- function(type, region, daily_or_weekly, directory) {
  # Check if type is valid
  if (!(type %in% c('simulated', 'real'))) {
    stop("Error: Invalid type. Type must be 'simulated' or 'real'.")
  }
  
  # Check if region is valid
  if (!(region %in% c('BE', 'GE'))) {
    stop("Error: Invalid region. Region must be 'BE' or 'GE'.")
  }
  
  # Check if daily_or_weekly is valid
  if (!(daily_or_weekly %in% c('daily', 'weekly'))) {
    stop("Error: Invalid value for daily_or_weekly. Must be 'daily' or 'weekly'.")
  }
  
  # Load data based on type
  if (type == 'simulated') {
    load(file.path(directory, "simulated_data.Rdata"))
    obs <- observations_simulate
    population <- readRDS(file.path(directory, "simulated_pop.Rds"))
    real_beta <- readRDS(file.path(directory, "simulated_real_beta.Rds"))
  } else {
    region_filename <- paste0("real_data_", region, "_", daily_or_weekly, ".Rdata")
    population_filename <- paste0("real_pop_", region, "_", daily_or_weekly, ".Rds")
    load(file.path(directory, region_filename))
    obs <- observations
    population <- readRDS(file.path(directory, population_filename))
    real_beta <- NULL 
  }
  
  # Return loaded data
  return(list(observations = obs, population = population, real_beta = real_beta))
}

save_matrices_as_Rdata <- function(X_values, U_values, P_X_values, P_U_values, directory) {
  # Save X_values, U_values, P_X_values, and P_U_values as .Rdata files
  save(X_values, file = paste0(directory, "/X_values.Rdata"))
  save(U_values, file = paste0(directory, "/U_values.Rdata"))
  save(P_X_values, file = paste0(directory, "/P_X_values.Rdata"))
  save(P_U_values, file = paste0(directory, "/P_U_values.Rdata"))
  
  # Calculate sigmoid of U_values[1, ]
  U_scaled <- round(sigmoid(U_values[1, ]),3)
  
  # Save sigmoid of U_values[1, ] as CSV
  U_scaled_df <- data.frame(U_scaled)
  write.csv(U_scaled_df, file = paste0(directory, "/U_scaled.csv"), row.names = FALSE)
}

load_and_process_data <- function(directory_res, time_grid) {
  # Load matrices from saved .Rdata files
  load(file.path(directory_res, "X_values.Rdata"))
  load(file.path(directory_res, "U_values.Rdata"))
  load(file.path(directory_res, "P_X_values.Rdata"))
  load(file.path(directory_res, "P_U_values.Rdata"))
  
  # Process U_values for visualization
  U_scaled <- sigmoid(U_values[1, ])
  U_plot <- data.frame(time = time_grid, U_scaled = U_scaled, row.names = NULL)
  
  # Process P_U_values for visualization
  P_U <- P_U_values[1, 1, ]
  
  # Generate values for error area
  ymin <- mapply(function(mu, sigma) sigmoid(qnorm(0.025, mean = mu, sd = sqrt(sigma))), U_scaled, P_U)
  ymax <- mapply(function(mu, sigma) sigmoid(qnorm(0.975, mean = mu, sd = sqrt(sigma))), U_scaled, P_U)
  
  # Create data frame for P_plot
  P_plot <- data.frame(time = time_grid, ymin = ymin, ymax = ymax, row.names = NULL)
  
  # Return processed data
  return(list(U_plot = U_plot, P_plot = P_plot, U_scaled = U_scaled))
}

plot_contact_rate <- function(type, U_plot, ymin, ymax, U_scaled, real_beta_df = NULL) {
  library(ggplot2)
  
  if (type == 'simulated') {
    ggplot() +
      geom_line(data = U_plot, aes(x = time, y = U_scaled, color = "Estimated Contact Rate"), linewidth = 1) +
      geom_line(data = real_beta_df, aes(x = time, y = real_beta, color = "Real Contact Rate"), linetype = "dashed") +
      geom_ribbon(data = data.frame(time = U_plot$time, ymin = ymin, ymax = ymax), aes(x = time, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      labs(x = "Time", y = "Contact rate", title = "Contact rate with 95%-confidence interval",
           color = "Legend") +  
      scale_color_manual(values = c("Estimated Contact Rate" = "darkgreen", "Real Contact Rate" = "lightblue4"),
                         labels = c("Estimated Contact Rate", "Real Contact Rate"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "lightgreen"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2))
  } else {
    ggplot() +
      geom_line(data = U_plot, aes(x = time, y = U_scaled, color = "Estimated Contact Rate"), linewidth = 1) +
      geom_ribbon(data = data.frame(time = U_plot$time, ymin = ymin, ymax = ymax), aes(x = time, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      labs(x = "Time", y = "Contact rate", title = "Contact rate with 95%-confidence interval",
           color = "Legend") +  
      scale_color_manual(values = c("Estimated Contact Rate" = "darkgreen"),
                         labels = "Estimated Contact Rate", name = "Line") +  
      scale_fill_manual(values = c("Error Area" = "lightgreen"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top")+
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2))
  }
}

plot_data <- function(obs, type) {
  # Adjust variable names and labels based on type
  if (type == 'simulated') {
    S_var <- 'S'
    I_var <- 'I'
    R_var <- 'R'
    legend_labels <- c("Susceptible", "Infected", "Recovered")
    plot_title <- "S, I and R counts"
  } else if (type == 'real') {
    S_var <- 'S'
    I_var <- 'I'
    R_var <- 'D'  # Assuming 'D' represents 'Deaths' in real data
    legend_labels <- c("Susceptible", "Infected", "Deaths")
    plot_title <- "S, I and D counts"
  } else {
    stop("Invalid type. Please specify 'simulated' or 'real'.")
  }
  
  obs[[S_var]] <- obs[[S_var]] + 0.01
  obs[[I_var]] <- obs[[I_var]] + 0.01
  obs[[R_var]] <- obs[[R_var]] + 0.01

  cols <- c("S" = "darkblue", 
            "I" = "coral3", 
            "R" = "darkgreen")
    
  ggplot() +
    geom_line(data = obs, aes(x = date, y = .data[[S_var]], color = "S"), linewidth = 1) +
    geom_line(data = obs, aes(x = date, y = .data[[I_var]], color = "I"), linewidth = 1) +
    geom_line(data = obs, aes(x = date, y = .data[[R_var]], color = "R"), linewidth = 1) +
    scale_y_log10() +  # Log scale for y-axis
    labs(x = "Date", y = "Count (log-scale)", title = plot_title) +
    theme_minimal() +
    scale_color_manual(values = cols,
                       name = "Compartment")
}
