#' Perform Random Search for Hyperparameter Optimization
#'
#' This function performs a random search over specified hyperparameter ranges for a given model and returns a summary of the results.
#' It evaluates the model using different sets of hyperparameters and computes mean and median scores (SPE, NLPD, CRPS) for each set.
#'
#' @param model The model to be used for inference ('SEIRD' and 'SEIR' available).
#' @param obs The observed data.
#' @param df_beta Data frame containing contact rate values.
#' @param noise Observation noise level.
#' @param lambda Latency rate.
#' @param gamma Recovery rate.
#' @param eta Fatality rate.
#' @param pop Population size.
#' @param beta0 Initial contact rate value. 
#' @param jit Boolean, set TRUE to add a jitter to the innovation covariance (optional).
#' @param seed Seed for random number generation. Default is 5.
#' @param num_param_sets Number of parameter sets to sample.
#' @param seq_l Sequence of values for the `l` hyperparameter. 
#' @param seq_wiener_X Sequence of values for the `noise_wiener_X` hyperparameter.
#' @param seq_wiener_U Sequence of values for the `noise_wiener_U` hyperparameter.
#' @param seq_beta0prime Sequence of values for the `beta0prime` hyperparameter. 
#'
#' @return A data frame containing the summary of results for each sampled set of hyperparameters. The data frame includes the hyperparameter values, mean and median scores for SPE, NLPD, and CRPS.
#' @export
run_random_search <- function(model, obs_to_use, df_beta, noise, 
                              lambda, gamma, eta, pop, beta0, jit = FALSE,
                              seed = 5, 
                              num_param_sets, 
                              seq_l, 
                              seq_wiener_X, 
                              seq_wiener_U, 
                              seq_beta0prime) {
 
  # Start counting computation time
  tic("Duration of the whole workflow")
  
  param_grid <- expand.grid(
    l = seq_l,
    noise_wiener_X = seq_wiener_X,
    noise_wiener_U = seq_wiener_U,
    beta0prime = seq_beta0prime
  )
  
  set.seed(seed) 
  sampled_grid <- param_grid %>% sample_n(num_param_sets)
  
  num_params <- dim(param_grid)[2]
  num_predictions <- length(obs$t)
  
  # Matrix for results
  results_grid_search <- matrix(
    nrow = num_param_sets,
    ncol = (num_params + 3 * num_predictions)
  )
  
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = nrow(sampled_grid), clear = FALSE, width = 60
  )
  
  # Generate time grids for inference
  grids <- generate_grid(obs_to_use, num_points_between = 0)
  
  for (i in 1:num_param_sets) {
    param_i <- sampled_grid[i, ]
    
    # Initialize with current set of parameters
    initial_params <- initialization(
      model = model,
      obs = obs_to_use, 
      beta0 = beta0, 
      beta0prime = param_i$beta0prime,
      lambda = lambda, 
      gamma = gamma, 
      eta = eta,
      l = param_i$l, 
      scale = 1, 
      noise_obs = noise,
      noise_X = sqrt(noise), 
      noise_U = 0.01,
      noise_wiener_X = param_i$noise_wiener_X, 
      noise_wiener_U = param_i$noise_wiener_U,
      pop = pop, 
      num_initial_values = 1
    )
    
    # Run inference
    inference_results <- inference(model = model, 
                                   grids = grids, 
                                   obs = obs_to_use, 
                                   jit = jit,
                                   initial_params = initial_params)
    
    # Process data for scoring and visualization
    processed_data <- process_data(inference_results = inference_results, 
                                   grids = grids)
    U_plot <- processed_data$U_plot
    
    # Compute the different scores
    score_SPE <- squared_prediction_error(U_plot, df_beta)
    score_NLPD <- negative_log_predictive_density(U_plot, df_beta)
    score_CRPS <- continuous_ranked_probability_score(U_plot, df_beta)
    
    # Store results (during grid search)
    results_grid_search[i, ] <- c(param_i$l, param_i$noise_wiener_X, param_i$noise_wiener_U, param_i$beta0prime, score_SPE, score_NLPD, score_CRPS)

    # Update progress bar
    pb$tick()
    
  } # End of grid search
  
  # Extract the parameter values
  params <- results_grid_search[, 1:num_params]
  
  # Extract the score vectors
  SPE_scores <- results_grid_search[, (num_params + 1):(num_params + num_predictions)]
  NLPD_scores <- results_grid_search[, (num_params + num_predictions + 1):(num_params + 2 * num_predictions)]
  CRPS_scores <- results_grid_search[, (num_params + 2 * num_predictions + 1):(num_params + 3 * num_predictions)]
  
  # Calculate the mean and median for each score
  mean_SPE <- rowMeans(SPE_scores)
  median_SPE <- apply(SPE_scores, 1, median)
  
  mean_NLPD <- rowMeans(NLPD_scores)
  median_NLPD <- apply(NLPD_scores, 1, median)
  
  mean_CRPS <- rowMeans(CRPS_scores)
  median_CRPS <- apply(CRPS_scores, 1, median)
  
  # Combine into a new data frame
  results_summary <- data.frame(
    l = params[, 1],
    noise_wiener_X = params[, 2],
    noise_wiener_U = params[, 3],
    beta0prime = params[, 4],
    mean_SPE = mean_SPE,
    median_SPE = median_SPE,
    mean_NLPD = mean_NLPD,
    median_NLPD = median_NLPD,
    mean_CRPS = mean_CRPS,
    median_CRPS = median_CRPS
  )

  # Total computation time
  toc()
  
  return(results_summary)
}

#' Plot Mean and Median Scores for Hyperparameter Search
#'
#' This function generates and returns ggplot objects for visualizing the mean and median scores
#' (SPE, NLPD, CRPS) across different hyperparameters.
#'
#' @param results_summary A data frame containing the summary of results from the hyperparameter search.
#' Each row represents a different set of hyperparameters and contains columns for the mean and median
#' scores (`mean_SPE`, `mean_NLPD`, `mean_CRPS`, `median_SPE`, `median_NLPD`, `median_CRPS`).
#' @param parameters_to_plot A character vector of parameter names for which the mean and median scores
#' should be plotted. Default is `c('l', 'noise_wiener_X', 'noise_wiener_U', 'beta0prime')`.
#' @return A list containing two lists: `mean_plots` and `median_plots`. Each list contains ggplot objects
#' for the specified parameters.
#' @export
plot_scores <- function(results_summary, 
                        parameters_to_plot = 
                          c('l', 'noise_wiener_X', 
                            'noise_wiener_U', 'beta0prime')) {
  
  # Function to plot mean scores for a given parameter
  plot_mean_scores <- function(df, parameter) {
    df %>%
      ggplot(aes_string(x = parameter)) +
      geom_point(aes(y = mean_SPE, color = "SPE")) +
      geom_point(aes(y = mean_NLPD, color = "NLPD")) +
      geom_point(aes(y = mean_CRPS, color = "CRPS")) +
      labs(title = paste("Mean scores for parameter", parameter), x = parameter, y = "Mean score") +
      scale_color_manual(values = c("SPE" = "#E69F00", "NLPD" = "#56B4E9", "CRPS" = "#009E73")) +
      theme_minimal()
  }
  
  # Function to plot median scores for a given parameter
  plot_median_scores <- function(df, parameter) {
    df %>%
      ggplot(aes_string(x = parameter)) +
      geom_point(aes(y = median_SPE, color = "SPE")) +
      geom_point(aes(y = median_NLPD, color = "NLPD")) +
      geom_point(aes(y = median_CRPS, color = "CRPS")) +
      labs(title = paste("Median scores for parameter", parameter), x = parameter, y = "Median score") +
      scale_color_manual(values = c("SPE" = "#E69F00", "NLPD" = "#56B4E9", "CRPS" = "#009E73")) +
      theme_minimal()
  }
  
  # Create a list to store the plots
  plot_list_mean <- list()
  plot_list_median <- list()
  
  # Generate and store plots for each parameter
  for (param in parameters_to_plot) {
    plot_list_mean[[param]] <- plot_mean_scores(results_summary, param)
    plot_list_median[[param]] <- plot_median_scores(results_summary, param)
  }
  
  # Return the lists of plots
  return(list(mean_plots = plot_list_mean, median_plots = plot_list_median))
}

#' Find the Best Set of Hyperparameters
#'
#' This function identifies the best set of hyperparameters from a summary of results based on mean scores for SPE, NLPD, and CRPS.
#' It first checks for a set that minimizes all three scores, then checks for sets minimizing at least two scores, and finally prioritizes the minimal CRPS score if no other criteria are met.
#'
#' @param results_summary A data frame containing the summary of results with mean scores (SPE, NLPD, CRPS) for each set of hyperparameters.
#'
#' @return A data frame containing the best set of hyperparameters based on the criteria outlined.
#'
#' @export
find_best_hyperparameters <- function(results_summary) {
  # Find the best set of hyperparameters for each score
  best_spe <- results_summary[which.min(results_summary$mean_SPE), ]
  best_nlpd <- results_summary[which.min(results_summary$mean_NLPD), ]
  best_crps <- results_summary[which.min(results_summary$mean_CRPS), ]
  
  # Check if there is a set minimizing all three scores
  best_all <- results_summary %>%
    filter(mean_SPE == min(mean_SPE) & mean_NLPD == min(mean_NLPD) & mean_CRPS == min(mean_CRPS))
  
  if (nrow(best_all) > 0) {
    best_params <- best_all[1, ]
  } else {
    # Find sets minimizing at least two of the scores
    best_two <- results_summary %>%
      filter((mean_SPE == min(mean_SPE) & mean_NLPD == min(mean_NLPD)) |
               (mean_SPE == min(mean_SPE) & mean_CRPS == min(mean_CRPS)) |
               (mean_NLPD == min(mean_NLPD) & mean_CRPS == min(mean_CRPS)))
    
    if (nrow(best_two) > 0) {
      best_params <- best_two[1, ]
    } else {
      # If no sets minimize at least two scores, prioritize minimal CRPS
      best_params <- best_crps
    }
  }
  
  return(best_params)
}
