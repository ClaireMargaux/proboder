#' Perform random search for hyperparameter optimization
#'
#' This function performs a random search over specified hyperparameter ranges for a given model and returns a summary of the results.
#' It evaluates the model using different sets of hyperparameters and computes mean scores (SPE, NLPD, CRPS) for each set.
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
#' @param jit_X Jitter on the variance of X (optional, default 1e-9).
#' @param jit_U Jitter on the variance of U (optional, default 1e-9).
#' @param seed Seed for random number generation. Default is 5.
#' @param num_param_sets Number of parameter sets to sample.
#' @param seq_l Sequence of values for the `l` hyperparameter. 
#' @param seq_wiener_X Sequence of values for the `noise_wiener_X` hyperparameter.
#' @param seq_wiener_U Sequence of values for the `noise_wiener_U` hyperparameter.
#' @param seq_beta0prime Sequence of values for the `beta0prime` hyperparameter. 
#' @param seq_jit Sequence of values for the `jit` hyperparameter. 
#'
#' @return A data frame containing the summary of results for each sampled set of hyperparameters. The data frame includes the hyperparameter values and mean scores for SPE, NLPD, and CRPS.
#' @export
run_random_search <- function(model, obs_to_use, df_beta, noise_obs, 
                              lambda, gamma, eta, pop, beta0,
                              jit_X = 1e-9, jit_U = 1e-9,
                              seed = 5, 
                              num_param_sets, 
                              seq_l, 
                              seq_wiener_X, 
                              seq_wiener_U, 
                              seq_beta0prime,
                              seq_jit) {
 
  # Start counting computation time
  tic("Duration of the whole workflow")
  
  param_grid <- expand.grid(
    l = seq_l,
    noise_wiener_X = seq_wiener_X,
    noise_wiener_U = seq_wiener_U,
    beta0prime = seq_beta0prime,
    jit = seq_jit
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
  
  # Initialize a counter for warnings
  warning_count <- 0
  
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
      noise_obs = noise_obs,
      noise_X = sqrt(noise_obs), 
      noise_U = 0.2,
      noise_wiener_X = param_i$noise_wiener_X, 
      noise_wiener_U = param_i$noise_wiener_U,
      pop = pop, 
      num_initial_values = 1
    )
    
    # Run inference with error handling
    inference_results <- tryCatch({
      inference(model = model, 
                grids = grids, 
                obs = obs_to_use, 
                jit = param_i$jit,
                jit_X = jit_X, jit_U = jit_U,
                initial_params = initial_params)
    },
    warning = function(w) {
      warning_count <<- warning_count + 1
      # Return NULL or an appropriate placeholder
      return(NULL)
    },
    error = function(e) {
      warning_count <<- warning_count + 1
      # Return NULL or an appropriate placeholder
      return(NULL)
    }
    )
    
    # Check if inference was successful
    if (!is.null(inference_results)) {
      # Process data for scoring and visualization
      processed_data <- process_data(inference_results = inference_results,
                                     grids = grids,
                                     model = model,
                                     pop = pop, 
                                     lambda = lambda, 
                                     gamma = gamma, 
                                     eta = eta)
      U_plot <- processed_data$U_plot
      
      # Compute the different scores
      score_CRPS <- continuous_ranked_probability_score(U_plot, df_beta)
      score_NLPD <- negative_log_predictive_density(U_plot, df_beta)
      score_SPE <- squared_prediction_error(U_plot, df_beta)
      
      # Store results (during grid search)
      results_grid_search[i, ] <- c(param_i$l, param_i$noise_wiener_X, param_i$noise_wiener_U, param_i$beta0prime, param_i$jit, param_i$clip_X, param_i$clip_U, score_SPE, score_NLPD, score_CRPS)
    } else {
      # Handle the case where inference_results is NULL
      results_grid_search[i, 1:num_params] <- c(param_i$l, param_i$noise_wiener_X, param_i$noise_wiener_U, param_i$beta0prime, param_i$jit, param_i$clip_X, param_i$clip_U)
    }
    
    # Update progress bar
    pb$tick()
    
  } # End of grid search
  
  # Print the number of warnings encountered
  message("Number of warnings encountered: ", warning_count)
  
  # Extract the parameter values
  params <- results_grid_search[, 1:num_params]
  
  # Extract the score vectors
  SPE_scores <- results_grid_search[, (num_params + 1):(num_params + num_predictions)]
  NLPD_scores <- results_grid_search[, (num_params + num_predictions + 1):(num_params + 2 * num_predictions)]
  CRPS_scores <- results_grid_search[, (num_params + 2 * num_predictions + 1):(num_params + 3 * num_predictions)]
  
  # Calculate the mean and median for each score
  mean_CRPS <- rowMeans(CRPS_scores)
  mean_NLPD <- rowMeans(NLPD_scores)
  mean_SPE <- rowMeans(SPE_scores)
  
  # Combine into a new data frame
  results_summary <- data.frame(
    l = params[, 1],
    noise_wiener_X = params[, 2],
    noise_wiener_U = params[, 3],
    beta0prime = params[, 4],
    jit = params[, 5],
    mean_CRPS = mean_CRPS,
    mean_NLPD = mean_NLPD,
    mean_SPE = mean_SPE
  )

  # Total computation time
  toc()
  
  return(results_summary)
}

#' Plot mean scores for hyperparameter search
#'
#' This function generates and returns ggplot objects for visualizing the mean scores
#' (SPE, NLPD, CRPS) across different hyperparameters. The function can optionally aggregate scores
#' for each hyperparameter value or plot raw scores with NA values removed.
#'
#' @param results_summary A data frame containing the summary of results from the hyperparameter search.
#' Each row represents a different set of hyperparameters and contains columns for the mean 
#' scores (`mean_SPE`, `mean_NLPD`, `mean_CRPS`).
#' @param parameters_to_plot A character vector of parameter names for which the mean scores
#' should be plotted. Default is `c('l', 'noise_wiener_X', 'noise_wiener_U', 'beta0prime', 'jit')`.
#' @param aggregate A logical value indicating whether to aggregate scores by the parameter values.
#' If `TRUE`, scores are aggregated by taking the mean for each parameter value. If `FALSE`, 
#' raw scores with NA values removed are plotted. Default is `FALSE`.
#' @return A list `mean_plots` containing ggplots of mean scores.
#' If `aggregate` is `TRUE`, the plots reflect aggregated scores; 
#' otherwise, they reflect raw scores with NA values removed.
#' @export
plot_scores <- function(results_summary, 
                        parameters_to_plot = 
                          c('l', 'noise_wiener_X', 
                            'noise_wiener_U', 'beta0prime', 'jit'),
                        aggregate = FALSE) {
  
  # Define the LaTeX expressions
  params_expressions <- c(
    "\u2113",  # Unicode for ell
    "$Q_X$",
    "$Q_U$",
    "$m_{U^{(1)}}$",
    "jit"
  )
  
  # Create a named vector for easier lookup
  names(params_expressions) <- parameters_to_plot
  
  # Define color scale
  color_scale <- scale_color_manual(
    values = c("SPE" = "#E69F00", 
               "NLPD" = "#56B4E9", 
               "CRPS" = "#009E73"),
             name = "Score")
  
  # Function to aggregate data
  aggregate_data_mean <- function(df, parameter) {
    df %>%
      group_by(across(all_of(parameter))) %>%  # Group by the parameter
      summarize(
        mean_SPE = mean(mean_SPE, na.rm = TRUE),
        mean_NLPD = mean(mean_NLPD, na.rm = TRUE),
        mean_CRPS = mean(mean_CRPS, na.rm = TRUE),
        .groups = 'drop'  # Drop the grouping after summarizing
      )
  }
  
  # Function to plot mean scores for a given parameter
  plot_mean_scores <- function(df, parameter, aggregate) {
    # Aggregate data if needed
    if (aggregate) {
      aggregated_df <- aggregate_data_mean(df, parameter)
      title_text <- TeX(paste("Aggregated mean scores for", params_expressions[parameter]))
    } else {
      aggregated_df <- df %>%
        filter(!is.na(mean_SPE) & !is.na(mean_NLPD) & !is.na(mean_CRPS))
      title_text <- TeX(paste("Mean scores for", params_expressions[parameter]))
    }
    
    p <- ggplot(aggregated_df, aes(x = !!sym(parameter))) +
      geom_point(aes(y = mean_SPE, color = "SPE")) +
      geom_point(aes(y = mean_NLPD, color = "NLPD")) +
      geom_point(aes(y = mean_CRPS, color = "CRPS")) +
      labs(title = title_text, 
           x = TeX(params_expressions[parameter]), 
           y = "Mean score",
           color = "Score") +
      color_scale +  # Use predefined color scale
      theme_minimal()
    
    # Adjust x-axis scale for non-linear parameters
    if (parameter %in% c('noise_wiener_X', 'jit')) {
      p <- p + scale_x_log10()  # Apply log scale to x-axis
    }
    
    return(p)
  }
  
  # Create a list to store the plots
  plot_list_mean <- list()

  # Generate and store plots for each parameter
  for (param in parameters_to_plot) {
    plot_list_mean[[param]] <- plot_mean_scores(results_summary, param, aggregate)
  }
  
  # Return the lists of plots
  return(list(mean_plots = plot_list_mean))
}

#' Return table with best sets of hyperparameters for each score.
#'
#' This function returns a table with the best set of hyperparameters for the mean scores for CRPS, NLPD and SPE.
#'
#' @param best_parameters A list of a
#' - kable table, and a
#' - data frame
#' containing the summary of results with mean scores (CRPS, NLPD, SPE) for each set of hyperparameters.
#'
#' @return A data frame containing the best sets of hyperparameters based on the criteria outlined, including:
#' \item{Best SPE (Mean)}{The hyperparameters that yield the minimum mean SPE.}
#' \item{Best NLPD (Mean)}{The hyperparameters that yield the minimum mean NLPD.}
#' \item{Best CRPS (Mean)}{The hyperparameters that yield the minimum mean CRPS.}
#'
#' @export
find_best_hyperparameters <- function(results_summary) {
  
  # Find the best set of hyperparameters for each score (mean)
  best_spe_mean <- results_summary[which.min(results_summary$mean_SPE), ]
  best_nlpd_mean <- results_summary[which.min(results_summary$mean_NLPD), ]
  best_crps_mean <- results_summary[which.min(results_summary$mean_CRPS), ]

  # Extract only parameter columns from results_summary
  param_cols <- c("l", "noise_wiener_X", "noise_wiener_U", "beta0prime", "jit")
  
  # Function to create parameter rows
  create_param_row <- function(type_label, params) {
    params$Metric <- type_label
    params <- params %>%
      select(Metric, all_of(param_cols))  # Select only parameter columns
    return(params)
  }
  
  # Combine results into a single data frame with specific order
  combined_results <- rbind(
    create_param_row("Best mean CRPS", best_crps_mean),
    create_param_row("Best mean NLPD", best_nlpd_mean),
    create_param_row("Best mean SPE", best_spe_mean)
  )
  
  # Reorder columns to have 'Metric' as the first column
  combined_results <- combined_results %>%
    select(Metric, everything())
  
  # Generate table using knitr and kableExtra
  best_params_table <- kable(combined_results, align = "c", 
                             caption = "Best hyperparameters for each score") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) %>%
    column_spec(1, bold = TRUE) 
  
  best_parameters <- list(best_params_table = best_params_table, 
                          best_params_df = combined_results)
  
  return(best_parameters)
}

#' Generate and plot scores for the best hyperparameter sets
#'
#' This function evaluates three sets of hyperparameters identified as the best based on
#' their mean Continuous Ranked Probability Score (CRPS), Negative Log Predictive Density (NLPD),
#' and Standard Prediction Error (SPE) scores. The function runs inference using these hyperparameters,
#' computes the corresponding scores, and generates plots comparing the CRPS, NLPD, and SPE scores
#' across the different hyperparameter sets.
#'
#' @param best_parameters A list containing the best hyperparameters, generated by the 
#'        `find_best_hyperparameters` function. This list includes a data frame (`best_params_df`) 
#'        of the best hyperparameter sets.
#' @param model The model object used for inference. It contains the structure and parameters of the model.
#' @param obs_to_use A data frame of observed data used in the inference process.
#' @param df_beta A data frame containing transmission rate values used in score computation.
#' @param noise_obs A numeric value representing the observational noise level.
#' @param lambda, gamma, eta Numeric values representing parameters of the model used in the inference process.
#' @param pop A numeric value representing the population size used in the model.
#' @param seed An optional seed value for reproducibility of results.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{best_scores_df}{A data frame with the CRPS, NLPD, and SPE scores for each hyperparameter set.}
#'   \item{best_scores_plot}{A list of three ggplot objects, one for each score (CRPS, NLPD, SPE), comparing 
#'        the scores across the different hyperparameter sets.}
#' }
#'
#' @details
#' The function first initializes the model with each set of hyperparameters and runs inference. 
#' It then processes the data to compute the CRPS, NLPD, and SPE scores. Finally, it creates plots 
#' for each score type, showing the scores for each hyperparameter set.
#' 
#' @export
scores_of_best_params <- 
  function(best_parameters, 
           model, obs_to_use, 
           df_beta, noise_obs, 
           lambda, gamma, eta, 
           pop, 
           seed) {

  best_params_df <- best_parameters$best_params_df
  
  scores_list <- data.frame(matrix(ncol = 4, nrow = 3))
  colnames(scores_list) <- c('Set chosen','CRPS','NLPD','SPE')
  
  for (i in 1:nrow(best_params_df)) {
    # Extract hyperparameters
    params <- best_params_df[i, ]
    
    # Initialize with the current set of hyperparameters
    initial_params <- initialization(model = model, obs = obs_to_use,
                                     beta0 = beta0, beta0prime = params$beta0prime,
                                     lambda = lambda, gamma = gamma, eta = eta,
                                     l = params$l, scale = 1, noise_obs = noise_obs,
                                     noise_X = sqrt(noise_obs), noise_U = 0.2,
                                     noise_wiener_X = params$noise_wiener_X, 
                                     noise_wiener_U = params$noise_wiener_U,
                                     pop = pop,
                                     start_S = NULL, start_E = NULL, 
                                     start_I = NULL, start_R = NULL,
                                     start_D = NULL)
    
    # Generate time grids for inference
    grids <- generate_grid(obs_to_use, num_points_between = 0)
    
    # Run inference
    inference_results <- inference(model = model, 
                                   grids = grids, 
                                   obs = obs_to_use, 
                                   jit = params$jit, 
                                   initial_params = initial_params)
    
    # Process data for scoring
    processed_data <- process_data(inference_results = inference_results,
                                   grids = grids,
                                   model = model,
                                   pop = pop, 
                                   lambda = lambda, 
                                   gamma = gamma, 
                                   eta = eta)
    
    U_plot <- processed_data$U_plot

    # Compute and store scores
    scores <- compute_scores_and_table(U_plot = U_plot, df_beta = df_beta)[[2]]
    scores_list[i,1] <- params$Metric
    scores_list[i,2:4] <- scores$Value
  }
  
  # Create individual plots for CRPS, NLPD, and SPE
  crps_plot <- ggplot(scores_list, aes(x = `Set chosen`, y = CRPS, color = `Set chosen`)) +
    geom_point(size = 4) +
    labs(title = "CRPS scores for different hyperparameter sets",
         y = "CRPS", x = "Set Chosen") +
    theme_minimal() +
    scale_color_manual(values = c("Best mean CRPS" = "#E69F00", 
                                  "Best mean NLPD" = "#E69F00", 
                                  "Best mean SPE" = "#E69F00"))
  
  nlpd_plot <- ggplot(scores_list, aes(x = `Set chosen`, y = NLPD, color = `Set chosen`)) +
    geom_point(size = 4) +
    labs(title = "NLPD scores for different hyperparameter sets",
         y = "NLPD", x = "Set Chosen") +
    theme_minimal() +
    scale_color_manual(values = c("Best mean CRPS" = "#56B4E9", 
                                  "Best mean NLPD" = "#56B4E9", 
                                  "Best mean SPE" = "#56B4E9"))
  
  spe_plot <- ggplot(scores_list, aes(x = `Set chosen`, y = SPE, color = `Set chosen`)) +
    geom_point(size = 4) +
    labs(title = "SPE scores for different hyperparameter sets",
         y = "SPE", x = "Set Chosen") +
    theme_minimal() +
    scale_color_manual(values = c("Best mean CRPS" = "#009E73", 
                                  "Best mean NLPD" = "#009E73", 
                                  "Best mean SPE" = "#009E73"))
  
  # Combine the three plots into a list
  combined_plot <- list(crps_plot, nlpd_plot, spe_plot)
  
  best_scores_output <- list(best_scores_df = scores_list, 
                             best_scores_plot = combined_plot)
  
  return(best_scores_output)
}
