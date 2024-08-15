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
#' @param clip_X Clip on the variance of X.
#' @param clip_U Clip on the variance of U.
#' @param jit_X Jitter on the variance of X.
#' @param jit_U Jitter on the variance of U.
#' @param seed Seed for random number generation. Default is 5.
#' @param num_param_sets Number of parameter sets to sample.
#' @param seq_l Sequence of values for the `l` hyperparameter. 
#' @param seq_wiener_X Sequence of values for the `noise_wiener_X` hyperparameter.
#' @param seq_wiener_U Sequence of values for the `noise_wiener_U` hyperparameter.
#' @param seq_beta0prime Sequence of values for the `beta0prime` hyperparameter. 
#' @param seq_jit Sequence of values for the `jit` hyperparameter. 
#'
#' @return A data frame containing the summary of results for each sampled set of hyperparameters. The data frame includes the hyperparameter values, mean and median scores for SPE, NLPD, and CRPS.
#' @export
run_random_search <- function(model, obs_to_use, df_beta, noise_obs, 
                              lambda, gamma, eta, pop, beta0,
                              clip_X, clip_U,
                              jit_X, jit_U,
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
      noise_U = 0.01,
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
                jit_X = jit_X,
                jit_U = jit_U,
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
      score_SPE <- squared_prediction_error(U_plot, df_beta)
      score_NLPD <- negative_log_predictive_density(U_plot, df_beta)
      score_CRPS <- continuous_ranked_probability_score(U_plot, df_beta)
      
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
    jit = params[, 5],
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
#' (SPE, NLPD, CRPS) across different hyperparameters. The function can optionally aggregate scores
#' for each hyperparameter value or plot raw scores with NA values removed.
#'
#' @param results_summary A data frame containing the summary of results from the hyperparameter search.
#' Each row represents a different set of hyperparameters and contains columns for the mean and median
#' scores (`mean_SPE`, `mean_NLPD`, `mean_CRPS`, `median_SPE`, `median_NLPD`, `median_CRPS`).
#' @param parameters_to_plot A character vector of parameter names for which the mean and median scores
#' should be plotted. Default is `c('l', 'noise_wiener_X', 'noise_wiener_U', 'beta0prime', 'jit')`.
#' @param aggregate A logical value indicating whether to aggregate scores by the parameter values.
#' If `TRUE`, scores are aggregated by taking the mean for each parameter value. If `FALSE`, 
#' raw scores with NA values removed are plotted. Default is `FALSE`.
#' @return A list containing two lists: `mean_plots` and `median_plots`. Each list contains ggplot objects
#' for the specified parameters. The `mean_plots` list contains plots of mean scores, and the `median_plots`
#' list contains plots of median scores. If `aggregate` is `TRUE`, the plots reflect aggregated scores; 
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
  color_scale <- scale_color_manual(values = c("SPE" = "#E69F00", "NLPD" = "#56B4E9", "CRPS" = "#009E73"), 
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
  
  aggregate_data_median <- function(df, parameter) {
    df %>%
      group_by(across(all_of(parameter))) %>%  # Group by the parameter
      summarize(
        median_SPE = median(median_SPE, na.rm = TRUE),
        median_NLPD = median(median_NLPD, na.rm = TRUE),
        median_CRPS = median(median_CRPS, na.rm = TRUE),
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
    
    ggplot(aggregated_df, aes_string(x = parameter)) +
      geom_point(aes(y = mean_SPE, color = "SPE")) +
      geom_point(aes(y = mean_NLPD, color = "NLPD")) +
      geom_point(aes(y = mean_CRPS, color = "CRPS")) +
      labs(title = title_text, 
           x = TeX(params_expressions[parameter]), 
           y = "Mean score",
           color = "Score") +
      color_scale +  # Use predefined color scale
      theme_minimal()
  }
  
  # Function to plot median scores for a given parameter
  plot_median_scores <- function(df, parameter, aggregate) {
    # Aggregate data if needed
    if (aggregate) {
      aggregated_df <- aggregate_data_median(df, parameter)
      title_text <- TeX(paste("Aggregated median scores for", params_expressions[parameter]))
    } else {
      aggregated_df <- df %>%
        filter(!is.na(median_SPE) & !is.na(median_NLPD) & !is.na(median_CRPS))
      title_text <- TeX(paste("Median scores for", params_expressions[parameter]))
    }
    
    ggplot(aggregated_df, aes_string(x = parameter)) +
      geom_point(aes(y = median_SPE, color = "SPE")) +
      geom_point(aes(y = median_NLPD, color = "NLPD")) +
      geom_point(aes(y = median_CRPS, color = "CRPS")) +
      labs(title = title_text, 
           x = TeX(params_expressions[parameter]), 
           y = "Median score",
           color = "Score") +
      color_scale +  # Use predefined color scale
      theme_minimal()
  }
  
  # Create a list to store the plots
  plot_list_mean <- list()
  plot_list_median <- list()
  
  # Generate and store plots for each parameter
  for (param in parameters_to_plot) {
    plot_list_mean[[param]] <- plot_mean_scores(results_summary, param, aggregate)
    plot_list_median[[param]] <- plot_median_scores(results_summary, param, aggregate)
  }
  
  # Return the lists of plots
  return(list(mean_plots = plot_list_mean, median_plots = plot_list_median))
}

#' Find the best set of hyperparameters
#'
#' This function identifies the best set of hyperparameters from a summary of results based on mean and median scores for SPE, NLPD, and CRPS.
#' It first checks for a set that minimizes all three scores (mean and median), then checks for sets minimizing at least two scores (mean and median),
#' and finally returns the best sets with respect to individual scores (mean and median) and overall.
#' If no sets meet the criteria for minimizing two or all three scores, `NA` values are included in the results.
#'
#' @param results_summary A data frame containing the summary of results with mean and median scores (SPE, NLPD, CRPS) for each set of hyperparameters.
#'
#' @return A data frame containing the best sets of hyperparameters based on the criteria outlined, including:
#' \item{Best SPE (Mean)}{The hyperparameters that yield the minimum mean SPE.}
#' \item{Best NLPD (Mean)}{The hyperparameters that yield the minimum mean NLPD.}
#' \item{Best CRPS (Mean)}{The hyperparameters that yield the minimum mean CRPS.}
#' \item{Best All (Mean)}{The hyperparameters that minimize all three or at least two mean scores, or `NA` if no such set exists.}
#' \item{Best SPE (Median)}{The hyperparameters that yield the minimum median SPE.}
#' \item{Best NLPD (Median)}{The hyperparameters that yield the minimum median NLPD.}
#' \item{Best CRPS (Median)}{The hyperparameters that yield the minimum median CRPS.}
#' \item{Best All (Median)}{The hyperparameters that minimize all three or at least two median scores, or `NA` if no such set exists.}
#'
#' @export
find_best_hyperparameters <- function(results_summary) {
  
  # Find the best set of hyperparameters for each score (mean)
  best_spe_mean <- results_summary[which.min(results_summary$mean_SPE), ]
  best_nlpd_mean <- results_summary[which.min(results_summary$mean_NLPD), ]
  best_crps_mean <- results_summary[which.min(results_summary$mean_CRPS), ]
  
  # Check if there is a set minimizing all three mean scores
  best_all_mean <- results_summary %>%
    filter(mean_SPE == min(mean_SPE) & mean_NLPD == min(mean_NLPD) & mean_CRPS == min(mean_CRPS))
  
  if (nrow(best_all_mean) > 0) {
    best_params_mean <- best_all_mean[1, ]
  } else {
    # Find sets minimizing at least two of the mean scores
    best_two_mean <- results_summary %>%
      filter((mean_SPE == min(mean_SPE) & mean_NLPD == min(mean_NLPD)) |
               (mean_SPE == min(mean_SPE) & mean_CRPS == min(mean_CRPS)) |
               (mean_NLPD == min(mean_NLPD) & mean_CRPS == min(mean_CRPS)))
    
    if (nrow(best_two_mean) > 0) {
      best_params_mean <- best_two_mean[1, ]
    } else {
      # No best_two found, set to NA
      best_params_mean <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(results_summary)))
      names(best_params_mean) <- names(results_summary)
      best_params_mean$Type <- "Best All (Mean)"
    }
  }
  
  # Find the best set of hyperparameters for each score (median)
  best_spe_median <- results_summary[which.min(results_summary$median_SPE), ]
  best_nlpd_median <- results_summary[which.min(results_summary$median_NLPD), ]
  best_crps_median <- results_summary[which.min(results_summary$median_CRPS), ]
  
  # Check if there is a set minimizing all three median scores
  best_all_median <- results_summary %>%
    filter(median_SPE == min(median_SPE) & median_NLPD == min(median_NLPD) & median_CRPS == min(median_CRPS))
  
  if (nrow(best_all_median) > 0) {
    best_params_median <- best_all_median[1, ]
  } else {
    # Find sets minimizing at least two of the median scores
    best_two_median <- results_summary %>%
      filter((median_SPE == min(median_SPE) & median_NLPD == min(median_NLPD)) |
               (median_SPE == min(median_SPE) & median_CRPS == min(median_CRPS)) |
               (median_NLPD == min(median_NLPD) & median_CRPS == min(median_CRPS)))
    
    if (nrow(best_two_median) > 0) {
      best_params_median <- best_two_median[1, ]
    } else {
      # No best_two found, set to NA
      best_params_median <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(results_summary)))
      names(best_params_median) <- names(results_summary)
      best_params_median$Type <- "Best All (Median)"
    }
  }
  
  # Extract only parameter columns from results_summary
  param_cols <- c("l", "noise_wiener_X", "noise_wiener_U", "beta0prime", "jit")
  
  # Function to create parameter rows
  create_param_row <- function(type_label, params) {
    params$Type <- type_label
    params <- params %>%
      select(Type, all_of(param_cols))  # Select only parameter columns
    return(params)
  }
  
  # Combine results into a single data frame with specific order
  combined_results <- rbind(
    create_param_row("Best SPE (Mean)", best_spe_mean),
    create_param_row("Best NLPD (Mean)", best_nlpd_mean),
    create_param_row("Best CRPS (Mean)", best_crps_mean),
    create_param_row("Best All (Mean)", best_params_mean),
    create_param_row("Best SPE (Median)", best_spe_median),
    create_param_row("Best NLPD (Median)", best_nlpd_median),
    create_param_row("Best CRPS (Median)", best_crps_median),
    create_param_row("Best All (Median)", best_params_median)
  )
  
  # Reorder columns to have 'Type' as the first column
  combined_results <- combined_results %>%
    select(Type, everything())
  
  # Generate table using knitr and kableExtra
  best_params_table <- kable(combined_results, align = "c", 
                             caption = "Best Hyperparameters for Each Score (Mean and Median)") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) %>%
    column_spec(1, bold = TRUE) 
  
  return(best_params_table)
}
