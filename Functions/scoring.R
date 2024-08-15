#' Calculate Squared Prediction Error (SPE)
#'
#' This function calculates the squared prediction error (SPE) between the inferred and real values.
#' SPE is defined as the square of the difference between the inferred value and the real value.
#'
#' @param U_plot Data frame containing the inferred contact rate values.
#' @param real_beta_df Data frame containing the real contact rate values.
#' @return The mean squared prediction error.
#' @export
squared_prediction_error <- function(U_plot, real_beta_df){
  SPE <- ((U_plot$U_scaled - real_beta_df$beta)^2)
  return(SPE)
}

#' Calculate Negative Log Predictive Density (NLPD)
#'
#' This function calculates the negative log predictive density (NLPD) between the inferred and real values.
#' NLPD is defined as the negative logarithm of the probability density function of the real value given the inferred value.
#'
#' @param U_plot Data frame containing the inferred contact rate values.
#' @param real_beta_df Data frame containing the real contact rate values.
#' @return The mean negative log predictive density.
#' @export
negative_log_predictive_density <- function(U_plot, real_beta_df){
  NLPD <- (-log(dnorm(x = real_beta_df$beta, mean = U_plot$U_scaled, sd = sqrt(U_plot$P_U_scaled))))
  return(NLPD)
}

#' Calculate Continuous Ranked Probability Score (CRPS)
#'
#' This function calculates the continuous ranked probability score (CRPS) between the inferred and real values.
#' CRPS is a proper score function measuring the accuracy of probabilistic forecasts.
#'
#' @param U_plot Data frame containing the inferred contact rate values.
#' @param real_beta_df Data frame containing the real contact rate values.
#' @return The mean continuous ranked probability score.
#' @export
continuous_ranked_probability_score <- function(U_plot, real_beta_df){
  mu <- U_plot$U_scaled
  sigma <- sqrt(U_plot$P_U_scaled)
  y <- real_beta_df$beta
  omega <- (y - mu) / sigma
  
  CRPS <- (sigma * (
    omega * (2 * pnorm(omega) - 1) 
    + 2 * dnorm(omega) 
    - 1/sqrt(pi)
  ))
  
  return(CRPS)
}

#' Compute Scores and Generate Styled Table
#' 
#' This function computes various scoring metrics based on the input data and returns a styled table using kableExtra.
#' 
#' @param U_plot Data frame or matrix containing the predicted values.
#' @param df_beta Data frame or matrix containing the true values for comparison.
#' @return A styled table displaying mean and median
#' - Squared Prediction Error (SPE),
#' - Continuous Ranked Probability Score (CRPS), and
#' - Negative Log Predictive Density (NLPD)
#' @export
compute_scores_and_table <- function(U_plot, df_beta) {
  # Compute the different scores
  mean_SPE <- mean(squared_prediction_error(U_plot, df_beta))
  mean_CRPS <- mean(continuous_ranked_probability_score(U_plot, df_beta))
  mean_NLPD <- mean(negative_log_predictive_density(U_plot, df_beta))
  
  med_SPE <- median(squared_prediction_error(U_plot, df_beta))
  med_CRPS <- median(continuous_ranked_probability_score(U_plot, df_beta))
  med_NLPD <- median(negative_log_predictive_density(U_plot, df_beta))
  
  # Create a data frame with the results
  results <- data.frame(
    Metric = c("Mean Squared Prediction Error", 
               "Median Squared Prediction Error",
               "Mean Continuous Ranked Probability Score",
               "Median Continuous Ranked Probability Score",
               "Mean Negative Log Predictive Density",
               "Median Negative Log Predictive Density"),
    Value = c(mean_SPE, med_SPE, 
              mean_CRPS, med_CRPS, 
              mean_NLPD, med_NLPD)
  )
  
  # Generate nice table using knitr and kableExtra
  table <- kable(results, align = "l", caption = "Scoring results")
  scores_table <- kableExtra::kable_styling(table, 
                                            bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                                            full_width = FALSE)
  
  return(scores_table)
}
