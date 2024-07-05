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
#' @return A styled table displaying Mean Squared Prediction Error (SPE), Mean Negative Log Predictive Density (NLPD), and Mean Continuous Ranked Probability Score (CRPS).
#' @export
compute_scores_and_table <- function(U_plot, df_beta) {
  # Compute the different scores
  SPE <- mean(squared_prediction_error(U_plot, df_beta))
  NLPD <- mean(negative_log_predictive_density(U_plot, df_beta))
  CRPS <- mean(continuous_ranked_probability_score(U_plot, df_beta))
  
  # Create a data frame with the results
  results <- data.frame(
    Method = c("Mean Squared Prediction Error", 
               "Mean Negative Log Predictive Density", 
               "Mean Continuous Ranked Probability Score"),
    Value = c(SPE, NLPD, CRPS)
  )
  
  # Generate nice table using knitr and kableExtra
  table <- kable(results, align = "c", caption = "Scoring Methods Results")
  styled_table <- kableExtra::kable_styling(table, 
                                            bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                                            full_width = FALSE)
  
  return(styled_table)
}