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