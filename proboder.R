#################################### PROBODER ##################################
################################ Claire Descombes ##############################

################################
############# DATA #############
################################

# Total population (Bern, 2021).
population <-  1047473

# Import dataset.
setwd("~/Nextcloud/Documents/Mathe/HS23/Master thesis/data_covid_dashboard/sources-csv/data")
csv_data_cases <- read.csv(file = "COVID19Cases_geoRegion.csv")
csv_data_death <- read.csv(file = "COVID19Death_geoRegion.csv")

#'library(rjson)
#'library(jsonlite)
#'setwd("~/Nextcloud/Documents/Mathe/HS23/Master thesis/data_covid_dashboard/sources-json/data")
#'json_data <- fromJSON(file = "COVID19Cases_geoRegion_w.json")
#'json_dataframe <- as.data.frame(json_data)

# Data selection.
library(dplyr)
bern_data_cases <- csv_data_cases %>%
  filter(geoRegion == "BE")
bern_data_death <- csv_data_death %>%
  filter(geoRegion == "BE")

cases_per_day <- bern_data_cases %>%
  rename(date = datum) %>%
  group_by(date) %>%
  summarize(cases_per_day = sum(entries, na.rm = TRUE))
cases_per_day <- mutate(cases_per_day, date = as.Date(date))

death_per_day <- bern_data_death %>%
  rename(date = datum) %>%
  group_by(date) %>%
  summarize(death_per_day = sum(entries, na.rm = TRUE))
death_per_day <- mutate(death_per_day, date = as.Date(date))

#'infections_14_days_ago <- bern_data_cases %>%
#'mutate(date = as.Date(datum)) %>%
#'mutate(infections_14_days_ago = lag(entries, 14, default = 0)) %>%
#'select(date, infections_14_days_ago) 

#'recovered_per_day <- infections_14_days_ago %>%
#'mutate(recovered_per_day = pmax(0, infections_14_days_ago - death_per_day$death_per_day)) %>%
#'select(date, recovered_per_day)

dates <- seq(as.Date("2020-02-24"), as.Date("2023-01-01"), by = "day")
susceptibles_vector <- numeric(length(dates))
susceptibles_vector[1] <- population
for (i in 2:length(dates)) {
  daily_cases <- filter(cases_per_day, date == as.character(dates[i]))
  if (nrow(daily_cases) > 0) {
    daily_cases <- daily_cases$cases_per_day[1]
    susceptibles_vector[i] <- susceptibles_vector[i-1] - daily_cases
  } else {
    susceptibles_vector[i] <- susceptibles_vector[i-1]
  }
}
susceptibles_per_day <- data.frame(date = dates, susceptibles_per_day = susceptibles_vector)

observations <- left_join(susceptibles_per_day, cases_per_day, by = "date") %>%
  #left_join(recovered_per_day, by = "date") %>%
  left_join(death_per_day, by = "date")

#####################################
########## INITIALIZATION ###########
#####################################

X <- matrix(data= c(as.numeric(observations[1,2]),rep(0,11)),nrow = 12, ncol = 1)
U <- 0
P <- population

library(Matrix)
library(numDeriv)
library(matrixcalc)

#' @param recovery_rate
#' @param fatality_rate
#' @param length_scale
#' @return contact_rate

# Fixed parameters.
gamma <- 0.06 # recovery_rate
eta <- 0.002 # fatality_rate
l <- 14 # length_scale

# Drift matrices.
F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
F_X <- sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12))
  
# Dispersion matrices.
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4))

# Observation matrix.
H <- sparseMatrix(i = c(1,2,3), j = c(1,2,4), x = 1, dims = c(3,4))

# Observation noise.
R <- diag(1, nrow = 3, ncol = 3)

# Link function.
sigmoid <- function(z){
  return(1/(1+exp(-z)))
}

# ODE.
f <- function(X0,U){
  # Arguments:
  #   arg1: X, matrix(nrow=4, ncol=1), solution of the ODE
  #   arg2: U, numeric with values in [0,1], latent parameter of the ODE
  #   arg3: P, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: sol, numeric, evaluation of the ODE f, vector of size 4
  
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  beta <- U
  
  sol <- c(
    -beta*S*I/P,
    beta*S*I/P - gamma*I - eta*I,
    gamma*I,
    eta*I
  )
  
  return(sol)
}

# Measurement model.
h <- function(X,U){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, numeric, latent parameter of the ODE
  #
  # Returns:
  #   output: evaluation of the measurement model h, matrix(nrow=4, ncol=1)
  
  U <- sigmoid(U) # rescaling of U to [0,1]
  
  derivative <- X[5:8]
  ODE <- f(X[1:4],U)
  
  sol <- matrix(derivative - ODE, nrow = 4, ncol = 1)
  
  return(sol)
}

# Jacobian.
jacobian_matrix <- function(X,U){
  h1 <- function(X) h(X,U=U); h2 <- function(U) h(X=X,U)
  out <- cbind(jacobian(h1, x = X),jacobian(h2, x = U))
  return(out)
}
  
# Prediction step.
prediction <- function(m_X,P_X,m_U,P_U,F_X,F_U,L_X,L_U){
  # Arguments:
  #   arg1: m_X, vector(nrow=13), (previous) mean of X
  #   arg2: P_X, matrix(nrow=12, ncol=12), (previous) covariance of X
  #   arg3: m_U, numeric, (previous) mean of U
  #   arg4: P_U, numeric, (previous) (co)variance of U
  #   arg5: F_X, matrix(nrow=12, ncol=12), drif matrix of X
  #   arg6: F_U, matrix(nrow=2, ncol=2), drif matrix of U
  #   arg7: L_X, matrix(nrow=12, ncol=4), dispersion matrix of X
  #   arg8: L_U, matrix(nrow=2, ncol=1), dispersion matrix of U
  #
  # Returns:
  #   output: predicted means and covariances of X and U
  
  exp_F_X <- expm(F_X)
  exp_F_U <- expm(F_U)
  
  m_X_out <- exp_F_X %*% m_X
  m_U_out <- exp_F_U %*% m_U
  
  top_row_X <- cbind(F_X,L_X %*% t(L_X))
  bottom_row_X <- cbind(matrix(0,nrow=12,ncol=12),-t(F_X))
  
  top_row_U <- cbind(F_U,L_U %*% t(L_U))
  bottom_row_U <- cbind(matrix(0,nrow=2,ncol=2),-t(F_U))

  Gamma_X <- rbind(top_row_X,bottom_row_X)
  Gamma_U <- rbind(top_row_U,bottom_row_U)
  
  exp_Gamma_X <- expm(Gamma_X)
  exp_Gamma_U <- expm(Gamma_U)
  
  M_1_X <- exp_Gamma_X[1:12,1:12]
  M_2_X <- exp_Gamma_X[1:12,13:24]
  
  M_1_U <- exp_Gamma_U[1:2,1:2]
  M_2_U <- exp_Gamma_U[1:2,3:4]

  B_X <- M_2_X %*% t(M_1_X)
  B_U <- M_2_U %*% t(M_1_U)
  
  P_X_out <- exp_F_X %*% P_X %*% t(exp_F_X) + B_X
  P_U_out <- exp_F_U %*% P_U %*% t(exp_F_U) + B_U
  
  out <- list(m_X_out,P_X_out,m_U_out,P_U_out)
  
  return(out)
}

# Update step on tau_OBS.
update_of_observations <- function(m,P,y,H,R){
  # Arguments:
  #   arg1: m, vector(nrow=13), predicted mean of U and X
  #   arg2: P, matrix(nrow=13, ncol=13), predicted covariance of U and X
  #   arg3: y, vector(nrow=y), observations
  #   arg4: H, matrix(nrow=3, ncol=4), observation matrix
  #   arg5: R, matrix(nrow=4, ncol=4), observation noise
  #
  # Returns:
  #   output: vector(ncol=13) and matrix(nrow=13, ncol=13), updated
  #           mean and covariances of U and X
  
  v <- y - H %*% m # residual
  S <- H %*% P %*% t(P) + R # innovation covariance
  S_inv <- svd.inverse(S)
  K <- P %*% t(H) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  out <- list(m_out,P_out)
  
  return(out)
}

# Update step on tau_ODE.
J <- jacobian_matrix(X,U)
h <- h(X,U)
update_of_states <- function(m,P,h,J){
  # Arguments:
  #   arg1: m, vector(nrow=13), predicted mean of U and X
  #   arg2: P, matrix(nrow=13, ncol=13), predicted covariance of U and X
  #   arg3, h, vector(nrow=4), measurement model
  #   arg4: J, matrix(nrow=4, ncol=13), Jacobian of the measurement model
  #
  # Returns:
  #   output: vector(ncol=13) and matrix(nrow=13, ncol=13), updated
  #           mean and covariances of U and X
  
  v <- -h # residual
  S <- J %*% P %*% t(J) # innovation covariance
  S_inv <- svd.inverse(S)
  K <- P %*% t(J) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  out <- list(m_out,P_out)
  
  return(out)
}

#####################################
############# ALGORITHM #############
#####################################

