#################################### PROBODER ##################################
################################ Claire Descombes ##############################

################################
####### SIMULATED DATA #########
################################

library(HETTMO)
params_unstratified = set_parameters()
params_unstratified$p_detect1 <- 1
params_unstratified$p_detect2 <- 1
population_simulate <- params_unstratified$popsize
simulate <- simulate_data(params = params_unstratified, ts = 1:45)
I <- simulate[[1]]
S <- rep(0,length(I))
for (i in 1:length(I)){
  S[i] <- P - sum(I[1:i])
}
R <- simulate[[6]]
real_beta <- simulate[[3]]
date <- 1:length(S)
observations_simulate <- cbind(date,S,I,R)

################################
########## REAL DATA ###########
################################

# Total population (Bern, 2021).
population <-  1047473

# Import dataset.
setwd("~/Nextcloud/Documents/Mathe/HS23/Master thesis/data_covid_dashboard/sources-csv/data")

daily_or_weekly <- 'weekly' # choose either 'daily' or 'weekly'

if (daily_or_weekly == 'daily'){
  csv_data_cases <- read.csv(file = "COVID19Cases_geoRegion.csv")
  csv_data_death <- read.csv(file = "COVID19Death_geoRegion.csv")
} else {
  csv_data_cases <- read.csv(file = "COVID19Cases_geoRegion_w.csv")
  csv_data_death <- read.csv(file = "COVID19Death_geoRegion_w.csv")
}

#'library(rjson)
#'library(jsonlite)
#'setwd("~/Nextcloud/Documents/Mathe/HS23/Master thesis/data_covid_dashboard/sources-json/data")
#'json_data <- fromJSON(file = "COVID19Cases_geoRegion_w.json")
#'json_dataframe <- as.data.frame(json_data)

# Data selection.
library(dplyr)
region <- "BE"

data_cases <- csv_data_cases %>%
  filter(geoRegion == region)
data_death <- csv_data_death %>%
  filter(geoRegion == region)

cases <- data_cases %>%
  rename(date = datum) %>%
  group_by(date) %>%
  summarize(cases = sum(entries, na.rm = TRUE))
if (daily_or_weekly == 'daily'){
  cases <- mutate(cases, date = as.Date(date))
}

deaths <- data_death %>%
  rename(date = datum) %>%
  group_by(date) %>%
  summarize(deaths = sum(entries, na.rm = TRUE))
if (daily_or_weekly == 'daily'){
  deaths <- mutate(deaths, date = as.Date(date))
}

#' infections_14_days_ago <- bern_data_cases %>%
#' mutate(date = as.Date(datum)) %>%
#' mutate(infections_14_days_ago = lag(entries, 14, default = 0)) %>%
#' select(date, infections_14_days_ago) 

#' recovered_per_day <- infections_14_days_ago %>%
#' mutate(recovered_per_day = pmax(0, infections_14_days_ago - deaths$deaths)) %>%
#' select(date, recovered_per_day)

if (daily_or_weekly == 'daily'){
  dates <- seq(as.Date("2020-02-24"), as.Date("2023-01-01"), by = "day")
} else {
  dates <- cases$date
}

susceptibles_vector <- numeric(length(dates))
susceptibles_vector[1] <- population
for (i in 2:length(dates)) {
  cases_by_date <- filter(cases, date == as.character(dates[i]))
  if (nrow(cases_by_date) > 0) {
    cases_by_date <- cases_by_date$cases[1]
    susceptibles_vector[i] <- susceptibles_vector[i-1] - cases_by_date
  } else {
    susceptibles_vector[i] <- susceptibles_vector[i-1]
  }
}
susceptibles <- data.frame(date = dates, susceptibles = susceptibles_vector)

observations <- left_join(susceptibles, cases, by = "date") %>%
  #left_join(recovered_per_day, by = "date") %>%
  left_join(deaths, by = "date")
observations$deaths[is.na(observations$deaths)] <- 0

n <- nrow(observations) # size of data_grid

#' Train and validation sets.
#' set.seed(123)
#' indices <- seq_len(n)
#' train_indices <- sample(indices, size = 0.8 * n, replace = FALSE)
#' train_set <- observations[train_indices, ]
#' validation_set <- observations[-train_indices, ]

#####################################
############# FUNCTIONS #############
#####################################

# Link function.
sigmoid <- function(z){
  return(1/(1+exp(-z)))
}

# ODE.
f <- function(X0,beta,P,gamma,eta){
  # Arguments:
  #   arg1: X0, matrix(nrow=4, ncol=1), solution of the ODE
  #   arg2: beta, numeric with values in [0,1], latent parameter of the ODE
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
  
  S_out <- -beta*S*I/P
  I_out <- beta*S*I/P - gamma*I - eta*I
  R_out <- gamma*I
  D_out <- eta*I
  
  sol <- c(S_out,I_out,R_out,D_out)
  
  return(sol)
}

# Measurement model.
h <- function(X,U){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, matrix(nrow=2, ncol=1), latent parameter of the ODE and it's first derivative
  #
  # Returns:
  #   output: evaluation of the measurement model h, matrix(nrow=4, ncol=1)
  
  beta <- U[1]
  beta <- sigmoid(beta) # rescaling of beta to [0,1]
  
  X0 <- X[1:4]
  X1 <- X[5:8]
  ODE <- f(X0,beta)
  
  sol <- matrix(X1 - ODE, nrow = 4, ncol = 1)
  
  return(sol)
}

#' Jacobian of measurement model using jacobian function.
jacobian_measurement <- function(X,U){
  h1 <- function(X) h(X,U=U); h2 <- function(U) h(X=X,U)
  out <- cbind(jacobian(h1, x = X),jacobian(h2, x = U))
  return(out)
}

# Jacobian of f 'by hand'.
jacobian_U <- function(X,U,P,gamma,eta){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, numeric with values in [0,1], latent parameter of the ODE
  #   arg3: P, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: matrix(nrow=2, ncol=2), Jacobian of ODE f wrt U
  
  X0 <- X[1:4]
  
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  beta <- U
  
  d_dbeta <- matrix(
    data=c(-S * I / P, S * I / P, 0, 0),
    nrow=1, ncol=4
  )
  
  sol <- d_dbeta
  
  return(sol)
}

jacobian_X <- function(X,U,P,gamma,eta){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, numeric with values in [0,1], latent parameter of the ODE
  #   arg3: P, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: matrix(nrow=4, ncol=4), Jacobian of ODE f wrt X0
  
  X0 <- X[1:4]
  
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  beta <- U
  
  d_dS <- matrix(
    data=c(-beta * I / P, -beta * S / P, 0, 0),
    nrow=1, ncol=4
  )
  
  d_dI <- matrix(
    data=c(beta * I / P, beta * S / P - gamma - eta, 0, 0),
    nrow=1, ncol=4
  )
  
  d_dR <- matrix(
    data=c(0, gamma, 0, 0),
    nrow=1, ncol=4
  )
  
  d_dD <- matrix(
    data=c(0, eta, 0, 0),
    nrow=1, ncol=4
  )
  
  sol <- rbind(d_dS,d_dI,d_dR,d_dD)
  
  return(sol)
}
  
# Prediction step.
prediction_U <- function(m_U,P_U,F_U,L_U){
  # Arguments:
  #   arg1: m_U, numeric, (previous) mean of U
  #   arg2: P_U, numeric, (previous) (co)variance of U
  #   arg3: F_U, matrix(nrow=2, ncol=2), drif matrix of U
  #   arg4: L_U, matrix(nrow=2, ncol=1), dispersion matrix of U
  #
  # Returns:
  #   output: predicted mean and covariance U
  
  exp_F_U <- expm(F_U)
  m_U_out <- exp_F_U %*% m_U

  top_row_U <- cbind(F_U,L_U %*% t(L_U))
  bottom_row_U <- cbind(matrix(0,nrow=2,ncol=2),-t(F_U))
  Gamma_U <- rbind(top_row_U,bottom_row_U)
  exp_Gamma_U <- expm(Gamma_U)
  
  M_1_U <- exp_Gamma_U[1:2,1:2]
  M_2_U <- exp_Gamma_U[1:2,3:4]
  
  B_U <- M_2_U %*% t(M_1_U)
  
  P_U_out <- exp_F_U %*% P_U %*% t(exp_F_U) + B_U
  
  out <- list(m_U_out,P_U_out)
  
  return(out)
}

prediction_X <- function(m_X,P_X,F_X,L_X){
  # Arguments:
  #   arg1: m_X, vector(nrow=13), (previous) mean of X
  #   arg2: P_X, matrix(nrow=12, ncol=12), (previous) covariance of X
  #   arg3: F_X, matrix(nrow=12, ncol=12), drif matrix of X
  #   arg4: L_X, matrix(nrow=12, ncol=4), drif matrix of X
  #
  # Returns:
  #   output: predicted mean and covariance of X
  
  exp_F_X <- expm(F_X)
  m_X_out <- exp_F_X %*% m_X

  top_row_X <- cbind(F_X,L_X %*% t(L_X))
  bottom_row_X <- cbind(matrix(0,nrow=12,ncol=12),-t(F_X))
  Gamma_X <- rbind(top_row_X,bottom_row_X)
  exp_Gamma_X <- expm(Gamma_X)

  M_1_X <- exp_Gamma_X[1:12,1:12]
  M_2_X <- exp_Gamma_X[1:12,13:24]
  
  B_X <- M_2_X %*% t(M_1_X)

  P_X_out <- exp_F_X %*% P_X %*% t(exp_F_X) + B_X

  out <- list(m_X_out,P_X_out)
  
  return(out)
}

# Update step on tau_OBS.
update_of_observations <- function(m,P,y,H,R){
  # Arguments:
  #   arg1: m, vector(nrow=14), predicted mean of U and X
  #   arg2: P, matrix(nrow=14, ncol=14), predicted covariance of U and X
  #   arg3: y, vector(nrow=y), observations
  #   arg4: H, matrix(nrow=3, ncol=14), observation matrix
  #   arg5: R, matrix(nrow=3, ncol=3), observation noise
  #
  # Returns:
  #   output: vector(ncol=14) and matrix(nrow=14, ncol=14), updated
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
J <- jacobian_measurement(X,U)
h <- h(X,U)
update_of_states <- function(m,P,h,J){
  # Arguments:
  #   arg1: m, vector(nrow=14), predicted mean of U and X
  #   arg2: P, matrix(nrow=14, ncol=14), predicted covariance of U and X
  #   arg3, h, vector(nrow=4), measurement model
  #   arg4: J, matrix(nrow=4, ncol=14), Jacobian of the measurement model
  #
  # Returns:
  #   output: vector(ncol=14) and matrix(nrow=14, ncol=14), updated
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

matrix_P <- function(P_X,P_U){
  # Arguments:
  #   arg1: P_X, matrix(nrow=12, ncol=12), predicted covariance of X
  #   arg2: P_U, matrix(nrow=2, ncol=2), predicted covariance of U
  #
  # Returns:
  #   output: matrix(nrow=14, ncol=14), predicted covariance of X and U
  
  P_X <- as.matrix(P_X)
  P_U <- as.matrix(P_U)
  
  nrow_P <- nrow(P_X) + nrow(P_U)
  ncol_P <- ncol(P_X) + ncol(P_U)
  P <- matrix(0, nrow = nrow_P, ncol = ncol_P)
  P[1:nrow(P_X), 1:ncol(P_X)] <- P_X  # Top left corner
  P[(nrow_P - nrow(P_U) + 1):nrow_P, (ncol_P - ncol(P_U) + 1):ncol_P] <- P_U  # Bottom right corner
  
  return(P)
}

#####################################
########## INITIALIZATION ###########
#####################################

# Choice of data (in any case: date-S-I-D data).
type <- 'simulated' # set 'real' for real data, 'simulated' for simulated data
if (type == 'simulated'){
  obs <- observations_simulate
}else{
  obs <- observations
}

# X: solution of SIRD-ODE and its two first derivatives
X <- c(data= c(obs[1,2],rep(0,11)))
# U: latent parameter (contact rate) and its first derivative
U <- c(0,0)
# P: total population
if (type == 'simulated'){
  P <- population_simulate
}else{
  P <- population
}

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
F_X <- as.matrix(F_X)

# Dispersion matrices.
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4))
L_X <- as.matrix(L_X)

# Observation matrix (for observation of S,I and D).
H <- sparseMatrix(i = c(1,2,3), j = c(1,2,4), x = 1, dims = c(3,14))
H <- as.matrix(H)

# Observation noise.
R <- matrix(0.001, nrow = 3, ncol = 3)

# Noise of priors.
P_X <- matrix(0.001, nrow = 12, ncol = 12)
P_U <- matrix(0.001, nrow = 2, ncol = 2)

#####################################
############# ALGORITHM #############
#####################################

# Data grid.
data_grid <- obs[,'date']

# ODE grid.
ode_grid <- data_grid # more points could be added

# Overall time grid.
time_grid <- sort(unique(c(data_grid, ode_grid)))

# 'Artificial' observations for ODE measurements
zero_data = rep(x=0, 4)

data_idx <- 0
ode_idx <- 0

# Arrays to store values of X, P_X and U, P_U.
X_values <- matrix(data = NA, nrow = 12, ncol = length(time_grid))
U_values <- matrix(data = NA, nrow = 2, ncol = length(time_grid))
P_X_values <- array(data = NA, dim = c(12, 12, length(time_grid)))
P_U_values <- array(data = NA, dim = c(2, 2, length(time_grid)))

for (loc in time_grid){
  X_values[,loc] <- as.vector(X)
  U_values[,loc] <- as.vector(U)
  P_X_values[,,loc] <- as.matrix(P_X)
  P_U_values[,,loc] <- as.matrix(P_U)
  
  # Prediction step.
  U <- prediction_U(m_U=U,P_U=P_U,F_U=F_U,L_U=L_U)[[1]]
  P_U <- prediction_U(m_U=U,P_U=P_U,F_U=F_U,L_U=L_U)[[2]]
  X <- prediction_X(m_X=X,P_X=P_X,F_X=F_X,L_X=L_X)[[1]]
  P_X <- prediction_X(m_X=X,P_X=P_X,F_X=F_X,L_X=L_X)[[2]]
  
  # Update of observations.
  if (any(data_grid == loc)){
    m <- c(X,U)
    P <- matrix_P(P_X,P_U)
    y <- obs[which(obs[, 1] == loc),2:4]
    X <- update_of_observations(m,P,y,H,R)[[1]][1:4]
    U <- update_of_observations(m,P,y,H,R)[[1]][5:6]
    P_X <- update_of_observations(m,P,y,H,R)[[2]][1:12,1:12]
    P_U <- update_of_observations(m,P,y,H,R)[[2]][13:14,13:14]
  }
  
  # Update of states.
  if (any(ode_grid == loc)){
    P <- matrix_P(P_X,P_U)
    J <- jacobian_measurement(X,U,P,gamma,eta)
    m <- c(X,U)
    X <- update_of_states(m,P,h,J)[[1]][1:4]
    U <- update_of_states(m,P,h,J)[[1]][5:6]
    P_X <- update_of_states(m,P,h,J)[[2]][1:12,1:12]
    P_U <- update_of_states(m,P,h,J)[[2]][13:14,13:14]
  }
}
