#####################################
############# WORKFLOW ##############
#####################################

# Import functions
source('~/Documents/GitHub/proboder/LSODA.R')
source('~/Documents/GitHub/proboder/initialization.R')
source('~/Documents/GitHub/proboder/inference.R')
source('~/Documents/GitHub/proboder/processing_saving_loading.R')
source('~/Documents/GitHub/proboder/scoring.R')
source('~/Documents/GitHub/proboder/plotting.R')

# Necessary packages
library(Matrix) # for sparseMatrix() and expm()
library(numDeriv) # for jacobian()
library(matrixcalc) # for svd.inverse()
library(knitr) # for nice tables
library(kableExtra) # for storing nice tables
library(ggplot2) # for ggplot()
library(gridExtra) # for multiple plots
library(tictoc) # for benchmarking
library(progress) # to track progress of grid search
library(dplyr) # for treating data
library(lubridate) # to deal with dates

# Start counting computation time
tic("Duration of the whole workflow")

#################################
############# DATA ##############
#################################

# Model ('SEIRD' or 'SEIR')
model <- 'SEIR'

# Time grid
steps <- 1
max_time <- 30

# Fixed parameters
lambda <- 1/(2.6) # latency, ref value 1/(2.6)
gamma <- 1/(2.6) # recovery, ref value 1/(2.6)
eta <- 0 # fatality, ref value 0.024*(1/15), set to 0 if model = 'SEIR'

# Beta function
beta <- function(t){2*cos(t/30)-1}

# Starting compartment counts
if (model == 'SEIRD') {
  xstart <- c(S = 499980, E = 5, I = 5, R = 5, D = 0)
} else if (model == 'SEIR') {
  xstart <- c(S = 499985, E = 5, I = 5, R = 5)
}
pop <- sum(xstart)

# If noisy observation wished
noise <- 1 # noise to add on the data, set to 0 for data without noise
seed <- 5 # any integer, for reproducibility (set NA for no seed)

# Data simulation
sim <- simulate_data_LSODA(
  model = model,
  noise = noise,
  seed = seed,
  steps = steps,
  max_time = max_time,
  lambda = lambda, gamma = gamma, eta = eta,
  pop = pop, 
  beta = beta,
  xstart = xstart)

# Get simulated data sets
obs <- sim$obs
df_beta <- sim$df_beta
if (noise > 0) {
  obs_with_noise <- sim$obs_with_noise
}

# Sanity check.
head(obs)

# Visualization
plots <- plotting_simulated_data_lsoda(
  model = model,
  sim = sim,
  latency_rate = lambda,
  recovery_rate = gamma,
  fatality_rate = eta,
  log = TRUE)

simulated_compartments <- plots$simulated_compartments
simulated_beta <- plots$simulated_beta

print(simulated_compartments)
print(simulated_beta)

#####################################
########## INITIALIZATION ###########
#####################################

l <- 9.7 # lengthscale

initial_params <-
  initialization(model = model,
                 obs = obs, # choose if obs or obs_with_noise
                 beta0 = df_beta$beta[1], 
                 beta0prime = 0,
                 lambda = lambda, 
                 gamma = gamma, 
                 eta = eta,
                 l = l, 
                 scale = 1, 
                 noise_obs = noise,
                 noise_X = sqrt(noise), 
                 noise_U = 0.01,
                 noise_wiener_X = 50, 
                 noise_wiener_U = 0.01,
                 pop = pop)

# Generate time grids for inference
grids <- generate_grid(obs, num_points_between = 0)

#####################################
############# INFERENCE #############
#####################################

# Run inference
inference_results <- inference(model = model, 
                               grids = grids, 
                               obs = obs, 
                               initial_params = initial_params)

# Process data for visualization
processed_data <- process_data(inference_results,grids)
U_plot <- processed_data$U_plot
X_plot <- processed_data$X_plot

###############################
########### SCORING ###########
###############################

# Compute the different scores
SPE <- mean(squared_prediction_error(U_plot,df_beta))
NLPD <- mean(negative_log_predictive_density(U_plot,df_beta))
CRPS <- mean(continuous_ranked_probability_score(U_plot,df_beta))

# Create a data frame with the results
results <- data.frame(
  Method = c("Mean Squared Prediction Error", "Mean Negative Log Predictive Density", "Mean Continuous Ranked Probability Score"),
  Value = c(SPE, NLPD, CRPS)
)

# Generate nice table using knitr and kableExtra
table <- kable(results, align = "c", caption = "Scoring Methods Results")
(styled_table <- kableExtra::kable_styling(table, bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE))

#####################################
########### VISUALIZATION ###########
#####################################

# --------
# Plotting
# --------

# Plot compartment counts inferred from simulated data
compartments <- plot_compartments(model = model,
                                  obs = obs, 
                                  obs_with_noise = obs_with_noise, # if available, otherwise set NULL
                                  X_plot = X_plot)

print(compartments)

# Plot compartment counts separately
plots_sep <- plot_compartments_separately(model = model,
                                          obs = obs,
                                          obs_with_noise = NULL, # if available, otherwise set NULL
                                          X_plot = X_plot)
for (i in 1:nchar(model)) {
  plot_name <- paste0('plot_', i)
  assign(plot_name, plots_sep[[i]])
  print(get(plot_name))
}

# Plot contact rate with 95% confidence interval
contact_rate_with_CI <- plot_contact_rate_with_CI(U_plot = U_plot, 
                                                  df_beta = df_beta, # if available, otherwise set NULL
                                                  latency_rate = lambda, 
                                                  recovery_rate = gamma, 
                                                  fatality_rate = eta, 
                                                  lengthscale = l)

print(contact_rate_with_CI)

# ------
# Saving
# ------

if (FALSE) {
  
  directory_save <- '~/Documents/GitHub/proboder/Results/'
  save_processed_data(# Data
    inference_results = inference_results, 
    processed_data = processed_data, 
    directory = directory_save,
    simulated_compartments = simulated_compartments, 
    simulated_beta = simulated_beta,
    compartments = compartments,
    contact_rate_with_CI = contact_rate_with_CI,
    plots_sep = plots_sep,
    styled_table = styled_table)
  
}

# To store the table as a png file, use the manual export option from the viewer.
# Save in size 400x200.

# ------------
# Benchmarking
# ------------

# Total computation time
toc()
