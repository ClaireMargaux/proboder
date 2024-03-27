#####################################
############# WORKFLOW ##############
#####################################

# Import functions
source('~/Documents/GitHub/proboder/initialization.R')
source('~/Documents/GitHub/proboder/functions_for_inference.R')
source('~/Documents/GitHub/proboder/saving_loading_plotting.R')

# Necessary packages
library(Matrix) # for sparseMatrix()
library(numDeriv) # for jacobian()
library(matrixcalc) # for svd.inverse()
library(ggplot2) # for ggplot()

# Choose data to be imported (in case 'real': date-S-I-D, in case 'simulated': date-S-I-R)
directory_data <- "~/Documents/GitHub/proboder/Data" # directory of data
type <- 'simulated' # set 'real' for real data, 'simulated' for simulated data
region <- 'BE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'weekly' # choose either 'daily' or 'weekly' (if 'real' data selected)

# Set model type
model <- if(type == 'simulated'){"SIR"}else{"SID"}

# Import data
data <- load_data(type,region,daily_or_weekly,directory_data)
obs <- data$observations
pop <- data$population
if(type == 'simulated'){
  real_beta <- data$real_beta
}

# Sanity check.
head(obs)
summary(obs)

#####################################
########## INITIALIZATION ###########
#####################################

initial_params <- 
  initialization(model, obs, 
                 beta0 = 0.99, beta0prime = -2.5, 
                 gamma = 0.4, eta = 0, 
                 l = 8, pop)

#####################################
############# ALGORITHM #############
#####################################

# Data grid
data_grid <- obs[,'date']

# ODE grid
ode_grid <- data_grid # no more points than observations

# Overall time grid
time_grid <- sort(unique(c(data_grid, ode_grid)))

# Run inference.
inference_results <- inference(time_grid, obs, initial_params)

X_values <- inference_results$X_values
U_values <- inference_results$U_values
P_X_values <- inference_results$P_X_values
P_U_values <- inference_results$P_U_values

# ------------
# Save results
# ------------

# Specify directory for results
directory_res = "~/Documents/GitHub/proboder/Results"
# Save results to the specified directory
save_results_as_Rdata(X_values, U_values, P_X_values, P_U_values, directory_res)

#####################################
########### VISUALIZATION ###########
#####################################

# ---------------------
# Extract relevant data
# ---------------------

# Load and process data from the specified directory
processed_data <- load_and_process_data(directory_res,time_grid)
U_plot <- processed_data$U_plot
P_plot <- processed_data$P_plot
ymin <- P_plot$ymin
ymax <- P_plot$ymax
U_scaled <- processed_data$U_scaled
Xval <- processed_data$Xval

# Save processed data to the specified directory
save_processed_data(U_plot, P_plot, ymin, ymax, U_scaled, Xval, directory_res)

# Create data frame for real beta values (if available)
if(type=='simulated'){
  real_beta_df <- data.frame(time = data_grid, real_beta = real_beta)
}

# --------
# Plotting
# --------

setwd(directory_res)

# Plot data
pdf("SIR-counts.pdf")
plot_data(obs,Xval,model)
dev.off()
plot_data(obs,Xval,model)

# Plot contact rate
pdf("contact-rate-with-CI.pdf")
plot_contact_rate(type, U_plot, ymin, ymax, U_scaled, real_beta_df, gamma, eta, l)
dev.off()
plot_contact_rate(type, U_plot, ymin, ymax, U_scaled, real_beta_df, gamma, eta, l)