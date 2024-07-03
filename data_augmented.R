################################
########## REAL DATA ###########
################################

directory <- "~/home/claire/Documents/data_covid_dashboard/sources-csv/data" # directory of data
region <- 'GE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'daily' # choose either 'daily' or 'weekly' (if 'real' data selected)

# Assumptions
incubation_period <- 2  # Days from E to I
infectious_period <- 3  # Days from I to R
time_to_death <- 15 # Days from I to D

# Import dataset.
setwd("~/Documents/GitHub/proboder/data_covid_dashboard/sources-csv/data")

if (daily_or_weekly == 'daily'){
  csv_data_cases <- read.csv(file = "COVID19Cases_geoRegion.csv")
  csv_data_death <- read.csv(file = "COVID19Death_geoRegion.csv")
} else {
  csv_data_cases <- read.csv(file = "COVID19Cases_geoRegion_w.csv")
  csv_data_death <- read.csv(file = "COVID19Death_geoRegion_w.csv")
}

# Data selection.
library(dplyr)

# Total population (Bern or Geneva, 2021).
if(region == 'BE'){
  population <-  1047473
}else if(region == 'GE'){
  population <- 511921
}else{
  print('Invalid region!')
}

data_cases <- csv_data_cases %>%
  filter(geoRegion == region)
data_death <- csv_data_death %>%
  filter(geoRegion == region)

cases <- data_cases %>%
  rename(date = datum) %>%
  group_by(date) %>%
  summarize(cases = sum(entries, na.rm = TRUE)) %>%
  arrange(date) %>%  
  mutate(cases = cumsum(cases))
if (daily_or_weekly == 'daily'){
  cases <- mutate(cases, date = as.Date(date))
}

deaths <- data_death %>%
  rename(date = datum) %>%
  mutate(cumulative_entries = cumsum(entries)) %>%
  group_by(date) %>%
  summarize(deaths = sum(cumulative_entries, na.rm = TRUE))
if (daily_or_weekly == 'daily'){
  deaths <- mutate(deaths, date = as.Date(date))
}

if (daily_or_weekly == 'daily'){
  dates <- seq(as.Date("2020-02-24"), as.Date("2023-01-01"), by = "day")
} else {
  dates <- cases$date
}

observations <- inner_join(cases, deaths, by = "date") 
observations$deaths[is.na(observations$deaths)] <- 0
colnames(observations) <- c('t','R','D')
rm(csv_data_cases, csv_data_death, data_cases, data_death, cases, deaths)

# Initialize compartments
num_dates <- length(dates)
S <- numeric(num_dates)
ER <- numeric(num_dates)
ED <- numeric(num_dates)
E <- numeric(num_dates)
IR <- numeric(num_dates)
ID <- numeric(num_dates)
I <- numeric(num_dates)

# Backward computation of IR, ER
for (i in num_dates:2) {
  # Initialize new infections
  new_infections <- 0
  new_infections <- observations$R[i] - observations$R[i - 1]
  
  # Update compartment IR
  if (i - infectious_period >= 1) {
    for (j in 0:(infectious_period-1)) {
      IR[i - infectious_period + j] <- IR[i - infectious_period + j] + new_infections
    }
  } 
  
  # Update compartment ER
  if (i - infectious_period - incubation_period >= 1) {
    for (j in 0:(incubation_period-1)) {
      ER[i - infectious_period - incubation_period + j] <- ER[i - infectious_period - incubation_period + j] + new_infections
    }
  } 
}

# Backward computation of ID, ED
for (i in num_dates:2) {
  # Initialize new infections
  new_infections <- 0
  new_infections <- observations$D[i] - observations$D[i - 1]
  
  # Update compartment IR
  if (i - time_to_death >= 1) {
    for (j in 0:(time_to_death-1)) {
      ID[i - time_to_death + j] <- ID[i - time_to_death + j] + new_infections
    }
  } 
  
  # Update compartment ER
  if (i - time_to_death - incubation_period >= 1) {
    for (j in 0:(incubation_period-1)) {
      ED[i - time_to_death - incubation_period + j] <- ED[i - time_to_death - incubation_period + j] + new_infections
    }
  } 
}

# Computation of I (merging ID and IR)
I <- IR + ID
E <- ER + ED

# Remove temp variables
rm(IR,ID,ER,ED,i,j)

# Forward computation of S
S[1] = population
for (i in 2:num_dates) {
  S[i] <- S[i-1] - max(0,(E[i] - E[i-1]))
}

# Create dataframe for S, E, I, R, D compartments
compartments <- data.frame(
  date = observations$t,
  S = S,
  E = E,
  I = I,
  R = observations$R,
  D = observations$D
)
colnames(compartments) <- c('t','S','E','I','R','D')

compartments_to_plot <- compartments %>%
  mutate(
    S = ifelse(S <= 0, 0.001, S),
    E = ifelse(E <= 0, 0.001, E),
    I = ifelse(I <= 0, 0.001, I),
    R = ifelse(R <= 0, 0.001, R),
    D = ifelse(D <= 0, 0.001, D)
  )

countsplot <- ggplot(compartments_to_plot, aes(x = t, y = S)) +
  geom_line(aes(color = 'S')) +
  geom_line(aes(x = t, y = E, color = 'E')) +
  geom_line(aes(x = t, y = I, color = 'I')) +
  geom_line(aes(x = t, y = R, color = 'R')) +
  geom_line(aes(x = t, y = D, color = 'D')) +
  theme_minimal() +
  scale_y_continuous(trans='log10') + 
  ylab("log(count)") +
  xlab("time") + 
  ggtitle(paste(daily_or_weekly,"data")) 

countsplot

n <- nrow(observations) # size of data_grid

################################
######## SAVE THE DATA #########
################################

directory_save <- "~/Documents/GitHub/proboder/Data/real" # directory for saving

region_filename <- paste0("real_data_augmented_", region, "_", daily_or_weekly, ".Rdata")
population_filename <- paste0("real_pop_augmented_", region, "_", daily_or_weekly, ".Rds")
save(compartments, file = file.path(directory_save, region_filename))
saveRDS(population, file = file.path(directory_save, population_filename))