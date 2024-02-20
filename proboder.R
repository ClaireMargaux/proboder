#################################### PROBODER ##################################
################################ Claire Descombes ##############################

################################
############# DATA #############
################################

# Total population (2021).
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

infections_14_days_ago <- bern_data_cases %>%
  mutate(date = as.Date(datum)) %>%
  mutate(infections_14_days_ago = lag(entries, 14, default = 0)) %>%  # Assuming 'entries' represents new cases per day
  select(date, infections_14_days_ago) 

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
############# ALGORITHM #############
#####################################

# Initialization.
X <- matrix(data= c(as.numeric(observations[1,2]),rep(0,11)),nrow = 12, ncol = 1)
U <- 0

library(Matrix)

#' @param recovery_rate
#' @param fatality_rate
#' @param lengthscale
#' @return 
#' 

# Fixed parameters.
recovery_rate <- 0.06
fatality_rate <- 0.002
lengthscale <- 14

# Drift matrices.
F_U <- matrix(c(0,-(sqrt(3)/lengthscale)^2,1,-2*sqrt(3)/lengthscale), nrow = 2, ncol = 2)
F_X <- sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12))
  
# Dispersion matrices.
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4))

# Observation noise.
R <- diag(1, nrow = 3, ncol = 3)

# Link function.
sigmoid <- function(z){
  return(1/(1+exp(-z)))
}

# Prediction step.


# Update step on tau_OBS.


# Update step on tau_ODE.
