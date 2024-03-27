################################
########## REAL DATA ###########
################################

directory <- "~/Documents/GitHub/proboder/Data" # directory for data
region <- 'GE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'daily' # choose either 'daily' or 'weekly' (if 'real' data selected)

# Import dataset.
setwd("~/Documents/GitHub/proboder/data_covid_dashboard/sources-csv/data")

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
colnames(observations) <- c('date','S','I','D')

n <- nrow(observations) # size of data_grid

################################
######## SAVE THE DATA #########
################################

region_filename <- paste0("real_data_", region, "_", daily_or_weekly, ".Rdata")
population_filename <- paste0("real_pop_", region, "_", daily_or_weekly, ".Rds")
save(observations, file = file.path(directory, region_filename))
saveRDS(population, file = file.path(directory, population_filename))