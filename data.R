################################
########## REAL DATA ###########
################################

directory <- "~/home/claire/Documents/data_covid_dashboard/sources-csv/data" # directory of data
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

countsplot <- ggplot(observations, aes(x = t, y = R)) +
  geom_line(aes(color = 'R')) +
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

region_filename <- paste0("real_data_", region, "_", daily_or_weekly, ".Rdata")
population_filename <- paste0("real_pop_", region, "_", daily_or_weekly, ".Rds")
save(observations, file = file.path(directory_save, region_filename))
saveRDS(population, file = file.path(directory_save, population_filename))