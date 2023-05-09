source("util.R")

start_date <- as.Date("2020-01-01")
end_date <- as.Date("2021-12-31")

## 1. Aggregate all data into time series
deaths <- read_csv("data/kabwe_all_deaths.csv") %>% mutate(date = as.Date(date))
serology <- read_csv("data/kabwe_seroprevalence.csv") %>% mutate(date = as.Date(date))
cases <- read_csv("data/kabwe_cases.csv") %>% mutate(date = as.Date(date))
historic_deaths <- get_baseline_deaths()

## Data will run from 2020-01-01 to 2021-12-31
# Note that:
# - Serology is weekly
# - PCR data is ignored before 2021, as official policy for testing was different
# - Deaths are "ALL DEATHS", not just Covid-19
data <- data.frame(date = seq.Date(start_date, end_date, "day")) %>%
  left_join(., serology, by = "date") %>%
  left_join(., cases, by = "date") %>%
  left_join(., deaths, by = "date")

data <- list(
  data = data,
  historic = historic_deaths)

saveRDS(data, "data_timeseries.rds")


png("data_to_fit.png", units = "in", width = 14, height = 10, res = 300)
plot_data(data$data)
dev.off()
