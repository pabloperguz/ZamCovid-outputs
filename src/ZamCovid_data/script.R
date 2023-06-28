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
# - PC R data is ignored before 2021, as official policy for testing was different
# - Deaths are "ALL DEATHS", not just Covid-19
data <- data.frame(date = seq.Date(start_date, end_date, "day")) %>%
  left_join(., serology, by = "date") %>%
  left_join(., cases, by = "date") %>%
  left_join(., deaths, by = "date")

# Mulenga et al. report Kabwe prevalence of 6% in July (all ages) n = 646
# https://doi.org/10.1016/S2214-109X(21)00053-X
data$sero_pos_all <- NA_real_
data$sero_tot_all <- NA_real_
data[data$date == "2020-07-15",
     c("sero_pos_all", "sero_tot_all")] <- c(39, 646)

data <- list(
  data = data,
  historic = historic_deaths)

saveRDS(data, "data_timeseries.rds")


png("data_to_fit.png", units = "in", width = 14, height = 10, res = 300)
plot_data(data$data)
dev.off()
