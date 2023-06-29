source("util.R")

cull_date <- as.Date("2020-09-30")
## 1. Aggregate all data into time series
data <- read_csv("data/data_multidistrict.csv") %>%
  mutate(date = as.Date(date)) %>%
  filter(date <= cull_date)

historic_deaths <- read_csv("data/historic_deaths.csv")

data_timeseries <- list(
  data = data,
  historic = historic_deaths
)

saveRDS(data_timeseries, "data_timeseries.rds")

png("data_to_fit.png", units = "in", width = 14, height = 10, res = 300)
plot_data(data)
dev.off()
