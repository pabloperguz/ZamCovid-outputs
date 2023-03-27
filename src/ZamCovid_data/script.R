source("util.R")

end_date <- "2021-12-31"

## 1. Aggregate all data into time series
deaths <- read_csv("data/kabwe_all_deaths.csv")
serology <- read_csv("data/kabwe_seroprevalence.csv")
cases <- read_csv("data/kabwe_cases.csv")

# Only keep data up to 2021-12-31
# Ignore PCR data before 2021 as official policy for testing was different
data <- data %>%
  filter(as.Date(date) <= as.Date(end_date))

cols_pcr <- c(paste0("pcr_positive_", pcr_age_bands))
data[which(as.Date(data$date) < "2021-01-01"), cols_pcr] <- NA_integer_

write.csv(data, "data_timeseries.csv", row.names = FALSE)



cases <- readxl::read_excel("data/Data sets.xlsx", sheet = "COVID-19 Cases") %>%
  `colnames<-`(tolower(gsub(" ", "_", colnames(.))))

cases1 <- cases %>%
  select(date = date_specimen_collected, age, locality) %>%
  filter(!is.na(date), !is.na(age)) %>%
  mutate(locality = tolower(gsub(" ", "_", locality)))
