source("util.R")

version_check("ZamCovid", "0.1.0")

## Define date at which the data is capped for analysis
date <- date

## Regions: at the moment we only have Kabwe
regions <- c("kabwe")

## After starting with a model without vaccination
## * early December 2020: vaccination starts, expand vaccine classes 
epoch_dates <- c("2020-12-05")

## To use historic deaths timeseries to infer baseline deaths in 2020-2021, we'll
## read in the data object from ZamCovid_data
data_timeseries <- readRDS("data/data_timeseries.rds")
historic_deaths <- data_timeseries$historic

## Save all other data for use in the fitting task
data_timeseries <- data_timeseries$data
write_csv(data_timeseries, "data_timeseries.csv")

## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(assumptions, deterministic)

baseline <- lapply(regions, create_baseline, date, 
                   epoch_dates, pars, assumptions,
                   historic_deaths = historic_deaths)
names(baseline) <- regions

message("Writing parameters_info.csv")
write_csv(pars$info, "parameters_info.csv")
message("Writing parameters_proposal.csv")
write_csv(pars$proposal, "parameters_proposal.csv")
message("Writing parameters_prior.csv")
write_csv(pars$prior, "parameters_prior.csv")

message("Writing parameters_base.rds")
saveRDS(baseline, "parameters_base.rds")

message("Writing parameters_transform.R")
fs::file_copy("R/transform.R",
              "parameters_transform.R", overwrite = TRUE)
