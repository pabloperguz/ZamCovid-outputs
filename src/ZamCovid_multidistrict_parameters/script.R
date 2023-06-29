source("util.R")

version_check("ZamCovid", "0.1.1")

data <- readRDS("data/data_timeseries.rds")
data_multidistrict <- data$data
districts <- unique(data_multidistrict$district)
date <- max(data_multidistrict$date)

## Cap data and save for use in the fitting task
stopifnot(max(data_multidistrict$date) >= as.Date(date))
data_multidistrict <- data_multidistrict[data_multidistrict$date <= as.Date(date), ]
write_csv(data_multidistrict, "data_multidistrict.csv")

## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(assumptions, deterministic)

baseline <- lapply(districts, create_baseline, date, pars, assumptions)
names(baseline) <- districts

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
