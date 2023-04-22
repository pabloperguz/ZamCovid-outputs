source("util.R")
version_check("ZamCovid", "0.1.0")


## 1. Prepare elements of particle filter
pars <- fit_pars_load("inputs", region, assumptions,
                      short_run, deterministic)

control <- set_control(short_run, deterministic)

data_full <- read_csv("data/data_timeseries.csv")
data_fit <- parse_data(data_full,
                       fit_sero = TRUE, fit_deaths = TRUE,
                       sero_by_age = TRUE, deaths_by_age = TRUE)


## 2. Build particle filter and run pMCMC
filter <- ZamCovid_particle_filter(data_fit, pars$mcmc,
                                   control$particle_filter,
                                   deterministic = deterministic)

samples <- fit_run(pars, filter, control$pmcmc)


## 3. Post-processing of model fits ----
dat <- ZamCovid_fit_process(samples, pars, data_full, data_fit)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat, "outputs/fit.rds")
write_csv(dat$fit$parameters$info, "outputs/info.csv")
write_csv(dat$fit$parameters$proposal, "outputs/proposal.csv")

dir.create("plots", FALSE, TRUE)
message("Creating plots")
png("plots/pmcmc_traceplots.png", units = "in", width = 10, height = 6, res = 300)
plot_fit_traces(samples)
dev.off()

png("plots/fits_serology.png", units = "in", width = 10, height = 6, res = 300)
plot_serology(samples, data_fit)
dev.off()

png("plots/fits_deaths.png", units = "in", width = 10, height = 6, res = 300)
plot_deaths(samples, data_fit, age = FALSE)
dev.off()