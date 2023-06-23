source("util.R")
version_check("ZamCovid", "0.1.0")

pkgload::load_all("~/R_Projects/ZamCovid", export_all = FALSE)
## 1. Prepare elements of particle filter
pars <- fit_pars_load("inputs", district, assumptions,
                      short_run, deterministic)

control <- set_control(short_run, deterministic)

data_full <- read_csv("data/data_timeseries.csv")

data_fit <- parse_data(data_full, district, fit_sero = TRUE, fit_deaths = TRUE,
                       fit_cases = TRUE)

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


## 4. Plot fitted trajectories ----
sero <- plot_serology(dat, data_fit, all = TRUE, 0.4)
pcr <- plot_pcr_positivity(dat, data_fit, all = TRUE, 0.2) +
  theme(axis.text.x = element_blank())
deaths <- plot_deaths(dat, data_fit)

png("plots/fits.png", units = "in", width = 6, height = 8, res = 300)
(pcr / sero) | deaths
dev.off()

png("plots/inferred.png", units = "in", width = 6, height = 8, res = 300)
(plot_rt(dat) + theme(legend.position = "top")) / plot_severity(dat, age = FALSE)
dev.off()

png("plots/traceplots.png", units = "in", width = 16, height = 12, res = 300)
plot_fit_traces(dat$fit$samples)
dev.off()
