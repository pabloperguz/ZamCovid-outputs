source("util.R")
version_check("ZamCovid", "0.1.0")


## 1. Prepare elements of particle filter
pars <- fit_pars_load("parameters", region, assumptions,
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
# dat <- ZamCovid_fit_process(samples, pars, data_full, data_fit)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat, "outputs/fit.rds")

message("Creating plots")
png("outputs/pmcmc_traceplots.png", units = "in", width = 8, height = 8, res = 300)
plot_fit_traces(samples)
dev.off()
  
png("outputs/fits_serology.png", units = "in", width = 16, height = 10, res = 300)
plot_serology(samples, data_fit)
dev.off()

png("outputs/fits_deaths.png", units = "in", width = 16, height = 10, res = 300)

dev.off()


tmp <- samples$trajectories$state
tmp_d <- tmp[grep("D_", rownames(tmp)), , ]
plot(colMeans(as.data.frame(tmp_d[1, , -1L])), type = "l", 
     ylab = "Deaths")
lines(colMeans(as.data.frame(tmp[16, , -1L])), col = "red")
rownames(tmp)