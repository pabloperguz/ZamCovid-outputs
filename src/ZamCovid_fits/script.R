source("util.R")
version_check("ZamCovid", "0.1.0")


## 1. Prepare elements of particle filter
pars <- fit_pars_load("inputs", region, assumptions,
                      short_run, deterministic)

control <- set_control(short_run, deterministic)

data_full <- read_csv("data/data_timeseries.csv")
data_fit <- parse_data(data_full,
                       fit_sero = (!assumptions == "fit_deaths_only"),
                       sero_by_age = TRUE,
                       fit_deaths = (!assumptions == "fit_serology_only"))

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

png("plots/fits_ifr.png", units = "in", width = 6, height = 6, res = 300)
plot_severity(dat)
dev.off()

png("plots/fits_rt.png", units = "in", width = 6, height = 6, res = 300)
plot_rt(dat)
dev.off()

png("plots/pmcmc_traceplots.png", units = "in", width = 16, height = 12, res = 300)
plot_fit_traces(dat$fit$samples)
dev.off()

png("plots/fits_serology.png", units = "in", width = 8, height = 8, res = 300)
plot_serology(dat$fit$samples, data_full)
dev.off()

png("plots/fits_deaths.png", units = "in", width = 10, height = 10, res = 300)
plot_deaths(dat, data_fit) | plot_deaths_disag(dat)
dev.off()

png("plots/infections_inc.png", units = "in", width = 8, height = 10, res = 300)
plot_infection_incidence(dat)
dev.off()

# Epidemics9 plot
col1 <- (
  (plot_serology(dat$fit$samples, data_full, over15_only = TRUE) +
     theme(axis.text.x = element_blank())) /
    (plot_rt(dat) + theme(axis.text.x = element_blank())) /
    plot_infection_incidence(dat)) +
  plot_layout(heights = c(0.3, 0.3, 1))

col2 <- ((plot_deaths(dat, data_fit, week_only = TRUE) +
  theme(axis.text.x = element_blank(),
        legend.position = "none")) /
  (plot_severity(dat, age = FALSE) +
     theme(axis.text.x = element_blank())) / 
  plot_deaths_disag(dat, plot_age = FALSE)) +
  plot_layout(heights = c(0.3, 0.3, 1))

png("plots/epidemics.png", units = "in", width = 12, height = 12, res = 300)
(col1 | col2) +
  plot_annotation(tag_levels = "A")
dev.off()


rmarkdown::render("fit_results.Rmd")
