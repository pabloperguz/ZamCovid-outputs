source("util.R")

version_check("ZamCovid", "0.1.1")
dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)

dat <- combined_load("regional_results")
districts <- names(dat$samples)

par_names <- unique(dat$parameters$proposal$name)
subset_misc <- grep("^beta", par_names, value = TRUE, invert = TRUE)
par_labels <- forest_plot_labels(dat)


png("figs/forest_plot_misc.png", width = 10, height = 10, units = "in", res = 200)
plot_forest(dat, plot_type = "non_betas", par_labels = par_labels)
dev.off()

png("figs/forest_plot_betas.png", width = 10, height = 10, units = "in", res = 200)
plot_forest(dat, plot_type = "betas", par_labels = par_labels)
dev.off()

png("figs/data_fits.png", width = 12, height = 10, units = "in", res = 200)
plot_trajectories(dat, districts, age_band = "all")
dev.off()

png("figs/status_effective_susceptible.png", width = 10, height = 10, units = "in", res = 200)
plot_effective_susceptible(dat, districts)
dev.off()

png("figs/status_infection.png", width = 10, height = 10, units = "in", res = 200)
plot_infection_status(dat, districts)
dev.off()

png("figs/cumulative_attack_rate.png", width = 10, height = 10, units = "in", res = 200)
plot_cumulative_attack_rate(dat, districts)

png("figs/incidence.png", width = 10, height = 10, units = "in", res = 200)
plot_incidence(dat, districts)
dev.off()

png("figs/Rt_beta.png", width = 10, height = 10, units = "in", res = 200)
plot_Rt(dat, districts)
dev.off()

pdf("plots/manuscript_figure_1.pdf", width = 16.5, height = 10)
manuscript_plot_1(dat)
dev.off()

pdf("plots/manuscript_figure_2.pdf", width = 16.5, height = 10)
manuscript_plot_2(dat)
dev.off()

png("manuscript_plots/suppl_outcomes_demography.png", units = "in", width = 10, height = 10, res = 300)
suppl_outcomes_demography(dat)
dev.off()

## Render rmd
rmarkdown::render("manuscript_numbers.Rmd")
