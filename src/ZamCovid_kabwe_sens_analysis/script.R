source("util.R")
version_check("ZamCovid", "0.1.0")

# Read scenario inputs
scenarios <- load_scenarios("fits/")

# Outputs comparative plots
png("scenarios_data.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, c("fit_deaths_only", "fit_serology_only"))
dev.off()

png("scenarios_imm_waning.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, c("imm_waning_fast", "imm_waning_slow"))
dev.off()

png("scenarios_ifr.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, c("ifr_high", "ifr_low"))
dev.off()

png("scenarios_seroreversion.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, c("serorev_slow", "serorev_fast"))
dev.off()

png("scenarios_serology_sens.png", units = "in", width = 8, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, c("sero_sens_low", "sero_sens_high"))
dev.off()

png("scenarios_all.png", units = "in", width = 16, height = 8, res = 300)
plot_trajectory_scenarios(scenarios, names(scenarios$trajectories)[-1])
dev.off()
