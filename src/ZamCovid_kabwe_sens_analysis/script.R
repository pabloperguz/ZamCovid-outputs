source("util.R")
version_check("ZamCovid", "0.1.0")

dir.create("outputs", FALSE, TRUE)

# Read scenario inputs
scenarios <- load_scenarios("fits/")

# Plot sensitivity analyses
rmarkdown::render("sens_analysis_results.Rmd")
