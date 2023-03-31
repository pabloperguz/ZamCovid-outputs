## Pipeline of tasks to run Kabwe analysis

# pkgload::load_all("~/R_Projects/ZamCovid", export_all = FALSE)

## ---------------------------
## Develop and test locally
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------
short_run <- TRUE
assumptions <- "central"
env_keep <- c("root_dir", "short_run", "assumptions", "env_keep")
# ---------------------------------


## 1. ZamCovid_data ----

# Develop
orderly::orderly_develop_start("ZamCovid_data", use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_data"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_data", use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## 2. ZamCovid_parameters ----

# Develop
orderly::orderly_develop_start(
  "ZamCovid_parameters", parameters = list(assumptions = assumptions),
  use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_parameters"))
file.edit("script.R")
# Tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("severity_parameters", parameters = list(
  deterministic = deterministic, assumptions = assumptions),
  use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----

## 2. severity_fits ----

# Develop
orderly::orderly_develop_start("severity_fits",
                               parameters = list(region = "london",
                                                 short_run = short_run,
                                                 deterministic = deterministic,
                                                 assumptions = assumptions),
                               use_draft = "newer")
setwd(paste0(root_dir, "severity_fits"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run all regions
lapply(regions,
       function(r) orderly::orderly_run("severity_fits",
                                        parameters = list(region = r,
                                                          short_run = short_run,
                                                          deterministic = deterministic,
                                                          assumptions = assumptions),
                                        use_draft = "newer"))
rm(list = setdiff(ls(), env_keep))

#----

# 1. Data to fit
