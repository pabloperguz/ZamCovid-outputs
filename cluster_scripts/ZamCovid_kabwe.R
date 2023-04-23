## Pipeline of tasks to run Kabwe analysis

# pkgload::load_all("~/R_Projects/ZamCovid", export_all = FALSE)

## ---------------------------
## Develop and test locally
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------
short_run <- TRUE
date <- "2021-03-01"
assumptions <- "central"
deterministic <- TRUE
env_keep <- c("root_dir", "short_run", "date", "assumptions",
              "deterministic", "env_keep")
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
  "ZamCovid_parameters",
  parameters = list(assumptions = assumptions, date = date,
                    deterministic = deterministic),
  use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_parameters"))
file.edit("script.R")
# Tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_parameters", parameters = list(
  assumptions = assumptions, date = date, deterministic = deterministic),
  use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## 3. ZamCovid_fits ----

# Develop
orderly::orderly_develop_start("ZamCovid_fits",
                               parameters = list(region = "kabwe",
                                                 date = date,
                                                 short_run = short_run,
                                                 assumptions = assumptions,
                                                 deterministic = deterministic),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_fits"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_fits",
                     parameters = list(region = "kabwe",
                                       date = date,
                                       short_run = FALSE,
                                       assumptions = assumptions,
                                       deterministic = deterministic),
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----

# To remove failed tasks from draft folder
orderly::orderly_cleanup(failed_only = TRUE)
