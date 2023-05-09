## Pipeline of tasks to run Kabwe analysis

# pkgload::load_all("~/R_Projects/ZamCovid", export_all = FALSE)

## ---------------------------
## Develop and test locally
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------
short_run <- TRUE
date <- "2021-08-01"
assumptions <- "central"
## can be: central, fit_no_deaths, low_ifr, base_deaths_low, base_deaths_high
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
                                       short_run = short_run,
                                       assumptions = assumptions,
                                       deterministic = deterministic),
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----

# To remove failed tasks from draft folder
orderly::orderly_cleanup(failed_only = TRUE)


############################################################
### Run fits task in cluster

## 0. Install ZamCovid ----
# Sys.unsetenv("GITHUB_PAT")
# remotes::install_github("pabloperguz/ZamCovid",
#                         ref = "baseline-deaths",
#                         force = TRUE)
# 
# remotes::install_github("pabloperguz/ZamCovid",
#                         ref = "baseline-deaths",
#                         force = TRUE,
#                         lib = "Z:/Pablo/ZamCovid-outputs/contexts/lib/windows/4.2")

#----


## 1. Context ----
setwd(orderly::orderly_config()$root)
packages <- c("lubridate", "coda", "tidyr", "ggplot2",
              "viridisLite", "orderly", 'vaultr', 'readxl', "ggtext",
              'abind', 'here', "mcstate", "dust", "purrr", "patchwork",
              "stringr", "desplot", "rmarkdown", "ZamCovid", "eigen1",
              "jtools", "DescTools", "car", "data.table", "matrixStats",
              "Hmisc", "scales")
src <- conan::conan_sources(NULL,
                            repos = c("https://ncov-ic.github.io/drat",
                                      "https://raphaels1.r-universe.dev",
                                      "https://pabloperguz/ZamCovid"))
ctx <- context::context_save("contexts",
                             packages = packages,
                             package_sources = src)
cfg <- didehpc::didehpc_config(cluster = "wpia-hn",
                               template = 'AllNodes',
                               cores = 8)
obj <- didehpc::queue_didehpc(ctx, config = cfg)

#----


fits <- obj$enqueue(orderly::orderly_run('ZamCovid_fits',
                                         parameters = list(region = "kabwe",
                                                           date = "2021-08-01",
                                                           short_run = FALSE,
                                                           deterministic = TRUE,
                                                           assumptions = "central"),
                                         use_draft = "newer"))
fits_result <- fits$result()
