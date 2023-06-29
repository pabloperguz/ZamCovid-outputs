# pkgload::load_all("~/R_Projects/ZamCovid", export_all = FALSE)

## ---------------------------
## Develop and test locally
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------
short_run <- TRUE
assumptions <- "central"
districts <- c("kabwe", "lusaka", "livingstone", "ndola", "solwezi")
deterministic <- TRUE
env_keep <- c("root_dir", "short_run", "assumptions", "districts",
              "deterministic", "env_keep")
# ---------------------------------


## I. ZamCovid_multidistrict_data -----

# Develop
orderly::orderly_develop_start("ZamCovid_multidistrict_data",
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_multidistrict_data"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_multidistrict_data",
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## II. ZamCovid_multidistrict_parameters ----

# Develop
orderly::orderly_develop_start("ZamCovid_multidistrict_parameters",
                               parameters = list(assumptions = assumptions,
                                                 deterministic = deterministic),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_multidistrict_parameters"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_multidistrict_parameters",
                     parameters = list(assumptions = assumptions,
                                       deterministic = deterministic),
                     use_draft = "newer")
rm(list = setdiff(ls(), env_keep))

#----


## III. ZamCovid_multidistrict_fits ----

# Develop
orderly::orderly_develop_start("ZamCovid_multidistrict_fits",
                               parameters = list(district = "lusaka",
                                                 short_run = short_run,
                                                 assumptions = assumptions,
                                                 deterministic = deterministic),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_multidistrict_fits"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
lapply(districts,
       function(r) orderly::orderly_run("ZamCovid_multidistrict_fits",
                                        parameters = list(district = r,
                                                          short_run = short_run,
                                                          assumptions = assumptions,
                                                          deterministic = deterministic),
                                        use_draft = "newer"))
rm(list = setdiff(ls(), env_keep))

#----


## IV. ZamCovid_kabwe_sens_analysis_multidistrict ----

# Develop
orderly::orderly_develop_start("ZamCovid_kabwe_sens_analysis",
                               parameters = list(date = date,
                                                 short_run = short_run,
                                                 deterministic = deterministic),
                               use_draft = "newer")
setwd(paste0(root_dir, "ZamCovid_kabwe_sens_analysis"))
file.edit("script.R")
# tidy up
orderly::orderly_develop_clean()
rm(list = setdiff(ls(), env_keep))

# Run
orderly::orderly_run("ZamCovid_kabwe_sens_analysis",
                     parameters = list(date = date,
                                       short_run = short_run,
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
# remotes::install_github("mrc-ide/ZamCovid",
#                         ref = "master",
#                         force = TRUE)
# 
# remotes::install_github("mrc-ide/ZamCovid",
#                         ref = "master",
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
districts <- c("kabwe", "lusaka", "livingstone", "ndola", "solwezi")

#----

fits <- 
  obj$lapply(X = districts,
             FUN = function (x) {
               orderly::orderly_run('ZamCovid_fits_multidistrict',
                                    parameters = list(district = r,
                                                      short_run = FALSE,
                                                      deterministic = TRUE,
                                                      assumptions = "central"),
                                    use_draft = "newer")})
batch <- fits$name

res <- obj$task_bundle_get(batch)$results()

fits_result <- fits$result()
