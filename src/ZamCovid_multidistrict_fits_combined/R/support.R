combined_load <- function(dir, get_onward = FALSE) {
  districts <- c("kabwe", "lusaka", "livingstone", "ndola", "solwezi")
  
  fits <- file.path(dir, districts, "fit.rds")
  
  msg <- !file.exists(fits)
  if (any(msg)) {
    msg <- sprintf(" - %s", file.path(districts[msg], "fit.rds"))
    stop(sprintf("Missing fits at '%s': \n%s",
                 path, paste(msg, collapse = "\n")),
         call. = FALSE)
  }
  
  read_rds <- function(filename, d) {
    message(sprintf("Reading fit for %s", d))
    readRDS(filename)$fit
  }
  dat <- Map(read_rds, fits, districts)
  names(dat) <- districts
  ret <- switch_levels(dat)
  
  message("Creating new info, proposal and prior files")
  dir.create("outputs/parameters", FALSE, TRUE)
  ret$parameters <- switch_levels(ret$parameters)
  pars_combine <- c("info", "prior", "proposal")
  ret$parameters[pars_combine] <- 
    lapply(ret$parameters[pars_combine], dplyr::bind_rows)
  write_csv(ret$parameters$info, "outputs/parameters/info.csv")
  write_csv(ret$parameters$proposal, "outputs/parameters/proposal.csv")
  write_csv(ret$parameters$prior, "outputs/parameters/prior.csv")
  
  if (get_onward) {
    message("Getting simulate object")
    ret$onward <- combined_onward(ret)
  }
  
  ret
}


combined_onward <- function(dat) {
  date <- dat$samples[[1]]$info$date
  steps_per_day <- dat$samples[[1]]$info$data$steps_per_day
  list(date = date,
       step = ZamCovid:::numeric_date(date) * steps_per_day,
       steps_per_day = steps_per_day,
       dt = 1 / steps_per_day,
       pars = lapply(dat$samples, "[[", "pars"),
       base = lapply(dat$parameters, function(x) x$base),
       state = lapply(dat$samples, "[[", "state"),
       data = lapply(dat$samples, function(x) x$predict$filter$data),
       transform = lapply(dat$samples, function(x) x$predict$transform),
       info = lapply(dat$samples, "[[", "info"),
       vaccine = lapply(dat$samples, "[[", "vaccine"),
       simulate = combined_onward_simulate(dat))
}


combined_onward_simulate <- function(dat) {
  browser()
  simulate <- list_transpose(dat$simulate)
  
  dates <- dat$simulate[[1]]$date
  idx_dates <- dat$samples[[1]]$trajectories$date %in% dates
  state_names <- unique(c(rownames(simulate$state[[1]]), "deaths_hosp",
                          "sero_pos_1", "sero_pos_2"))
  state_names <- intersect(state_names,
                           rownames(dat$samples[[1]]$trajectories$state))
  
  state <- lapply(dat$samples, function(x)
    x$trajectories$state[state_names, , idx_dates])
  state <- aperm(abind_quiet(state, along = 4), c(1, 2, 4, 3))
  
  state_by_age <- lapply(list_transpose(simulate$state_by_age),
                         abind_quiet, along = 3)
  
  ret <- list(date = dates,
              state = state,
              state_by_age = state_by_age)
  
  ## This is not terrible:
  rt <- list_transpose(dat$rt)[c("Rt_general", "eff_Rt_general")]
  ## Rt_general and eff_Rt_general will have dimensions:
  ## [n particles x n regions x n dates]
  rt_combined <- lapply(rt, function(x)
    aperm(abind_quiet(x, along = 3), c(2, 3, 1))[, , idx_dates])
  
  idx_variant_dates <- dat$variant_rt[[1]]$date[, 1] %in% dates
  
  variant_rt <-
    list_transpose(dat$variant_rt)[c("Rt_general", "eff_Rt_general")]
  variant_rt_combined <- lapply(variant_rt, function(x)
    aperm(abind_quiet(x, along = 4), c(3, 4, 1, 2))[, , idx_variant_dates, ])
  
  rt_combined <- Map(function(rt, weighted_rt) {
    x <- abind_quiet(rt, weighted_rt, along = 4)
    dimnames(x)[[4]] <- c("strain_1", "strain_2", "both")
    x}, variant_rt_combined, rt_combined)
  
  ret <- c(ret, rt_combined)
  
  idx_dates_mv_rt <- dat$variant_rt[[1]]$date[, 1] %in% dates
  
  ## multivariant_Rt_general and multivariant_eff_Rt_general will have
  ## dimensions: [n particles x n regions x n variants x n dates]
  ## TODO: we are putting "multivariant" in the name here for clarity,
  ## so we may want to rename variant_rt wherever it appears to
  ## multivariant_rt
  mv_rt <-
    list_transpose(dat$variant_rt)[c("Rt_general", "eff_Rt_general")]
  mv_rt_combined <- lapply(mv_rt, function(x)
    aperm(abind_quiet(x, along = 4), c(3, 4, 2, 1))[, , , idx_dates_mv_rt])
  names(mv_rt_combined) <- paste0("multivariant_", names(mv_rt_combined))
  
  ret <- c(ret, mv_rt_combined)
  
  ret
}
