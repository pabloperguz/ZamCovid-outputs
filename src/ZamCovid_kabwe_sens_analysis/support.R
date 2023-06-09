load_scenarios <- function(dir) {
  
  scenarios <- list.files(dir)
  
  out <- purrr::map(scenarios, function(x)  extract_trajectories(dir, x))
  
  nms <- scenarios %>% gsub("_fits.rds", "", .)
  names(out) <- nms
  
  switch_levels(out)
}


extract_trajectories <- function(dir, file) {
  
  scen <- readRDS(file = paste0(dir, file))$fit
  traj_nms <- c("deaths_all_inc", "sero_pos_over15")
  
  pop <- sum(scen$parameters$base$population[, "n"][4:16])
  tmp <- lapply(traj_nms, function(x) summarise_trajectory(scen$samples, x, pop))
  tmp <- data.table::rbindlist(tmp)
  trajectories <- tmp
  
  severity <- summarise_trajectory(scen$severity, "ifr")
  
  eff_Rt_general <- summarise_trajectory(scen$rt, "eff_Rt_general")
  Rt_general <- summarise_trajectory(scen$rt, "Rt_general")
  rt <- rbind(eff_Rt_general, Rt_general)

  list(
    trajectories = trajectories,
    severity = severity,
    rt = rt
  )
}


summarise_trajectory <- function(sample, nm, pop = NULL) {
  
  if (nm == "ifr") {
    dates <- sample$date[-1]
    state <- sample$ifr[-1, ]
  } else if (nm %in% c("eff_Rt_general", "Rt_general")) {
    dates <- sample$date[-1, 1]
    state <- sample[[nm]][-1, ]
  } else if (nm == "sero_pos_over15") {
    dates <- sample$trajectories$date[-1]
    state <- t(sample$trajectories$state[nm, , -1]) / pop
  } else {
    dates <- sample$trajectories$date[-1]
    state <- t(sample$trajectories$state[nm, , -1])
  }
  
  data.frame(
    date = ZamCovid:::numeric_date_as_date(dates),
    state = nm,
    mean = rowMeans(state),
    lb = matrixStats::rowQuantiles(state, probs = 0.025),
    ub = matrixStats::rowQuantiles(state, probs = 0.975)
  )
}


switch_levels <- function(x) {
  nms <- names(x[[1]]) %||% seq_along(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}


unlist_scenarios <- function(x, what) {
  out <- purrr::map_df(x[[what]], ~as.data.frame(.x), .id = "scenario")
  out$scenario <- relevel(factor(out$scenario), "central")
  out
}
