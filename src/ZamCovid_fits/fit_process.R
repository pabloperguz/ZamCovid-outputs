ZamCovid_fit_process <- function(samples, parameters, data_full, data_fit) {
  
  region <- samples$info$region
  
  samples$trajectories$date <-
    samples$trajectories$time / samples$trajectories$rate
  date <- ZamCovid:::numeric_date_as_date(tail(samples$trajectories$date, 1))
  
  message("Computing Rt")
  rt <- calculate_ZamCovid_Rt(samples)
  
  message("Computing severity")
  severity <- calculate_severity(samples)
  
  message("Computing parameter MLE and covariance matrix")
  parameters_new <- fit_parameters(samples, parameters)
  
  
  list(fit = list(
    samples = samples,
    rt = rt,
    severity = severity,
    parameters = parameters_new,
    data = list(fitted = data_fit, full = data_full))
  )
}


calculate_severity <- function(samples) {
  
  time <- samples$trajectories$time
  dates <- samples$trajectories$date
  state <- samples$trajectories$state
  
  ret <- calculate_severity_region(state, time, dates)
  ret
}


calculate_severity_region <- function(state, time, date) {
  
  # String vectors to formulate severity trajectory names needed
  sev_traj <- grep("^ifr", rownames(state), value = TRUE)
  
  severity <- list()
  
  for (s in sev_traj) {
    tmp <- t(state[s, , ])
    severity[[s]] <- tmp
  }
  
  severity$time <- time
  severity$date <- date
  class(severity) <- "IFR_t_trajectories"
  
  severity
}


calculate_ZamCovid_Rt <- function(samples, weight_Rt = FALSE) {
  time <- samples$trajectories$time
  info <- samples$info$info
  epoch_dates <- samples$info$epoch_dates
  
  state <- samples$trajectories$state
  pars <- samples$pars
  transform <- samples$predict$transform
  
  
  ret <- calculate_ZamCovid_Rt_region(pars, state, transform,
                                      time, info, epoch_dates, weight_Rt)
  ret
}


calculate_ZamCovid_Rt_region <- function(pars, state, transform,
                                         time, info, epoch_dates, weight_Rt,
                                         keep_strains_Rt) {
  
  index_S <- grep("^S_", rownames(state))
  index_R <- grep("^R_", rownames(state))
  index_ps <- NULL
  
  S <- state[index_S, , , drop = FALSE]
  R <- state[index_R, , , drop = FALSE]
  
  dates <- time / 4
  
  n_pars <- nrow(pars)
  
  rt <- list(time = numeric(0),
             date = numeric(0),
             beta = numeric(0),
             eff_Rt_general = numeric(0),
             Rt_general = numeric(0))
  
    
    pars_model <- lapply(seq_rows(pars), function(j)
      transform(pars[j, ]))
    
    n_strains <- 1
    n_strains_R <- 1
    n_vacc_classes <- pars_model[[1]]$n_vacc_classes
    
    suffix <- paste0("_", c(ZamCovid:::model_age_bins()$start))
    S_nms <- get_names("S", list(n_vacc_classes), suffix)
    S1 <- S[S_nms, , , drop = FALSE]
    
    rt <- ZamCovid::ZamCovid_Rt_trajectories(
      time, S1, pars_model, R = NULL)

  class(rt) <- c("Rt_trajectories", "Rt")
  rt
}


fit_parameters <- function(samples, parameters) {

  keep <- function(x, region) {
    is.na(x) | x %in% region
  }
  
  region <- samples$info$region
  info <- parameters$raw$info[keep(parameters$raw$info$region, region), ]
  rownames(info) <- NULL
  
  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]
  info$initial[match(names(initial), info$name)] <- unname(initial)
  
  covariance <- cov(samples$pars)
  rownames(covariance) <- NULL
  proposal <- data_frame(region = samples$info$region,
                         name = colnames(covariance),
                         covariance)
  
  prior <- parameters$prior
  prior$region <- samples$info$region
  prior <- prior[, c("region", setdiff(names(prior), "region"))]
  rownames(prior) <- NULL
  
  
  parameters$info <- info
  parameters$prior <- prior
  parameters$proposal <- proposal
  
  parameters
}


get_names <- function(state_name, suffix_list, suffix0 = NULL) {
  if (is.null(suffix0)) {
    suffixes <- list()
  } else {
    suffixes <- list(suffix0)
  }
  for (i in seq_along(suffix_list)) {
    nm <- names(suffix_list)[[i]]
    if (length(nm) == 0) {
      nm <- ""
    }
    suffixes <- c(suffixes,
                  list(c("", sprintf("_%s%s", nm,
                                     seq_len(suffix_list[[i]] - 1L)))))
  }
  suffixes <- expand.grid(suffixes)
  nms <- apply(suffixes, 1,
               function(x) sprintf("%s%s",
                                   state_name, paste0(x, collapse = "")))
  nms
}
