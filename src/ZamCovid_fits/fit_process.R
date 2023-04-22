
ZamCovid_fit_process <- function(samples, parameters, data_full, data_fit,
                                 Rt = TRUE, simulate = FALSE, severity = FALSE) {
  
  region <- samples$info$region
  message("Computing restart information")
  restart <- fit_process_restart(samples, parameters)
  samples$restart <- NULL
  
  samples$trajectories$date <-
    samples$trajectories$time / samples$trajectories$rate
  date <- ZamCovid:::numeric_date_as_date(tail(samples$trajectories$date, 1))
  
  ## The Rt calculation is slow and runs in serial
  # if (Rt) {
  #   message("Computing Rt")
  #   rt <- calculate_ZamCovid_Rt(samples, TRUE)
  # } else {
  #   rt <- NULL
  # }
  
  message("Computing parameter MLE and covariance matrix")
  parameters_new <- fit_parameters(samples, parameters)
  
  
  ## Output simulate object
  # if (simulate) {
  #   message("Preparing onward simulation object")
  #   simulate <- create_simulate_object(samples, date, samples$info$date)
  # } else {
  #   simulate <- NULL
  # }
  # 
  # 
  # ## Extract demography
  # if (severity) {
  #   message("Extracting demography")
  #   model_demography <- extract_demography(samples)
  #   message("Extracting severity outputs")
  #   severity <- extract_severity(samples)
  #   message("Calculating intrinsic severity")
  #   intrinsic_severity <- calculate_intrinsic_severity(samples, parameters$base)
  # } else {
  #   model_demography <- NULL
  #   severity <- NULL
  #   intrinsic_severity <- NULL
  # }
  
  
  ## This can possibly merge in with the initial restart processing if
  ## we're careful now.  We need to have computed Rt first is the only
  ## trick.
  # if (!is.null(restart)) {
  #   ## When adding the trajectories, we might as well strip them down
  #   ## to the last date in the restart
  #   restart_date <- max(restart$state$time)
  #   i <- samples$trajectories$date <= restart_date
  #   
  #   if (Rt) {
  #     parent_rt <- rt_filter_time(rt, i, samples$info$multiregion)
  #   } else {
  #     parent_rt <- NULL
  #   }
  #   
  #   restart$parent <- list(
  #     trajectories = trajectories_filter_time(samples$trajectories, i),
  #     rt = parent_rt,
  #     data = data,
  #     prior = parameters$raw$prior)
  # }
  # 
  
  list(
    fit = list(samples = samples,
               rt = NULL,
               severity = NULL,
               simulate = NULL,
               parameters = parameters_new,
               data = list(fitted = data_fit, full = data_full)),
    restart = restart)
}


fit_process_restart <- function(samples, parameters) {
  if (is.null(samples$restart)) {
    return(NULL)
  }
  
  ## TODO: this is now done here and also later, but it's fairly fast
  pars <- fit_parameters(samples, parameters)
  ## TODO: could pass in samples$info$pars, but somehow feels more
  ## awkward to do so
  pars$prior <- fit_process_restart_priors(samples$pars, pars)
  pars$sample <- samples$pars
  class(pars) <- "pars_pmcmc"
  
  pars$base <- parameters$base
  
  list(state = samples$restart,
       info = samples$info,
       pars = pars)
}


## TODO: not tested yet
create_simulate_object <- function(samples, start_date_sim, date) {
  start_date_sim <- ZamCovid:::numeric_date(start_date_sim)
  fit_dates <- samples$trajectories$date
  idx_dates <- (fit_dates >= start_date_sim) &
    (fit_dates <= ZamCovid:::numeric_date(date))
  date <- fit_dates[idx_dates]
  
  state_keep <- c("deaths", "admitted",
                  "diagnoses", "infections", "hosp", "icu")
  state_full <- samples$trajectories$state
  
  if (samples$info$multiregion) {
    region <- samples$info$region
    state_by_age <- lapply(region, function(r)
      extract_age_class_state(state_full[, r, , idx_dates]))
    names(state_by_age) <- region
    
    ## thin trajectories, but here we already have a regional dimension
    state <- state_full[state_keep, , , idx_dates]
    ## However, we have to push it out
  } else {
    state_by_age <- extract_age_class_state(state_full[, , idx_dates])
    
    ## thin trajectories and reshape to add a regional dimension:
    state <- state_full[state_keep, , idx_dates]
    n_samples <- ncol(state_full)
    ## TODO: this needs to be moved to be (1, n_samples) (swapping
    ## dimensions) to match the multiregion filter.
    state <- mcstate::array_reshape(state, i = 2, c(n_samples, 1))
  }
  
  list(date = date, state = state, state_by_age = state_by_age)
}


calculate_ZamCovid_Rt <- function(samples, weight_Rt) {
  browser()
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


## TODO: not tested yet
calculate_ZamCovid_Rt_region <- function(pars, state, transform,
                                         time, info, epoch_dates, weight_Rt,
                                         keep_strains_Rt) {
  browser()
  index_S <- grep("^S_", rownames(state))
  index_R <- grep("^R_", rownames(state))
  index_ps <- grep("^prob_strain", rownames(state))
  
  S <- state[index_S, , , drop = FALSE]
  R <- state[index_R, , , drop = FALSE]
  prob_strain <- state[index_ps, , , drop = FALSE]
  
  dates <- time / 4
  
  n_pars <- nrow(pars)
  
  type <- c("eff_Rt_general", "Rt_general")
  
  
  rt <- list(time = numeric(0),
             date = numeric(0),
             beta = numeric(0),
             eff_Rt_general = numeric(0),
             Rt_general = numeric(0))
  
  for (i in seq_len(length(epoch_dates) + 1L)) {
    if (i == 1) {
      dates1 <- which(dates <= epoch_dates[1])
      initial_time_from_parameters <- TRUE
    } else  if (i <= length(epoch_dates)) {
      dates1 <- which(dates > epoch_dates[i - 1] & dates <= epoch_dates[i])
      initial_time_from_parameters <- FALSE
    } else {
      dates1 <- which(dates > epoch_dates[i - 1])
      initial_time_from_parameters <- FALSE
    }
    
    if (length(dates1) == 0) {
      next
    }
    
    pars_model <- lapply(seq_rows(pars), function(j)
      transform(pars[j, ])[[i]]$pars)
    
    n_strains <- pars_model[[1]]$n_strains
    n_strains_R <- pars_model[[1]]$n_strains_R
    n_vacc_classes <- pars_model[[1]]$n_vacc_classes
    
    suffix <- paste0("_", c(ZamCovid:::model_age_bins()$start))
    S_nms <- get_names("S", list(n_vacc_classes), suffix)
    
    time1 <- time[dates1]
    S1 <- S[S_nms, , dates1, drop = FALSE]
    
    if (n_strains == 1) {
      R1 <- NULL
      prob_strain1 <- NULL
    } else {
      R_nms <- get_names("R",
                         list(S = n_strains_R, V = n_vacc_classes),
                         suffix)
      R1 <- R[R_nms, , dates1, drop = FALSE]
      prob_strain1 <- prob_strain[, , dates1, drop = FALSE]
    }
    
    rt1 <- sircovid::lancelot_Rt_trajectories(
      time1, S1, pars_model, type = type,
      initial_time_from_parameters = initial_time_from_parameters,
      shared_parameters = FALSE, R = R1, prob_strain = prob_strain1,
      weight_Rt = weight_Rt, keep_strains_Rt = keep_strains_Rt)
    
    reshape_rt <- function(r) {
      if (weight_Rt) {
        tmp <- array(NA, c(length(time1), 3, n_pars))
        tmp[, 3, ] <- r
      } else {
        tmp <- array(NA, c(length(time1), 2, n_pars))
      }
      tmp[, 1, ] <- r
      tmp
    }
    
    for (nm in names(rt)) {
      if (n_strains == 1 && !(weight_Rt && !keep_strains_Rt)) {
        if (nm %in% type) {
          rt1[[nm]] <- reshape_rt(rt1[[nm]])
        }
      }
      if (length(rt[[nm]]) == 0) {
        rt[[nm]] <- rt1[[nm]]
      } else if (length(dim(rt1[[nm]])) == 2) {
        rt[[nm]] <- rbind(rt[[nm]], rt1[[nm]])
      } else {
        rt[[nm]] <- abind1(rt[[nm]], rt1[[nm]])
      }
    }
  }
  
  if (weight_Rt && keep_strains_Rt) {
    for (nm in type) {
      colnames(rt[[nm]]) <- c("strain_1", "strain_2", "weighted")
    }
  }
  
  class(rt) <- c("Rt_trajectories", "Rt")
  rt
}


calculate_intrinsic_severity <- function(samples, base_pars) {
  
  pars <- samples$pars
  transform <- samples$predict$transform
  
  multiregion <- samples$info$multiregion
  
  what <- c("IFR", "IHR", "HFR")
  
  if (multiregion) {
    ## pars, state, time, info, epoch_dates
    ret <- lapply(samples$info$region, function(r)
      calculate_intrinsic_severity_region(
        pars[, , r], transform[[r]], what, base_pars[[r]]))
    names(ret) <- samples$info$region
  } else {
    ret <-
      calculate_intrinsic_severity_region(pars, transform, what, base_pars)
  }
  ret
}


calculate_intrinsic_severity_region <- function(pars, transform, what,
                                                base_pars) {
  
  dates <- base_pars$intrinsic_severity_dates
  strain_epochs <- base_pars$strain_epochs
  
  ## We calculate the intrinsic severity in pairs, so we split the strains
  ## into a list of pairs here
  strains <- split(strain_epochs, ceiling(seq_along(strain_epochs) / 2))
  
  time_vect <- dates * 4
  
  sev_vector <- function(x) {
    mean <- mean(x)
    lb <- quantile(x, 0.025)
    ub <- quantile(x, 0.975)
    
    out <- c(mean, lb, ub)
    names(out) <- c("mean", "lb", "ub")
    out
  }
  
  calc_instrinsic_severity_strains <- function(x) {
    # j will be the first stage in which the strains appear together
    j <- max(x) + 1
    pars_model <- lapply(spimalot:::seq_rows(pars),
                         function(i) transform(pars[i, ])[[j]]$pars)
    
    sev <- sircovid::lancelot_ifr_excl_immunity(time_vect, pars_model)
    sev$time <- NULL
    sev
  }
  
  ## Calculate the intrinsic severity for each of the pairs of strains
  intrinsic_severity_strains <- lapply(strains,
                                       calc_instrinsic_severity_strains)
  
  get_what <- function(w) {
    sev_variant <- function(x) {
      y <- t(apply(x, 1, sev_vector))
      data.frame(period = names(dates), y) %>%
        pivot_longer(!period, names_to = "estimate")
    }
    
    variants <- list()
    
    for (i in seq_along(strains)) {
      if (length(strains[[i]]) == 1) {
        ## If we have an odd number of strains, there will be one strain on
        ## its own, this deals with that case
        n_cols <- ncol(intrinsic_severity_strains[[i]][[w]])
        variants[[names(strains[[i]])[1]]] <-
          sev_variant(intrinsic_severity_strains[[i]][[w]][, n_cols, ])
      } else {
        variants[[names(strains[[i]])[1]]] <-
          sev_variant(intrinsic_severity_strains[[i]][[w]][, 1, ])
        variants[[names(strains[[i]])[2]]] <-
          sev_variant(intrinsic_severity_strains[[i]][[w]][, 2, ])
      }
    }
    
    dplyr::bind_rows(variants, .id = "name")
  }
  
  ret <- lapply(what, get_what)
  names(ret) <- what
  
  ret <- dplyr::bind_rows(ret, .id = "source")
  
  ret$period <- factor(ret$period, levels = unique(names(dates)))
  
  ret
}


extract_age_class_state <- function(state) {
  n_groups <- length(ZamCovid:::model_age_bins()$start)
  names_index <- c("cum_infections_disag", "diagnoses_admitted", "D_all")
  
  ## output cumulative states by
  ## age / vaccine class / sample / region / time
  arrays <- lapply(names_index, function(x) state[grep(x, rownames(state)), , ])
  names(arrays) <- c("infections", "diagnoses_admitted", "deaths")
  strata <- nrow(arrays$deaths) / n_groups
  
  f <- function(array) {
    x <- mcstate::array_reshape(array, 1L, c(n_groups, strata))
    vacc_classes <- c("unvaccinated", "partial_protection", "full_protection",
                      "waned_protection", "booster_protection",
                      "booster_waned_protection", "booster2_protection",
                      "booster2_waned_protection")
    colnames(x) <- vacc_classes[seq_len(ncol(x))]
    
    
    ## aggregate age groups
    groups <- list(age_0 = 1:6, # 0-4, 5-9, 10-14, 15-19, 20-24, 25-29
                   age_30 = 7:10,  # 30-34, 35-39, 40-44, 45-49
                   age_50 = 11:15, # 50-54, 55-59, 60-64, 65-69, 70-74
                   age_75 = 16:17, # 75-79, 80+
                   chw = 18, chr = 19)
    
    res <- lapply(groups,
                  function(i) apply(x[i, , , , drop = FALSE], 2:4, sum))
    
    # distribute CHW between 30-49 and 50-74 age groups
    # distribute CHR between 50-74 and 75+ age groups
    res$age_30 <- res$age_30 + 0.75 * res$chw
    res$age_50 <- res$age_50 + 0.25 * res$chw + 0.1 * res$chr
    res$age_75 <- res$age_75 + 0.9 * res$chr
    res$chw <- NULL
    res$chr <- NULL
    
    # take mean across particles
    ret <- apply(abind_quiet(res, along = 4), c(1, 3, 4), mean)
    
    # [age, vaccine status, region, time]
    ret <- round(aperm(ret, c(3, 1, 2)))
    ret <- mcstate::array_reshape(ret, 3, c(1, dim(ret)[3]))
    
    ret
  }
  
  lapply(arrays, f)
  
}

reduce_trajectories <- function(samples, severity) {
  ## Remove unused trajectories for predict function in combined
  remove_strings <- c("prob_strain", "S_", "R_", "I_weighted_", "D_hosp_",
                      "D_all_", "diagnoses_admitted_", "cum_infections_disag_",
                      "cum_n_vaccinated", "vacc_uptake_", "cum_admit_", "ifr",
                      "ihr", "hfr")
  
  if (severity) {
    remove_strings <- remove_strings[!remove_strings %in%
                                       c("diagnoses_admitted_", "vacc_uptake_",
                                         "D_hosp_")]
    
    samples <- get_deaths_admissions_by_vacc_class(samples)
    samples <- calc_weighted_protected(samples)
  }
  
  re <- sprintf("^(%s)", paste(remove_strings, collapse = "|"))
  
  samples <- summarise_states(samples)
  
  multiregion <- samples$info$multiregion
  state <- samples$trajectories$state
  keep <- !grepl(re, rownames(state))
  
  if (multiregion) {
    samples$trajectories$state <- samples$trajectories$state[keep, , , ]
  } else {
    samples$trajectories$state <- samples$trajectories$state[keep, , ]
  }
  
  ## TODO: this should be moved into a function with a more
  ## appropriate name
  
  ## Calculate Pillar 2 positivity and cases
  if (samples$info$model_type == "BB") {
    samples <- calculate_negatives(samples)
  } else {
    samples <- calculate_cases(samples)
  }
  
  samples
}


trajectories_filter_time <- function(trajectories, i) {
  trajectories$time <- trajectories$time[i]
  trajectories$date <- trajectories$date[i]
  trajectories$predicted <- trajectories$predicted[i]
  if (length(dim(trajectories$state)) == 3) {
    trajectories$state <- trajectories$state[, , i, drop = FALSE]
  } else {
    trajectories$state <- trajectories$state[, , , i, drop = FALSE]
  }
  trajectories
}


rt_filter_time <- function(rt, i, multiregion) {
  
  if (multiregion) {
    lapply(rt, rt_filter_time, i, FALSE)
  } else {
    rt_filter_time1 <- function(x) {
      if (length(dim(x)) == 2) {
        ret <- x[i, , drop = FALSE]
      } else {
        ret <- x[i, , , drop = FALSE]
      }
      ret
    }
    
    ret <- lapply(rt, rt_filter_time1)
    class(ret) <- class(rt)
    ret
  }
}


deaths_filter_time <- function(x, restart_date) {
  i <- x$date < restart_date
  for (v in c("output_t", "lower_bound", "upper_bound")) {
    x[[v]] <- x[[v]][i, , drop = FALSE]
  }
  x$date <- x$date[i]
  x$data <- x$data[ZamCovid:::numeric_date(x$data$date) < restart_date, ]
  x
}


calculate_vaccination <- function(state, vaccine_efficacy, cross_immunity) {
  
  de <- dim(vaccine_efficacy[[1]])
  if (length(de) == 2L) {
    ## We've changed the parameter preparation so that this branch
    ## should never be used
    message("WARNING: you are using old versions of tasks")
    n_groups <- de[[1]]
    n_strain <- 1
    n_vacc_classes <- de[[2]]
    multistrain <- FALSE
  } else {
    n_groups <- de[[1]]
    n_strain <- de[[2]]
    n_vacc_classes <- de[[3]]
    multistrain <- n_strain > 1
  }
  n_days <- dim(state)[[3]]
  
  # extract array of  mean by age / vaccine class / region == 1 / time
  get_mean_avt <- function(nm, state, strain = TRUE) {
    idx <- grep(nm, rownames(state))
    d <- c(n_groups, if (strain) n_strain else 1, n_vacc_classes)
    x <- mcstate::array_reshape(state[idx, , ], i = 1, d = d)
    out <- apply(x, c(1, 2, 3, 5), mean)
    ## TODO: this aperm should be removed, and code below here updated
    aperm(out, c(1, 3, 2, 4))
  }
  
  ## mean R by vaccine class / region == 1, time
  R <- get_mean_avt("^R_", state, TRUE)
  ## sum out age
  R <- apply(R, c(2, 3, 4), sum)
  
  ## mean cumulative vaccinations by age / vaccine class / region == 1 / time
  n_vaccinated <- get_mean_avt("^cum_n_vaccinated_", state, FALSE)
  
  vp <- get_vaccine_protection(vaccine_efficacy)
  
  # Methodology: calculate incidence of first / second doses,
  # number in each strata in total,
  # number in each strata who are recovered, use these to calculate proportion
  # protected as shown in main fig A / B.
  # to calculate proportion protected we need
  # s = stratum, r = region, t = time
  # let V be number in each vaccination stage V[a, s, r, t]
  # ler R be number recovered R[s, r, t]
  # from top to bottom of figure 1B
  # - number vaccinated sum_sr(V[s, r, t])
  # - number protected after vaccination against severe disease:
  #   i.e., sum_asr(V[a, s, r, t] * eff_severe[a, s])
  # - number protected after vaccination against infection
  #   i.e., sum_asr(V[a, s, r, t] * eff_inf[a, s])
  # - number protected after infection (this is added to all of the above):
  #   i.e., sum_sr(R[s, r, t])
  # - number protected after infection only:
  #   i.e. sum_r(R[1, r, t]
  
  
  ## create array of number in each vaccine stratum
  V <- array(0, dim = dim(n_vaccinated))
  V[, n_vacc_classes, , ] <- n_vaccinated[, n_vacc_classes - 1L, , ]
  for (i in seq(2, n_vacc_classes - 1)) {
    V[, i, , ] <- n_vaccinated[, i - 1, , ] - n_vaccinated[, i, , ]
  }
  
  sum_sr <- function(x) c(apply(x, c(2, 3), sum))
  sum_asr <- function(x) c(apply(x, c(3, 4), sum))
  
  if (multistrain) {
    
    ever_vaccinated <- colSums(n_vaccinated[, 1, , ])
    R_strain_1 <- apply(R[, -2, ], c(1, 3), sum) + R[, 2, ] * cross_immunity[2]
    R_strain_2 <-  apply(R[, -1, ], c(1, 3), sum) + R[, 1, ] * cross_immunity[1]
    
    n_protected <- list(
      strain_1 = rbind(
        ever_vaccinated = ever_vaccinated,
        protected_against_infection = sum_asr(c(vp$infection[, 1, ]) * V),
        protected_against_severe_disease =
          sum_asr(c(vp$severe_disease[, 1, ]) * V),
        protected_against_death = sum_asr(c(vp$death[, 1, ]) * V),
        ever_infected = colSums(R_strain_1),
        ever_infected_unvaccinated = R_strain_1[1, ]),
      strain_2 = rbind(
        ever_vaccinated = ever_vaccinated,
        protected_against_infection = sum_asr(c(vp$infection[, 2, ]) * V),
        protected_against_severe_disease =
          sum_asr(c(vp$severe_disease[, 2, ]) * V),
        protected_against_death = sum_asr(c(vp$death[, 2, ]) * V),
        ever_infected = colSums(R_strain_2),
        ever_infected_unvaccinated = R_strain_2[1, ]))
    
  } else {
    
    n_protected <- list(strain_1 = rbind(
      ever_vaccinated = colSums(n_vaccinated[, 1, , ]),
      protected_against_infection = sum_asr(c(vp$infection) * V),
      protected_against_severe_disease = sum_asr(c(vp$severe_disease) * V),
      protected_against_death = sum_asr(c(vp$death) * V),
      ever_infected = sum_sr(R),
      ever_infected_unvaccinated = R[1, , , drop = FALSE]))
  }
  
  
  ## calculate n_doses
  
  # Output number of first and second doses
  doses <- n_vaccinated[, c(1, 3), , , drop = FALSE]
  colnames(doses) <- c("first_dose", "second_dose")
  doses_inc <- aperm(apply(doses, c(1, 2, 3), diff), c(2, 3, 4, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-4], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")
  
  n_doses <- abind_quiet(doses, doses_inc, along = 2)
  
  list(n_protected = lapply(n_protected, mcstate::array_reshape, i = 2,
                            d = c(1, ncol(n_protected[[1]]))),
       n_doses = n_doses)
}


get_vaccine_protection <- function(vaccine_efficacy, booster_efficacy = NULL) {
  if (!is.null(booster_efficacy)) {
    stopifnot(identical(names(vaccine_efficacy), names(booster_efficacy)))
    vaccine_efficacy <- Map(cbind, vaccine_efficacy, booster_efficacy)
  }
  
  efficacy_infection <- 1 - vaccine_efficacy$rel_susceptibility
  efficacy_disease <- efficacy_infection + (1 - efficacy_infection) *
    (1 - vaccine_efficacy$rel_p_sympt)
  efficacy_severe_disease <- efficacy_disease + (1 - efficacy_disease) *
    (1 - vaccine_efficacy$rel_p_hosp_if_sympt)
  efficacy_death <- efficacy_severe_disease + (1 - efficacy_severe_disease) *
    (1 - vaccine_efficacy$rel_p_death)
  
  list(infection = efficacy_infection,
       disease = efficacy_disease,
       severe_disease = efficacy_severe_disease,
       death = efficacy_death)
}


extract_age_class_outputs <- function(samples) {
  ## Get index of states with age/vacc outputs
  names_index <- c("cum_infections_disag", "diagnoses_admitted", "D_all")
  index <- grep(paste(names_index, collapse = "|"),
                dimnames(samples$trajectories$state)[[1]])
  
  ## Save object for the age_vacc_class states
  if (length(dim(samples$trajectories$state)) == 4) {
    samples$trajectories$state[index, , , ]
  } else {
    samples$trajectories$state[index, , ]
  }
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


calculate_negatives <- function(samples) {
  date <- ZamCovid:::numeric_date_as_date(samples$trajectories$date)
  state <- samples$trajectories$state
  transform <- samples$predict$transform
  multiregion <- samples$info$multiregion
  
  if (multiregion) {
    region <- samples$info$region
    pars_model <- lapply(region, function(r)
      lapply(seq_len(nrow(samples$pars)), function(i)
        last(transform[[r]](samples$pars[i, , r]))$pars))
    negatives <- lapply(seq_along(region), function(i)
      calculate_negatives_region(state[, i, , ], pars_model[[i]], date))
    
    ## Add an extra dimension, then bind together:
    negatives <- lapply(negatives, mcstate::array_reshape,
                        2, c(1, dim(state)[[3]]))
    negatives <- abind_quiet(negatives, along = 2)
  } else {
    pars_model <- lapply(seq_len(nrow(samples$pars)), function(i)
      last(transform(samples$pars[i, ]))$pars)
    negatives <- calculate_negatives_region(state, pars_model, date)
  }
  
  ## Then bind into the main object:
  samples$trajectories$state <- abind_quiet(state, negatives, along = 1)
  
  samples
}


calculate_negatives_region <- function(state, pars_model, date) {
  
  pillar2_age_bands <- c("_under15", "_15_24", "_25_49",
                         "_50_64", "_65_79", "_80_plus")
  over25_age_bands <- c("_25_49", "_50_64", "_65_79", "_80_plus")
  p_NC_names <- c(paste0("p_NC", pillar2_age_bands),
                  paste0("p_NC_weekend", pillar2_age_bands))
  
  pars_base <- pars_model[[1]]
  pars <- t(vapply(pars_model, function(p) unlist(p[p_NC_names]),
                   numeric(length(p_NC_names))))
  
  calc_negs <- function(group) {
    neg1 <- pars_base[[paste0("N_tot", group)]] -
      state[paste0("sympt_cases", group, "_inc"), , ]
    
    neg1[, grepl("^S", weekdays(date))] <-
      neg1[, grepl("^S", weekdays(date))] *
      pars[, paste0("p_NC_weekend", group)]
    neg1[, !grepl("^S", weekdays(date))] <-
      neg1[, !grepl("^S", weekdays(date))] *
      pars[, paste0("p_NC", group)]
    
    neg1
  }
  
  ## Calculate the negatives for each age band
  neg <-
    vapply(pillar2_age_bands, calc_negs, array(0, dim = dim(state)[c(2, 3)]))
  neg <- aperm(neg, c(3, 1, 2))
  rownames(neg) <- pillar2_age_bands
  
  aggregate_age_bands <- function(x, name) {
    agg <- apply(x, c(2, 3), sum)
    agg <- array(agg, c(1, dim(agg)))
    rownames(agg) <- name
    agg
  }
  
  ## Calculate the negatives for aggregated age bands
  neg <- abind1(neg, aggregate_age_bands(neg[pillar2_age_bands, , ], ""))
  neg <- abind1(neg, aggregate_age_bands(neg[over25_age_bands, , ], "_over25"))
  rownames(neg) <- paste0("pillar2_negs", rownames(neg))
  
  neg
}


calculate_cases <- function(samples) {
  date <- ZamCovid:::numeric_date_as_date(samples$trajectories$date)
  state <- samples$trajectories$state
  transform <- samples$predict$transform
  multiregion <- samples$info$multiregion
  
  if (multiregion) {
    region <- samples$info$region
    pars_model <- lapply(region, function(r)
      lapply(seq_len(nrow(samples$pars)), function(i)
        last(transform[[r]](samples$pars[i, , r]))$pars))
    cases <- lapply(seq_along(region), function(i)
      calculate_cases_region(state[, i, , ], pars_model[[i]], date))
    
    ## Add an extra dimension, then bind together:
    cases <- lapply(cases, mcstate::array_reshape,
                    2, c(1, dim(state)[[3]]))
    cases <- abind_quiet(cases, along = 2)
  } else {
    pars_model <- lapply(seq_len(nrow(samples$pars)), function(i)
      last(transform(samples$pars[i, ]))$pars)
    cases <- calculate_cases_region(state, pars_model, date)
  }
  
  ## Then bind into the main object:
  samples$trajectories$state <- abind_quiet(state, cases, along = 1)
  
  samples
}


calculate_cases_region <- function(state, pars_model, date) {
  cases_names <-
    c("phi_pillar2_cases_under15", "phi_pillar2_cases_15_24",
      "phi_pillar2_cases_25_49", "phi_pillar2_cases_50_64",
      "phi_pillar2_cases_65_79", "phi_pillar2_cases_80_plus",
      "phi_pillar2_cases_weekend_under15", "phi_pillar2_cases_weekend_15_24",
      "phi_pillar2_cases_weekend_25_49", "phi_pillar2_cases_weekend_50_64",
      "phi_pillar2_cases_weekend_65_79", "phi_pillar2_cases_weekend_80_plus")
  
  
  pars <- t(vapply(pars_model, function(p) unlist(p[cases_names]),
                   numeric(length(cases_names))))
  
  calc_cases <- function(group) {
    cases1 <- state[paste0("sympt_cases_", group, "_inc"), , ]
    cases1[, grepl("^S", weekdays(date))] <-
      cases1[, grepl("^S", weekdays(date))] *
      pars[, paste0("phi_pillar2_cases_weekend_", group)]
    cases1[, !grepl("^S", weekdays(date))] <-
      cases1[, !grepl("^S", weekdays(date))] *
      pars[, paste0("phi_pillar2_cases_", group)]
    
    array(cases1, c(1, dim(cases1)))
  }
  
  cases_under15 <- calc_cases("under15")
  cases_15_24 <- calc_cases("15_24")
  cases_25_49 <- calc_cases("25_49")
  cases_50_64 <- calc_cases("50_64")
  cases_65_79 <- calc_cases("65_79")
  cases_80_plus <- calc_cases("80_plus")
  
  cases_over25 <- cases_25_49 + cases_50_64 + cases_65_79 + cases_80_plus
  cases_all <- cases_under15 + cases_15_24 + cases_over25
  
  cases <- abind1(cases_all, cases_over25)
  cases <- abind1(cases, cases_under15)
  cases <- abind1(cases, cases_15_24)
  cases <- abind1(cases, cases_25_49)
  cases <- abind1(cases, cases_50_64)
  cases <- abind1(cases, cases_65_79)
  cases <- abind1(cases, cases_80_plus)
  
  row.names(cases) <- paste0("pillar2_cases",
                             c("", "_over25", "_under15", "_15_24",
                               "_25_49", "_50_64", "_65_79", "_80_plus"))
  
  cases
}

summarise_states <- function(samples) {
  date <- ZamCovid:::numeric_date_as_date(samples$trajectories$date)
  state <- samples$trajectories$state
  transform <- samples$predict$transform
  multiregion <- samples$info$multiregion
  
  if (multiregion) {
    region <- samples$info$region
    pars_model <- lapply(region, function(r)
      lapply(seq_len(nrow(samples$pars)), function(i)
        last(transform[[r]](samples$pars[i, , r]))$pars))
    extra_state <- lapply(seq_along(region), function(i)
      summarise_states_region(state[, i, , ], pars_model[[i]]))
    
    ## Add an extra dimension, then bind together:
    extra_state <- lapply(extra_state, mcstate::array_reshape,
                          2, c(1, dim(state)[[3]]))
    extra_state <- abind_quiet(extra_state, along = 2)
  } else {
    pars_model <- lapply(seq_len(nrow(samples$pars)), function(i)
      last(transform(samples$pars[i, ]))$pars)
    extra_state <- summarise_states_region(state, pars_model)
  }
  
  ## Then bind into the main object:
  samples$trajectories$state <- abind_quiet(state, extra_state, along = 1)
  
  samples
}


summarise_states_region <- function(state, pars_model) {
  
  pars <- pars_model[[length(pars_model)]]
  n_groups <- pars$n_groups
  n_strains <- pars$n_strains
  n_strains_R <- pars$n_strains_R
  n_vacc_classes <- pars$n_vacc_classes
  
  i_recovered <- match(paste0("recovered", seq_len(n_strains_R)),
                       rownames(state))
  recovered <- state[i_recovered, , ]
  state <- state[-i_recovered, , ]
  recovered[is.na(recovered)] <- 0
  recovered[c(1, 2), , ] <- recovered[c(1, 2), , ] + recovered[c(4, 3), , ]
  recovered <- recovered[c(1, 2, 5), , ]
  row.names(recovered) <- c("recovered_1", "recovered_2", "recovered_historic")
  
  get_vaccine_status <- function(cum_n_vaccinated) {
    v <- c(pars$N_tot_all - cum_n_vaccinated[1],
           -diff(cum_n_vaccinated))
    v
  }
  
  vaccine_status <- state[grep("^cum_n_vaccinated_", rownames(state)), , ]
  
  vacc_uptake <- which(endsWith(rownames(vaccine_status), "_1"))[c(1:17)]
  vacc_uptake <- vaccine_status[vacc_uptake, , ]
  rownames(vacc_uptake) <- paste0("vacc_uptake_", seq_len(nrow(vacc_uptake)))
  extra_states <- abind1(recovered, vacc_uptake)
  
  vaccine_status[is.na(vaccine_status)] <- 0
  vaccine_status <- array(vaccine_status, c(n_groups, n_vacc_classes,
                                            dim(vaccine_status)[2:3]))
  vaccine_status <- apply(vaccine_status, c(2, 3, 4), sum)
  vaccine_status <- apply(vaccine_status, c(2, 3), get_vaccine_status)
  row.names(vaccine_status) <-
    paste0("vaccine_status_", seq_len(n_vacc_classes))
  extra_states <- abind1(extra_states, vaccine_status)
  
  extra_states
}

## adapted from sircovid:::calculate_index
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

extract_demography <- function(samples) {
  
  multiregion <- samples$info$multiregion
  
  if (multiregion) {
    extract1 <- function(i) {
      samples$trajectories$state <- samples$trajectories$state[, i, , ]
      samples$predict$transform <- samples$predict$transform[[i]]
      samples$pars <- samples$pars_full[, , i]
      
      extract_demography_region(samples)
    }
    
    regions <- samples$info$region
    ret <- lapply(seq_along(regions), extract1)
    names(ret) <- regions
    
  } else {
    ret <- extract_demography_region(samples)
  }
  
  ret
}

extract_demography_region <- function(samples) {
  # remove the first date as it is more than one day before the second date
  trajectories <- samples$trajectories$state[, , -1L]
  # remove the second date here also as we will lose a date doing diff below
  date <- samples$trajectories$date[-c(1, 2)]
  
  admissions <- trajectories[grep("^cum_admit", rownames(trajectories)), , ]
  deaths_hosp <- trajectories[grep("^D_hosp_", rownames(trajectories)), , ]
  D_all <- trajectories[grep("^D_all_", rownames(trajectories)), , ]
  
  ## We obtain community deaths by subtracting deaths_hosp from D_all
  deaths_comm <- D_all - deaths_hosp
  rownames(deaths_comm) <- gsub("all", "comm", rownames(deaths_comm))
  
  process <- function(x) {
    total <- x[, , dim(x)[3]]
    prop_total <- t(total) / colSums(total) * 100
    
    daily <- apply(x, c(1, 2), diff)
    mean <- apply(daily, c(1, 2), mean)
    
    list(total = total,
         prop_total = prop_total,
         mean_prop_total = colMeans(prop_total),
         mean_t = mean,
         date = date)
  }
  
  output <- list(
    admissions = process(admissions),
    deaths_hosp = process(deaths_hosp),
    deaths_comm = process(deaths_comm))
}


extract_severity <- function(samples) {
  time <- samples$trajectories$time
  date <- samples$trajectories$date
  state <- samples$trajectories$state
  
  multiregion <- samples$info$multiregion
  
  if (multiregion) {
    ## state, time, date
    ret <- lapply(samples$info$region, function(r)
      extract_severity_region(state[, r, , ], time, date))
    names(ret) <- samples$info$region
  } else {
    ret <- extract_severity_region(state, time, date)
  }
  ret
}


get_r_vaccinated <- function(samples) {
  
  multiregion <- samples$info$multiregion
  state <- samples$trajectories$state
  dim <- dim(state)
  info <- samples$info$info
  dim_R <- info[[length(info)]]$dim$R
  
  R_names <- grep("^R", rownames(state), value = TRUE)
  R <- array(state[R_names, , ], dim = c(dim_R, dim[-1]))
  R[is.na(R)] <- 0
  
  
  ## Sum over ages and strains (first two dimensions)
  ## Then drop unvaccinated
  if (multiregion) {
    recovered_V <- apply(R, c(3, 4, 5, 6), sum)
    recovered_V <- recovered_V[-1, , , ]
  } else {
    recovered_V <- apply(R, c(3, 4, 5), sum)
    recovered_V <- recovered_V[-1, , ]
  }
  
  rownames(recovered_V) <- paste0("recovered_V", seq_len(nrow(recovered_V)))
  
  
  ## Get total recovered vaccinated
  if (multiregion) {
    recovered_V_all <- apply(recovered_V, c(2, 3, 4), sum)
  } else {
    recovered_V_all <- apply(recovered_V, c(2, 3), sum)
  }
  recovered_V_all <- array(recovered_V_all, c(1, dim(recovered_V_all)))
  rownames(recovered_V_all) <- "recovered_V_all"
  
  state <- abind1(state, recovered_V_all)
  state <- abind1(state, recovered_V)
  samples$trajectories$state <- state
  
  samples
}

get_deaths_admissions_by_vacc_class <- function(samples) {
  
  multiregion <- samples$info$multiregion
  state <- samples$trajectories$state
  dim <- dim(state)
  info <- samples$info$info
  
  disag_names <- c("diagnoses_admitted", "D_hosp")
  for (i in disag_names) {
    
    dim_i <- info[[length(info)]]$dim[[i]]
    
    names_i <- grep(paste0("^", i), rownames(state), value = TRUE)
    disag_i <- mcstate::array_reshape(state[names_i, , ], 1, dim_i)
    disag_i[is.na(disag_i)] <- 0
    
    if (multiregion) {
      vacc_i <- apply(disag_i, c(2, 3, 4, 5), sum)
      state <- state[!rownames(state) %in% names_i, , , ]
    } else {
      vacc_i <- apply(disag_i, c(2, 3, 4), sum)
      state <- state[!rownames(state) %in% names_i, , ]
    }
    
    rownames(vacc_i) <- paste(i, "vacc", seq(0, nrow(vacc_i) - 1), sep = "_")
    state <- abind1(state, vacc_i)
  }
  
  samples$trajectories$state <- state
  samples
}


calc_weighted_protected <- function(samples) {
  
  state <- samples$trajectories$state
  multiregion <- samples$info$multiregion
  
  protected_names <- c("protected_S_vaccinated",
                       "protected_R_unvaccinated",
                       "protected_R_vaccinated")
  
  if (multiregion) {
    weight <- state[c("prob_strain", "prob_strain_1"), , , ]
    weight[is.na(weight)] <- 0
    
    wt_protected <- function(nm) {
      protected <- state[paste(nm, c(1, 2), sep = "_"), , , ]
      protected[is.na(protected)] <- 0
      
      apply(protected * weight, c(2, 3, 4), sum)
    }
    
    weighted_state <- vapply(protected_names, wt_protected,
                             array(0, dim(state)[c(2, 3, 4)]))
    weighted_state <- aperm(weighted_state, c(4, 1, 2, 3))
    rownames(weighted_state) <- paste0(protected_names, "_weighted")
  } else {
    weight <- state[c("prob_strain", "prob_strain_1"), , ]
    weight[is.na(weight)] <- 0
    
    wt_protected <- function(nm) {
      protected <- state[paste(nm, c(1, 2), sep = "_"), , ]
      protected[is.na(protected)] <- 0
      
      apply(protected * weight, c(2, 3), sum)
    }
    
    weighted_state <- vapply(protected_names, wt_protected,
                             array(0, dim(state)[c(2, 3)]))
    weighted_state <- aperm(weighted_state, c(3, 1, 2))
    rownames(weighted_state) <- paste0(protected_names, "_weighted")
  }
  
  samples$trajectories$state <- abind1(state, weighted_state)
  samples
  
}


extract_severity_region <- function(state, time, date) {
  
  # String vectors to formulate severity trajectory names needed
  sev_traj <- grep("^ifr|^ihr|^hfr", rownames(state), value = TRUE)
  
  # ignore CHR and CHW
  sev_traj <- grep("CHR|CHW", sev_traj, invert = TRUE, value = TRUE)
  
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