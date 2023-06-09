shift_doses <- function(vaccine_schedule, vaccine_days_to_effect) {
  ## This function shifts, for each dose j, all jth doses in vaccine_schedule by
  ## vaccine_days_to_effect[j]
  n_groups <- dim(vaccine_schedule$doses)[1]
  n_doses <- dim(vaccine_schedule$doses)[2]
  n_days <- dim(vaccine_schedule$doses)[3]
  
  shift1 <- function(i, j) {
    c(rep(0, vaccine_days_to_effect[j]),
      vaccine_schedule$doses[i,j,],
      rep(0, max(vaccine_days_to_effect) - vaccine_days_to_effect[j]))
  }
  
  schedule_effect <-
    vapply(seq_len(n_groups),
           function (i) {
             vapply(seq_len(n_doses), function(j) shift1(i, j),
                    numeric(n_days + max(vaccine_days_to_effect)))},
           array(0, c(n_days + max(vaccine_days_to_effect), n_doses)))
  schedule_effect <- aperm(schedule_effect, c(3, 2, 1))
  
  vaccine_schedule$doses <- schedule_effect
  vaccine_schedule
}

calculate_average_vacc_efficacy <- function(vaccine_efficacy, prop_pfizer) {
  
  # unvax / 1st dose no eff / 1st dose full eff / 2nd dose full eff / waned
  vaccine_efficacy$week_wane <- NULL
  ret <- vaccine_efficacy %>%
    tidyr::pivot_longer(-c(type, vaccine, dose), names_to = "analysis") %>%
    tidyr::expand_grid(group = seq_len(16)) %>%
    # code dose numbers according to their vaccine strata
    dplyr::mutate(prop_pfizer = prop_pfizer[group],
                  stratum = forcats::fct_recode(as.character(dose),
                                                "stratum_3" = "1",
                                                "stratum_4" = "2",
                                                "stratum_5" = "waned")) %>%
    dplyr::select(-dose) %>%
    tidyr::pivot_wider(names_from = stratum) %>%
    # add in strata 1 in which there is 0 vaccine protection
    dplyr::mutate(stratum_1 = 0, .before = stratum_3) %>%
    # add in strata 3 in which there is 0 vaccine protection for 3 weeks after dose 1
    dplyr::mutate(stratum_2 = 0, .before = stratum_3) %>%
    tidyr::pivot_longer(dplyr::starts_with("stratum"), names_to = "stratum") %>%
    tidyr::pivot_wider(names_from = vaccine) %>%
    # calculate combined vaccine efficacy based on % of PF
    dplyr::mutate(efficacy = AZ * (1 - prop_pfizer) + PF * prop_pfizer) %>%
    dplyr::select(-c(AZ, PF)) %>%
    tidyr::pivot_wider(names_from = type, values_from = efficacy) %>%
    dplyr::arrange( analysis, stratum, group)
  
  ## split by analysis type
  split(ret, ret$analysis)
}

get_vaccine_conditional_prob <- function(eff_death,
                                         eff_severe_disease,
                                         eff_disease, eff_infection,
                                         eff_onwards_transmission = NULL) {
  
  n_group <- 16
  rel_susceptibility <- matrix(1 - eff_infection, n_group)
  rel_p_sympt <- matrix(1 - eff_disease, n_group) / rel_susceptibility
  rel_p_hosp_if_sympt <-
    matrix(1 - eff_severe_disease, n_group) / (rel_susceptibility * rel_p_sympt)
  rel_p_death <-
    matrix(1 - eff_death, n_group) / (rel_susceptibility * rel_p_sympt * rel_p_hosp_if_sympt)
  res <- list(rel_susceptibility = rel_susceptibility,
              rel_p_sympt = rel_p_sympt,
              rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
              rel_p_death = rel_p_death)
  
  if(!is.null(eff_onwards_transmission)) {
    res$rel_infectivity <- matrix(1 - eff_onwards_transmission, n_group)
  } else {
    res$rel_infectivity <- res$rel_susceptibility
    res$rel_infectivity[] <- 1
  }
  
  ## This caused problems because of rounding.
  ## never interested in efficacy beyond 0.1% hence the rounding
  if (any(signif(unlist(res), 3) < 0) |
      any(signif(unlist(res), 3) > 1)) {
    stop("incompatible vaccine efficacy parameters")
  }
  
  res
}


vaccination_schedule <- function(date, region, uptake, days_to_effect, data,
                                 n_doses, dose_start_dates) {
  
  dose_cols <- paste0("dose", seq_len(n_doses))
  if (!all(dose_cols %in% names(data))) {
    stop(sprintf("n_doses = %s so expected dose column names: %s",
                 n_doses, paste(squote(dose_cols), collapse = ", ")))
  }
  if (length(dose_start_dates) != n_doses) {
    stop(sprintf("n_doses = %s so expected length of dose_start_dates to be %s",
                 n_doses, n_doses))
  }
  
  ## Remove dose columns above n_doses
  keep <- c("date", "region", "vaccine", "age_band_min", "age_band_max",
            dose_cols)
  data <- data[, keep]
  
  ## Remove all data before earliest dose start date
  data <- data[data$date >= min(dose_start_dates), ]
  
  ## Ignore dose data before each corresponding start date
  for (i in seq_len(n_doses)) {
    data[[paste0("dose", i)]][data$date < dose_start_dates[i]] <- 0
  }
  
  ## Remove any data after date parameter
  data <- data[data$date <= date, ]
  
  ## Make sure there are no funny region names
  data$region <- tolower(data$region)
  data$region <- gsub(" ", "_", data$region)
  data <- data[data$region == region, ]
  
  ## Summing over vaccine types
  data <- data %>%
    tidyr::pivot_longer(dose_cols, names_to = "dose", values_to = "count") %>%
    dplyr::group_by(.data$date, .data$region,
                    .data$age_band_min, .data$age_band_max, .data$dose) %>%
    dplyr::summarise(count = sum(.data$count, na.rm = TRUE))
  
  ## Fill in any missing dates
  missing_dates <-
    setdiff(as.character(seq.Date(as.Date(min(data$date)), date, by = 1)),
            unique(data$date))
  if (length(missing_dates) > 0) {
    data_missing <- data.frame(date = missing_dates,
                               region = region,
                               age_band_min = NA,
                               age_band_max = NA,
                               dose = "dose1",
                               count = 0)
    data <- rbind(data, data_missing)
  }
  data <- data %>%
    dplyr::arrange(date, .data$age_band_min)
  
  ## remove any trailing days with zero doses if they are within
  ## min(days_to_effect) of the date parameter
  agg_data <- data %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(count = sum(count))
  dates_remove <-
    agg_data$date[agg_data$date > max(agg_data$date[agg_data$count > 0])]
  dates_remove <-
    dates_remove[dates_remove > date - min(days_to_effect)]
  data <- data[!(data$date %in% dates_remove), ]
  
  data <- data %>%
    tidyr::pivot_wider(names_from = dose, values_from = count)
  
  # This and related functions are WIP
  schedule <- NULL # ZamCovid::vaccine_schedule_from_data(data, region, uptake)
  
  ret <- list(
    data = data,
    schedule = schedule)
  class(ret) <- "vaccination_data" # soon
  ret
}


shift_doses <- function(vaccine_schedule, vaccine_days_to_effect) {
  ## This function shifts, for each dose j, all jth doses in vaccine_schedule
  # by vaccine_days_to_effect[j]
  n_groups <- dim(vaccine_schedule$doses)[1]
  n_doses <- dim(vaccine_schedule$doses)[2]
  n_days <- dim(vaccine_schedule$doses)[3]
  
  shift1 <- function(i, j) {
    c(rep(0, vaccine_days_to_effect[j]),
      vaccine_schedule$doses[i,j,],
      rep(0, max(vaccine_days_to_effect) - vaccine_days_to_effect[j]))
  }
  
  schedule_effect <-
    vapply(seq_len(n_groups),
           function (i) {
             vapply(seq_len(n_doses), function(j) shift1(i, j),
                    numeric(n_days + max(vaccine_days_to_effect)))},
           array(0, c(n_days + max(vaccine_days_to_effect), n_doses)))
  schedule_effect <- aperm(schedule_effect, c(3, 2, 1))
  
  vaccine_schedule$doses <- schedule_effect
  vaccine_schedule
}
