compute_severity <- function(pars, severity_data) {
  dt <- 0.25 # TODO: tidy this up
  ## WARNING: this is a hack
  list2env(pars, environment())
  
  # Overall probability of death (mechanistically of dying outside hospital)
  p_G_D_value <- p_G_D
  
  # Simplifying assumptions
  p_star_value <- 0.1 # Probability of being admitted with positive PCR
  p_H_value <- 0 # Probability of hospitalisation
  p_H_D_value <- 0 # Probability death in hospital
  
  # e. Simplifying assumption - note sero_sensitivity is 0.831
  severity_data[severity_data$Name == "p_sero_pos", 2:17] <- 1
  
  severity <- ZamCovid::ZamCovid_parameters_severity(
    dt,
    severity_data,
    p_G_D = list(value = p_G_D_value),
    p_H = list(value = p_H_value),
    p_H_D = list(value = p_H_D_value),
    p_star = list(value = p_star_value))
  
  severity
}


compute_progression <- function(pars, progression_data) {
  dt <- 0.25
  ## WARNING: this is a hack
  list2env(pars, environment())
  
  k_parameters <-
    progression_data[grep("^k_", progression_data$parameter), ]
  gammas <-
    progression_data[grep("^gamma_", progression_data$parameter),
                     "value"]
  names(gammas) <-
    progression_data[grep("^gamma_", progression_data$parameter),
                     "parameter"]
  gammas <- as.list(gammas)
  
  gamma_E <- gammas$gamma_E
  gamma_H_R_value <- gammas$gamma_H_R
  gamma_H_D_value <- gammas$gamma_H_D
  gamma_PCR_pre_value <- 1 / 2
  gamma_PCR_pos_value <- 0.083
  
  # Time to diagnosis if admitted without test
  gamma_U_value <- 1 / 3
  
  progression <- ZamCovid::ZamCovid_parameters_progression(
    dt,
    gamma_E = list(value = gamma_E),
    gamma_H_D = list(value = gamma_H_D_value),
    gamma_H_R = list(value = gamma_H_R_value),
    gamma_U = list(value = gamma_U_value),
    gamma_PCR_pre = list(value = gamma_PCR_pre_value),
    gamma_PCR_pos = list(value = gamma_PCR_pos_value)
  )
  progression[k_parameters$parameter] <- k_parameters$value
  
  progression$gamma_R <- gammas$gamma_R
  progression$k_sero_pre <- 1
  progression$gamma_sero_pre <- 1 / 13
  progression$k_PCR_pre <- 1
  progression$k_PCR_pos <- 1
  progression$k_sero_pos <- 1
  progression$gamma_sero_pos <- gammas$gamma_sero_pos
  
  progression
}


compute_observation <- function(pars) {
  
  ## WARNING: this is a hack
  list2env(pars, environment())
  
  observation <- ZamCovid::ZamCovid_parameters_observation()
  observation$kappa_death_all <- 1 / alpha_D
  
  observation
}


make_transform <- function(baseline) {
  
  expected <- c("date", "region", "population", "epoch_dates",
                "beta_date", "beta_names",
                "severity_data", "progression_data",
                "sens_and_spec", "seed_size", "seed_pattern",
                "rel_gamma_wildtype", "base_death_date", "base_death_value",
                ## Lots of vaccination things
                "rel_severity", "vaccine_eligibility_min_age",
                "vaccine_progression_rate",
                "vaccine_schedule", "vaccine_schedule_effect",
                "vaccine_uptake", "vaccine_mean_days_between_doses",
                "vaccine_index_dose2", "vaccine_days_to_effect",
                "n_doses")
  stopifnot(setequal(expected, names(baseline)))
  
  epoch_dates <- baseline$epoch_dates
  
  
  expected <- c("start_date", baseline$beta_names, "p_G_D", "alpha_D")
  
  function(pars) {
    
    stopifnot(setequal(expected, names(pars)))
    beta_value <- unname(pars[baseline$beta_names])
    pars <- as.list(pars) # using list access below
    
    progression <- compute_progression(pars, baseline$progression_data)
    severity <- compute_severity(pars, baseline$severity_data)
    observation <- compute_observation(pars)
    
    vaccine_schedule <- baseline$vaccine_schedule_effect
    
    stage_parameters <- function(vaccine_doses) {
      
      rel_severity <- baseline$rel_severity
      
      if (vaccine_doses == 0) {
        
        rel_severity <- lapply(
          rel_severity, function(x) x[, 1, drop = FALSE])
        vaccine_progression_rate <- 0
        vaccine_schedule <- NULL
        vaccine_index_dose2 <- NULL
        ## this is just the default value
        n_doses <- 2L

      } else if (vaccine_doses == 2) {

        rel_severity <- lapply(
          rel_severity, function(x) x[, , 1:5, drop = FALSE])
        vaccine_progression_rate <- baseline$vaccine_progression_rate[1:5]
        vaccine_schedule$doses <- vaccine_schedule$doses[, 1:2, , drop = FALSE]
        vaccine_schedule$n_doses <- 2L
        vaccine_index_dose2 <- baseline$vaccine_index_dose2
        n_doses <- 2L

      } else {
        stop("vaccine_doses must be 0 or 2")
      }

      ZamCovid::ZamCovid_parameters(
        start_date = pars$start_date,
        # region = baseline$region,
        population = baseline$population,
        beta_date = baseline$beta_date,
        beta_value = beta_value,
        
        base_death_date = baseline$base_death_date,
        base_death_value = baseline$base_death_value,
        
        severity = severity,
        progression = progression,
        observation = observation,
        sens_and_spec = baseline$sens_and_spec,
        
        initial_seed_size = baseline$seed_size,
        initial_seed_pattern = baseline$seed_pattern,
        
        rel_susceptibility = rel_severity$rel_susceptibility,
        rel_p_sympt = rel_severity$rel_p_sympt,
        rel_p_hosp_if_sympt = rel_severity$rel_p_hosp_if_sympt,
        rel_p_death = rel_severity$rel_p_death,
        rel_infectivity = rel_severity$rel_infectivity,
        
        vaccine_progression_rate = vaccine_progression_rate,
        vaccine_schedule = vaccine_schedule,
        vaccine_index_dose2 = vaccine_index_dose2,
        n_doses = n_doses
      )
      
    }
    stage_parameters(0)
  }
}
