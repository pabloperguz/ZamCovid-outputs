
create_baseline <- function(region, date, epoch_dates, pars, assumptions) {
  
  pars_info <- pars$info
  pars_vcv <- pars$proposal
  
  message(sprintf("Creating baseline for '%s'", region))
  pars_info <- pars_info[pars_info$region == region | is.na(pars_info$region), ]
  pars_info <- setNames(pars_info$initial, pars_info$name)
  
  ## 1. Set up basic variables ----
  date <- as.Date(date)
  dt <- 0.25
  vaccine_days_to_effect_dose1 <- 0
  vaccine_days_to_effect_dose2 <- 7
  vaccine_days_to_effect <- c(vaccine_days_to_effect_dose1,
                              vaccine_days_to_effect_dose2)
  
  ## 2. Read in data ----
  uptake_by_age <- read.csv("data/vaccine_uptake.csv", row.names = "group")
  data_vaccination <- read_csv("data/data_vaccination.csv")
  vaccine_efficacy <- read_csv("data/vaccine_efficacy.csv")
  progression_data <- read_csv("data/progression_data.csv")
  severity_data <- read_csv("data/support_severity.csv")
  population <- read_csv("data/population.csv")[, c("age", region)] %>%
    `colnames<-`(., c("age", "n"))
  
  if (region == "kabwe") {
    # This is a hack!!!
    # We don't have population breakdown from Kabwe!
    # n corresponds to national population from 2018, at a point the population
    # in Kabwe District was reported to be 230,802
    
    # New 2022 census indicates a national population of 19,610,769
    # We'll assume the same growth for Kabwe and same age distribution
    n_zambia <- 19610769
    g <-  n_zambia / sum(population$n)
    n_kabwe <- ceiling(230802 * g)
    population$n <- round((population$n / sum(population$n)) * n_kabwe)
  }
  
  ## 3. Set-up basic model parameters ----
  # Beta change points - A vector of date (strings) for the beta parameters.
  # We are agnostic as to the effect of official NPI dates at specific dates
  # This is thus just a vector of monthly dates
  beta_date <- as.character(
    seq.Date(as.Date("2020-03-01"), as.Date(date), "month"))
  beta_names <- sprintf("beta%d", seq_along(beta_date))
  
  # Set of parameters that will be fitted for each model type
  to_fit_all <- c(
    # direct
    "start_date", beta_names,
    # severity
    "p_G_D"
  )
  stopifnot(setequal(to_fit_all, names(pars_info)))
  
  
  ## Initial seeding: seed 10 per million over 1 day (4 steps)
  seed_size <- 10e-6 * sum(population[, "n"])
  seed_pattern <- c(1, 1, 1, 1)
  
  
  # For internal check: what is the assumed serial interval of SARS-CoV-2?
  # Calculation adapted from Svensson
  # https://www2.math.su.se/matstat/reports/seriea/2005/rep14/report.pdf
  # so SI = E + (P^2 + P·C_1 + C_1^2) / (P + C_1)
  # We ensure that mean[T_I_C_1 + T_I_C_2] is unchanged
  rel_si_wildtype <- 1
  mean_E <- 2 / 0.865
  mean_P <- 1.68
  mean_C_1 <- 2.14
  mean_C_2 <- 1.86
  rel_C_2_wildtype <- (mean_C_2 + (1 - rel_si_wildtype) * mean_C_1) / mean_C_2
  mean_SI <- mean_E + (mean_P^2 + mean_C_1 + mean_C_2^2) / (mean_P + mean_C_1)
  rel_gamma_wildtype <- list(E = 1 / rel_si_wildtype,
                             A = 1 / rel_si_wildtype,
                             P = 1 / rel_si_wildtype,
                             C_1 = 1 / rel_si_wildtype,
                             C_2 = 1 / rel_C_2_wildtype)
  
  
  ## 5. Set-up vaccination parameters and assumptions ----
  vaccine_eligibility_min_age <- 5
  mean_days_between_doses <- 7 * 11 # second dose starts 12 weeks after first
  mean_time_to_waned <- 24 * 7 # assume exponential with mean 24 weeks
  time_to_dose_1_effect <- 7 * 3 # assume exponential with mean 3 weeks
  
  # Compartment progression rates
  # These are progression rates OUT of the specified compartment, where
  # zeroes mean movement is controlled by vaccination schedule rather than
  # a rate parameter and waned is an absorbing state
  vaccine_progression_rate <- c(0,                         # unvaccinated 
                                1 / time_to_dose_1_effect, # first dose no effect
                                0,                         # first dose full effect
                                1/ mean_time_to_waned,     # second dose
                                0)                         # waned

  # Proportion of vaccines by type and age
  vacc_doses_by_age_date_vaccine <-
    data_vaccination %>%
    dplyr::mutate(age_band_min = replace(age_band_min, age_band_min == 16, 15)) %>%
    dplyr::filter(!is.na(age_band_min),
                  age_band_min != 80) %>%
    dplyr::group_by(vaccine, age_band_min) %>%
    dplyr::summarise(n = sum(first_dose, na.rm = TRUE)) %>%
    dplyr::group_by(age_band_min) %>%
    dplyr::mutate(freq = n / sum(n, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = c(age_band_min),
                       values_from = freq, names_from = vaccine)
  
  # check all proportions sum to 1
  stopifnot(
    all(abs(rowSums(
      vacc_doses_by_age_date_vaccine[3:16, -1], na.rm = TRUE) - 1) < 1e-8))
  
  # check proportion that are Pfizer/Moderna by age group
  prop_pfizer <- apply(cbind(vacc_doses_by_age_date_vaccine$Pfizer,
                             vacc_doses_by_age_date_vaccine$Moderna), 1, sum, na.rm = TRUE)
  
  # if vaccinated under the age of 10 should be Pfizer/Moderna
  prop_pfizer[ZamCovid:::model_age_bins()$end < 10] <- 1
  
  
  ## Calculate average vaccine efficacy and rel_severity
  average_vacc_efficacy <-
    calculate_average_vacc_efficacy(vaccine_efficacy, prop_pfizer)
  rel_severity <- lapply(average_vacc_efficacy, function(e)
    get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                 e$infection, e$transmission))
  
  ## Set right VE assumption values across sensitivity analyses
  if (assumptions == "ve_high") {
    rel_severity <- rel_severity$ve_high
  } else if (assumptions == "ve_low") {
    rel_severity <- rel_severity$ve_low
  } else {
    rel_severity <- rel_severity$central
  }
  
  # If required for SA, change Se&Sp parameters here
  sens_and_spec <- ZamCovid::ZamCovid_parameters_sens_and_spec()
  
  ## Note that vaccine_uptake[i, j] is proportional uptake of dose j for group i 
  vaccine_uptake <- 
    array(uptake_by_age$central, c(length(uptake_by_age$central), 2))
  
  data_vaccination <- data_vaccination %>%
    dplyr::rename(dose1 = first_dose,
                  dose2 = second_dose,
                  dose3 = third_dose,
                  dose4 = fourth_dose,
                  dose5 = fifth_dose)
  n_doses <- 3
  dose_start_dates <- c("2020-12-08",
                        "2020-12-08",
                        "2021-09-15")
  
  # vaccination <- 
  #   vaccination_schedule(date, region, vaccine_uptake, 
  #                        vaccine_days_to_effect, data_vaccination,
  #                        n_doses, dose_start_dates)
  # 
  # vaccine_schedule_real <- vaccination$schedule
  
  ## shift doses to account for time between dose and effect of dose
  ## TODO: WIP
  vaccine_schedule_effect <- NULL
    # shift_doses(vaccine_schedule_real, vaccine_days_to_effect)
  
  n_doses <- NULL # vaccination$schedule$n_doses
  vaccine_index_dose2 <- 3L
  
  ret <- list(
    date = date,
    region = region,
    population = population,
    epoch_dates = ZamCovid:::numeric_date(epoch_dates),
    
    beta_date = ZamCovid:::numeric_date(beta_date),
    beta_names = beta_names,
    severity_data = severity_data,
    progression_data = progression_data,
    sens_and_spec = sens_and_spec,
    seed_size = seed_size,
    seed_pattern = seed_pattern,
    
    rel_severity = rel_severity,
    rel_gamma_wildtype = rel_gamma_wildtype,
    
    vaccine_eligibility_min_age = vaccine_eligibility_min_age,
    vaccine_progression_rate = vaccine_progression_rate,
    vaccine_schedule = NULL, #vaccine_schedule_real,
    vaccine_schedule_effect = vaccine_schedule_effect,
    vaccine_uptake = vaccine_uptake,
    vaccine_mean_days_between_doses = mean_days_between_doses,
    vaccine_index_dose2 = vaccine_index_dose2,
    vaccine_days_to_effect = vaccine_days_to_effect,
    n_doses = n_doses)
  
  message("  - Creating transformation function")
  tr <- make_transform(ret)
  message("  - Testing transformation function")
  p <- tr(pars_info)
  message("  - Testing creating model with transformed parameters")
  m <- ZamCovid::ZamCovid$new(p, 0, 1)
  ret
}