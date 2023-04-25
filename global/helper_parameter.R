setwd(orderly::orderly_config()$root)
source("global/util.R")

add_parameter <- function(name, initial, min, max, proposal,
                          assumptions = "central", model = "deterministic",
                          integer = FALSE, include = TRUE) {
  
  info_filename <-
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename) 
  
  regions <- unique(parameters_info$region)
  regions <- regions[!is.na(regions)]
  
  new_par <- data.frame(
    region = regions,
    name = name,
    initial = initial,
    min = min,
    max = max,
    integer = integer,
    include = include
  )
  
  parameters_info <- rbind(parameters_info, new_par)
  parameters_info <- dplyr::arrange(parameters_info,
                                    region, name)

  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  parameters_proposal <- read_csv(proposal_filename)
  
  parameters_proposal[[name]] <- 0
  k <- parameters_proposal$name == parameters_proposal$name[1]
  
  new_prop <- parameters_proposal[k, ]
  new_prop$name <- name
  new_prop[, -c(1, 2)] <- 0
  
  if (!is.numeric(proposal)) {
    avg_prop <- NULL
    
    for (i in regions) {
      
      avg_prop[i] <- parameters_proposal[
        parameters_proposal$region == i &
          parameters_proposal$name == proposal, proposal]
    }
    proposal = abs(mean(avg_prop))
    stopifnot(is.numeric(proposal))
  }
  new_prop[[name]] <- proposal
  
  parameters_proposal <- rbind(parameters_proposal, new_prop)
  parameters_proposal <- dplyr::arrange(parameters_proposal, region, name)
  col_order <- c(1, 2, 2 + order(names(parameters_proposal)[-c(1,2)]))
  parameters_proposal <- parameters_proposal[, col_order]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}


remove_parameter <- function(name, assumptions = "central", model = "deterministic") {
  
  info_filename <-
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste("src/ZamCovid_parameters/pars",assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename)
  parameters_info <- parameters_info[parameters_info$name != name, ]
  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  
  parameters_proposal <- read_csv(proposal_filename)
  parameters_proposal[[name]] <- NULL
  parameters_proposal <- parameters_proposal[parameters_proposal$name != name, ]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
  
}


rename_parameter <- function(old_name, new_name,
                             assumptions = "central", model = "deterministic") {
  
  info_filename <-
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "proposal.csv", sep = "/")
  
  
  parameters_info <- read_csv(info_filename)
  parameters_info$name[parameters_info$name == old_name] <- new_name
  parameters_info <- dplyr::arrange(parameters_info,
                                    region, name)

  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  parameters_proposal <- read_csv(proposal_filename)
  parameters_proposal$name[parameters_proposal$name == old_name] <- new_name
  parameters_proposal <- dplyr::arrange(parameters_proposal, region, name)
  names(parameters_proposal)[names(parameters_proposal) == old_name] <- new_name
  col_order <- c(1, 2, 2 + order(names(parameters_proposal)[-c(1,2)]))
  parameters_proposal <- parameters_proposal[, col_order]

  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}


add_beta <- function(beta_name, beta_initial, min = 0, max = 1,
                     proposal = NULL, factor = NULL,
                     assumptions = "central", model = "deterministic") {
  
  info_filename <-
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste("src/ZamCovid_parameters/pars", assumptions, model,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename)
  
  regions <- unique(parameters_info$region)
  regions <- regions[!is.na(regions)]
  
  if (!is.null(factor)) {
    f <- factor
  } else {
    f <- 1
  }
  
  new_par <- data.frame(
    region = regions,
    name = beta_name,
    initial = parameters_info$initial[parameters_info$name == beta_initial] * f,
    min = min,
    max = max,
    integer = FALSE,
    include = TRUE
  )
  
  parameters_info <- rbind(parameters_info, new_par)
  parameters_info <- dplyr::arrange(parameters_info,
                                    region, name)
  
  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  
  parameters_proposal <- read_csv(proposal_filename)
  
  if (is.null(proposal)) {
    avg_prop <- NULL
    regions <- unique(parameters_proposal$region)
    for (i in regions[!is.na(regions)]) {
      
      avg_prop[i] <- parameters_proposal[
        parameters_proposal$region == i &
          parameters_proposal$name == beta_initial, beta_initial]
    }
    proposal = abs(mean(avg_prop))
  }
  
  parameters_proposal[[beta_name]] <- 0
  k <- parameters_proposal$name == "beta1"
  new_prop <- parameters_proposal[k, ]
  new_prop$name <- beta_name
  new_prop[, -c(1, 2)] <- 0
  new_prop[[beta_name]] <- proposal
  
  parameters_proposal <- rbind(parameters_proposal, new_prop)
  parameters_proposal <- dplyr::arrange(parameters_proposal, region, name)
  col_order <- c(1, 2, 2 + order(names(parameters_proposal)[-c(1,2)]))
  parameters_proposal <- parameters_proposal[, col_order]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
}


swap_regions_proposal_and_info <- function(region_edit, region_source,
                                           assumptions = "central", 
                                           deterministic = TRUE) {
  
  if (deterministic) {
    deterministic <- "deterministic"
  } else {
    deterministic <- "stochastic"
  }
  
  info_filename <-
    paste("src/ZamCovid_parameters/pars", assumptions, deterministic,
          "info.csv", sep = "/")
  proposal_filename <- 
    paste("src/ZamCovid_parameters/pars", assumptions, deterministic,
          "proposal.csv", sep = "/")
  
  parameters_info <- read_csv(info_filename)
  
  parameters_info[parameters_info$region == region_edit, "initial"] <-
    parameters_info[parameters_info$region == region_source, "initial"]
  
  write.csv(parameters_info, info_filename, row.names = FALSE)
  
  
  parameters_proposal <- read_csv(proposal_filename)
  
  par_cols <- 3:ncol(parameters_proposal)
  
  parameters_proposal[parameters_proposal$region == region_edit, par_cols] <-
    parameters_proposal[parameters_proposal$region == region_source, par_cols]
  
  write.csv(parameters_proposal, proposal_filename, row.names = FALSE)
  
}


nudge_info <- function(region, pars, factor = NULL, value = NULL,
                       assumptions = "central", deterministic = TRUE) {
  
  if (deterministic) {
    deterministic <- "deterministic"
  } else {
    deterministic <- "stochastic"
  }
  
  info_filename <-
    paste("src/ZamCovid_parameters/pars", assumptions, deterministic,
          "info.csv", sep = "/")
  
  info <- read_csv(info_filename)
  
  for (p in pars) {
    
    init <- info[info$region == region & info$name == p, "initial"]
    new <- NULL
    if (!is.null(factor)) {
      new <- init * factor
      info[info$region == region & info$name == p, "initial"] <- new
    } else if (!is.null(value)) {
      new <- value
      info[info$region == region & info$name == p, "initial"] <- new
    } else {
      stop("Need to specify either a scalar factor or value")
    }
    print(paste(p, "for", region, "changed from", init, "to", new))
  }
  
  write.csv(info, info_filename, row.names = FALSE)
  
}


## TODO: This function hasn't been used for a while; update or remove
zero_a_proposal <- function(name, regions,
                            assumptions = "central", deterministic = TRUE) {
  
  if (deterministic) {
    deterministic <- "deterministic"
  } else {
    deterministic <- "stochastic"
  }
  proposal_filename <- 
    paste("src/ZamCovid_parameters/pars", assumptions, deterministic,
          "proposal.csv", sep = "/")
  
  
  proposal <- read_csv(proposal_filename)
  
  proposal[proposal$region %in% regions, name] <- 0
  proposal[proposal$region %in% regions & proposal$name == name,
           !names(proposal) %in% c("region", "name")] <- 0
  proposal <- dplyr::arrange(proposal, region, name)
  col_order <- c(1, 2, 2 + order(names(proposal)[-c(1,2)]))
  proposal <- proposal[, col_order]
  
  write.csv(proposal, proposal_filename, row.names = FALSE)
}
