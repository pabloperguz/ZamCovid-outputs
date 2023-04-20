load_mcmc_parameters <- function(assumptions, deterministic) {
  
  message(sprintf("Assumptions are with '%s' parameters", assumptions))
  
  deterministic <- ifelse(deterministic, "deterministic", "stochastic")
  path_pars <- file.path("pars", assumptions, deterministic)
  
  info <- read_csv(file.path(path_pars, "info.csv"))
  proposal <- read_csv(file.path(path_pars, "proposal.csv"))
  
  prior <- create_priors(info)
  proposal <- update_proposal(info, proposal)
  
  ret <- list(info = info,
              proposal = proposal,
              prior = prior)
  
  ret
}


update_proposal <- function(info, proposal) {
  msg <- setdiff(info$name, proposal$name)
  if (length(msg) > 0) {
    message(sprintf("Adding %d parameters to proposal matrix", length(msg)))
    extra <- expand.grid(region = unique(proposal$region), name = msg)
    extra[names(proposal)[-(1:2)]] <- 0
    proposal <- rbind(proposal, extra)
    proposal[msg] <- 0
    proposal <- proposal[order(proposal$region, proposal$name), ]
    col_order <- c(1, 2, 2 + order(names(proposal)[-c(1, 2)]))
    proposal <- proposal[, col_order, drop = FALSE]
  }
  proposal
}
