plot_fit_traces <- function(samples) {
  
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    n_chains <- length(unique(samples$chain))
  }
  
  cols <- rev(viridisLite::viridis(n_chains))
  
  reorder_beta <- function(nms) {
    i <- grep("^beta[0-9]+$", nms)
    j <- order(as.integer(sub("^beta", "", nms[i])))
    k <- seq_along(nms)
    k[i] <- i[j]
    k
  }
  
  i <- reorder_beta(colnames(samples$pars))
  pars <- samples$pars[, i]
  nms <- colnames(pars)
  probs <- samples$probabilities
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  new_grid <- function(n, title) {
    par(mfrow = rep(ceiling(sqrt(n + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1 + as.integer(title), 1))
  }
  
  plot_traces1 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    ess <- coda::effectiveSize(coda::as.mcmc(traces))
    
    if (name == "log_likelihood") {
      main <- ""
    } else {
      main <- sprintf("ess = %s", round(sum(ess)))
    }
    matplot(traces, type = "l", lty = 1,
            xlab = "Iteration", bty = "n",
            ylab = name, col = cols,
            main = main,
            font.main = 1)
    rug(samples$iteration[samples$chain == 1], ticksize = 0.1)
  }
  
  new_grid(length(nms), FALSE)
  for (nm in nms) {
    plot_traces1(samples$pars[, nm], nm)
  }
  plot_traces1(probs[, "log_likelihood"], "log_likelihood")
}


plot_serology <- function(samples, data_fit,
                          which = c("over15", "15_19", "20_29",
                                    "30_39", "40_49", "50_plus")) {
  
  sero_cols <- grep("sero", colnames(data_fit), value = TRUE)
  data_fit <- data_fit[, c("date", "date_string", sero_cols)]
  
  data_plot <- data_fit[, c("date", "date_string")]
  for (w in which) {
    keep <- grep(w, colnames(data_fit), value = TRUE)
    tmp <- data_fit[, keep]
    
    pos <- tmp[, grep("pos", colnames(tmp))]
    tot <- tmp[, grep("tot", colnames(tmp))]
    out <- data.frame(pos / tot)
    colnames(out) <- paste0("sero_pos_", w)
    
    data_plot <- cbind(data_plot, out)
  }
  
  states <- samples$trajectories$state
  nms <- rownames(states)
  
  p <- samples$predict$transform(samples$pars[1, ])
  sens <- p$sero_sensitivity
  spec <- p$sero_specificity
  
  
  # Get model trajectories and calculate seropositivity
  traj_to_long <- function(t, tmp) {
    data.frame(
      date = as.Date(data_plot$date_string),
      state = t,
      mean = colMeans(tmp),
      lb = matrixStats::colQuantiles(tmp, probs = 0.025),
      ub = matrixStats::colQuantiles(tmp, probs = 0.975),
      data = data_plot[, t]
    )
  }
  
  traj <- paste0("sero_pos_", which)
  sero_pos <- samples$trajectories$state[traj, , ]
  
  df <- NULL
  for (t in traj) {
    n <- p[[paste0("N_tot_", gsub("sero_pos_", "", t))]]
    sero_pos <- states[t, , -1]
    tmp <- ((sens *  sero_pos + (1 - spec) * (n - sero_pos)) / n) - 0.01
    ret <- traj_to_long(t, tmp)
    df <- rbind(df, ret)
  }
  
  # res_infs <- sample$trajectories$state["infections", , ] / p$N_tot_all * 100
  df$state <- factor(df$state, levels = unique(df$state))
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean, col = state)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = state), alpha = 0.4) +
    geom_point(aes(y = data, col = state), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 0.3),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    facet_wrap(~state, scales = "free_y") +
    labs(x = "", y = "Proportion sero-positive") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
}


get_new_pars <- function(samples, priors) {
  
  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]
  
  vcv <- cov(samples$pars)
  dimnames(vcv) <- NULL
  
  for (i in 1:length(priors)) {
    nm <- priors[[i]]$name
    priors[[i]]$initial <- initial[[nm]]
  }
  
  list(vcv = vcv,
       priors = priors)
}
