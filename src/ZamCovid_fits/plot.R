# Get model trajectories into long format tibble
traj_to_long <- function(t, tmp, dates, data) {
  data.frame(
    date = dates,
    state = t,
    mean = colMeans(tmp),
    lb = matrixStats::colQuantiles(tmp, probs = 0.025),
    ub = matrixStats::colQuantiles(tmp, probs = 0.975),
    data = data
  )
}


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
  
  i <- reorder_beta(colnames(samples$pars_full))
  pars <- samples$pars_full[, i]
  nms <- colnames(pars)
  probs <- samples$probabilities_full
  
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
  }
  
  new_grid(length(nms), FALSE)
  for (nm in nms) {
    plot_traces1(samples$pars_full[, nm], nm)
  }
  plot_traces1(probs[, "log_likelihood"], "log_likelihood")
}


plot_serology <- function(samples, data_fit,
                          which = c("over15", "15_19", "20_29",
                                    "30_39", "40_49", "50_plus")) {
  
  sero_cols <- grep("sero", colnames(data_fit), value = TRUE)
  data_fit <- data_fit[, c("date", "date_string", sero_cols)]
  
  data_plot <- NULL
  for (w in which) {
    keep <- grep(w, colnames(data_fit), value = TRUE)
    tmp <- data_fit[, keep]
    
    pos <- tmp[, grep("pos", colnames(tmp))]
    tot <- tmp[, grep("tot", colnames(tmp))]
    pos[is.na(pos)] <- 0
    tot[is.na(tot)] <- 0
    
    out <- data.frame(as.Date(data_fit$date_string)) %>%
      cbind(., as.data.frame(Hmisc::binconf(pos, tot))) %>%
      `colnames<-`(c("date", "data_mean", "data_lb", "data_ub")) %>%
      mutate(state = paste0("sero_pos_", w))
    
    data_plot <- rbind(data_plot, out)
  }
  
  states <- samples$trajectories$state
  nms <- rownames(states)
  
  p <- samples$predict$transform(samples$pars[1, ])
  sens <- p$sero_sensitivity
  spec <- p$sero_specificity
  
  traj <- paste0("sero_pos_", which)
  sero_pos <- samples$trajectories$state[traj, , ]
  
  dates_vect <- unique(data_plot$date)
  df <- NULL
  for (t in traj) {
    dates <- unique(as.Date(data_plot$date))
    data <- data_plot %>%
      filter(state == t)
    n <- p[[paste0("N_tot_", gsub("sero_pos_", "", t))]]
    sero_pos <- states[t, , -1]
    tmp <- ((sens *  sero_pos + (1 - spec) * (n - sero_pos)) / n) - 0.01
    ret <- traj_to_long(t, tmp, dates, data)
    df <- rbind(df, ret)
  }
  
  df <- left_join(df, data_plot)
  df$state <- factor(df$state, levels = unique(df$state))
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean, col = state)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = state), alpha = 0.4) +
    geom_point(aes(y = data_mean, col = state), size = 0.9, alpha = 0.5) +
    geom_errorbar(aes(ymin = data_lb, ymax = data_ub, 
                        col = state), linewidth = 0.3, alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    facet_wrap(~state, scales = "free_y") +
    labs(x = "", y = "Proportion sero-positive") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
}


plot_deaths <- function(dat, data_fit) {
  
  which <- paste0("deaths_", c("all", "hosp"))
  values <- c("cadetblue", "darkgoldenrod3", "chartreuse4")
  breaks <- c("Model", "Data (fitted)", "Data (not fitted)")
  
  data_plot <- data_fit[, c("date", "date_string", which)]
  data_plot$deaths_comm_inc <- data_plot$deaths_all - data_plot$deaths_hosp
  colnames(data_plot) <- c("date", "date_string", "deaths_all_inc",
                           "deaths_hosp_inc", "deaths_comm_inc")
  
  states <- dat$fit$samples$trajectories$state
  nms <- rownames(states)
  
  traj <- c("deaths_all_inc")
  dates <- as.Date(data_plot$date_string)
  df <- NULL
  for (t in traj) {
    data <- data_plot[, t]
    ret <- traj_to_long(t, states[t, , -1], dates, data)
    df <- rbind(df, ret)
  }
  
  df$state <- factor(df$state, levels = unique(df$state))
  df <- df %>%
    filter(!is.na(data))
  
  ylim <- max(df$ub) * 1.25
  
  daily <- ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data, col = "Data (fitted)"), size = 0.7, alpha = 0.9) +
    geom_point(aes(y = NA_real_, col = "Data (not fitted)"), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_color_manual(values = values, breaks = breaks,
                       guide =
                         guide_legend(override.aes = list(linetype = c(1, NA, NA),
                                                          size = rep(2, 3),
                                                          shape = c(NA, 16, 16)))) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Daily (fitted)") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  ylim <- max(cumsum(df$ub), cumsum(df$data)) + 100
  
  cumulative <- ggplot(df, aes(x = date)) +
    geom_line(aes(y = cumsum(mean), col = "Model")) +
    geom_ribbon(aes(ymin = cumsum(lb), ymax = cumsum(ub), fill = "Model"), alpha = 0.4) +
    geom_point(aes(y = cumsum(data), col = "Data (not fitted)"), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Cumulative") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  weekly <- df %>%
    mutate(date = lubridate::floor_date(date, "week", week_start = 3)) %>%
    select(!state) %>%
    pivot_longer(!c(date)) %>%
    group_by(date, name) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider()
    
  weekly <- ggplot(weekly, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4) +
    geom_point(aes(y = data, col = "Data (not fitted)"), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Weekly") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  daily / weekly / cumulative +
    plot_annotation(title = "All-cause deaths (all ages)")
  
}


plot_rt <- function(dat) {

  rt <- dat$fit$rt$Rt_general
  eff_rt <- dat$fit$rt$eff_Rt_general
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$rt$time[, 1] / 4)
  
  df_to_long <- function(df, which) {
    data.frame(
      date = dates,
      mean = rowMeans(df),
      lb = matrixStats::rowQuantiles(df, probs = 0.025),
      ub = matrixStats::rowQuantiles(df, probs = 0.975),
      type = which
    )
  }
  
  df <- rbind(df_to_long(rt, "Rt_general"),
              df_to_long(eff_rt, "eff_Rt_general"))
  
  ggplot(df, aes(date)) +
    geom_line(aes(y = mean, col = type)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = type), alpha = 0.3) +
    geom_hline(yintercept = 1, linewidth = 0.5, col = "red", linetype = 3) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 20, 1)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(x = "", y = "Reproduction number") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
}


plot_severity <- function(dat, age = TRUE, xmin = "2020-04-01") {
  
  ifr <- get_severity(dat$fit$severity, "ifr")
  ifr[ifr$date < as.Date(xmin), c("mean", "lb", "ub")] <- NA_real_
  
  p1 <- ggplot(ifr, aes(date, mean)) +
    geom_line(col = "blue4") +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.3, fill = "blue4") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "Effective IFR", x = "") +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  ifr_age <- get_severity(dat$fit$severity, "ifr", TRUE) %>%
    filter(date >= as.Date(xmin)) %>%
    pivot_longer(!c(date, age), names_to = "estimate") %>%
    group_by(age, estimate) %>%
    summarise(value = mean(value)) %>%
    pivot_wider(names_from = estimate) %>%
    mutate(age = factor(age, levels = seq(0, 75, 5)))
  
  p2 <- ggplot(ifr_age, aes(y = age, col = age)) +
    geom_pointrange(aes(x = mean, xmin = lb, xmax = ub)) +
    scale_x_continuous(labels = scales::percent_format(),
                       trans = "log",
                       breaks = c(0, 0.001, 0.01, 0.05, 0.1, 0.2)) +
    labs(x = "Effective IFR", y = "Age group") +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          legend.position = "none")
  
  if (age) {
    p1 + p2
  } else {
    p1
  }
}


plot_infection_incidence <- function(dat) {
  
  
  sample <- dat$fit$samples$trajectories$state
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$samples$trajectories$date)
  data <- rep(NA_real_, length(dates))
  
  inf_traj <- c( "infections_inc", "reinfections_inc")
  
  inf <- traj_to_long("infections_inc", 
                      sample["infections_inc", ,], 
                      dates, data)
  
  reinf <- traj_to_long("reinfections_inc", 
                        sample["reinfections_inc", ,], 
                        dates, data)
  
  pop <- sum(dat$fit$parameters$base$population$n)
  
  p_inf <- ggplot(inf, aes(x = date)) +
    geom_line(aes(y = (mean / pop) * 1e3)) +
    geom_ribbon(aes(ymin = (lb / pop) * 1e3,
                    ymax = (ub / pop) * 1e3), alpha = 0.3) +
    scale_color_manual(values = "blue4") +
    scale_fill_manual(values = "blue4") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "Infection incidence (x 1,000 population)", x = "") +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  p_reinf <- ggplot(reinf, aes(x = date)) +
    geom_ribbon(aes(ymax = 1, ymin = mean / inf$mean, fill = "1st infections")) +
    geom_ribbon(aes(ymin = 0, ymax = mean / inf$mean, fill = "Re-infections")) +
    scale_fill_manual(values = c("royalblue4", "orange2")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "% daily infections", x = "") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "top",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))

  p_inf / p_reinf  
}


get_severity <- function(severity, what, by_age = FALSE) {
  
  age_bands <- ZamCovid:::model_age_bins()$start
  dates <- ZamCovid:::numeric_date_as_date(severity$date[-1])
  
  if (by_age) {
    ret <- NULL
    get <- paste0(what, "_age_", age_bands)
  } else {
    ret <- data.frame(date = dates)
    get <- what
  }
  
  for (g in get) {
    sev <- severity[[g]][-1, ]
    
    tmp <- data.frame(t(apply(sev, 1, quantile, c(0.5, 0.025, 0.975),
                              na.rm = TRUE))) %>%
      setNames(., c("mean", "lb", "ub"))
    
    if (by_age) {
     
      tmp <- tmp %>%
       mutate(date = dates,
              age = gsub("ifr_age_", "", g)) %>%
        select(date, age, mean, lb, ub)
      
     ret <- rbind(ret, tmp)
     
    } else {
      ret <- cbind(ret, tmp)
    }
  }
  
  ret
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
