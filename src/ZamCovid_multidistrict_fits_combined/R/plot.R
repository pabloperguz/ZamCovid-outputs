# Get model trajectories into long format tibble
traj_to_long <- function(t, tmp, dates, data = NULL) {
  
  if (all(is.null(data))) {data = rep(NA_real_, dim(tmp)[2])}
  
  data.frame(
    date = dates,
    state = t,
    mean = colMeans(tmp),
    lb = matrixStats::colQuantiles(tmp, probs = 0.025),
    ub = matrixStats::colQuantiles(tmp, probs = 0.975),
    data = data
  )
}


calc_seropos <- function(tmp, dates, name) {
  
  pos <- tmp[, grep("pos", colnames(tmp))]
  tot <- tmp[, grep("tot", colnames(tmp))]
  pos[is.na(pos)] <- 0
  tot[is.na(tot)] <- 0
  
  data.frame(date = dates) %>%
    cbind(., as.data.frame(Hmisc::binconf(pos, tot))) %>%
    `colnames<-`(c("date", "data_mean", "data_lb", "data_ub")) %>%
    mutate(state = name)
}


forest_plot_labels <- function(dat) {
  par_names <- unique(dat$parameters$proposal$name)
  beta_date <- ZamCovid:::numeric_date_as_date(dat$samples[[1]]$info$beta_date)
  
  n_betas <- length(beta_date)
  
  labels <- lapply(seq_len(n_betas), 
                   function(x) {
                     date_x <- as.character(beta_date[x])
                     bquote(beta[.(x)] ~ (.(date_x)))
                   })
  names(labels) <- paste0("beta", seq_len(n_betas))
  
  plus <- "+"
  labels <- c(labels,
              alpha_D = expression(alpha[D]),
              p_G_D = expression(p[paste(G[D], ",", 1)]^max),
              rho_pcr_tests = expression(rho[PCR[test]]),
              start_date = expression(t["Seed"])
  )
  labels
}

get_par_labels <- function(dat) {
  
  par_names <- unique(dat$fit$parameters$proposal$name)
  beta_date <- ZamCovid:::numeric_date_as_date(dat$fit$samples$info$beta_date)
  
  n_betas <- length(beta_date)
  
  labels <- lapply(seq_len(n_betas), 
                   function(x) {
                     date_x <- as.character(beta_date[x])
                     bquote(beta[.(x)] ~ (.(date_x)))
                   })
  names(labels) <- paste0("beta", seq_len(n_betas))
  
  plus <- "+"
  labels <- c(labels,
              alpha_D = expression(alpha["D"]),
              mu_D_1 = expression(mu["D,1"]),
              mu_D_2 = expression(mu["D,2"]),
              p_G_D = expression(p[paste(G[D])]^max),
              start_date = expression(t["seed"])
  )
  labels
}

plot_forest <- function(dat, plot_type = "all", nrow = 4, par_labels = list()) {
  
  stopifnot(plot_type %in% c("all", "betas", "non_betas"))
  
  info <- dat$samples[[1]]$info
  samples <- dat$samples
  date <- info$date
  par_names <- colnames(samples[[1]]$pars)
  
  beta_date <- ZamCovid:::numeric_date_as_date(info$beta_date)
  beta_names <- par_names[substr(par_names, 1, 4) == "beta"]
  beta_names <- beta_names[order(as.numeric(gsub("beta", "", beta_names)))]
  
  hps <- dat$parameters$prior
  pars_info <- dat$parameters$info
  
  
  ## Can set the xmax for specific parameters here. For any left as NA, it will
  ## instead just use the maximum from the parameters info
  par_max <- rep(NA, length(par_names))
  names(par_max) <- par_names
  par_max[beta_names] <- 0.25
  par_max["rho_pcr_tests"] <- 0.02
  
  n_districts <- length(samples)
  
  extract_sample <- function(par_name) {
    lapply(samples, function(x) as.numeric(x$pars[, par_name]))
  }
  
  ylim <- c(0.5, n_districts + 0.5)
  labels <- c(kabwe = "Kab", lusaka = "Lus", livingstone = "Liv", 
              ndola = "Ndo", solwezi = "Sol")
  
  plot_axis <- function() {
    yax_pos <- 0.9
    plot(0, 0, type = "n",
         ylab = "",
         xlab = "",
         xlim = c(0, 0),
         ylim = ylim,
         axes = FALSE
    )
    axis(side = 2, at =  at, labels = labels, las = 1, pos = yax_pos)
  }
  
  col_line <- "red3"
  
  plot_par <- function(par_name) {
    if (par_name == "start_date") {
      plot_start_date()
    } else {
      # browser()
      par <- extract_sample(par_name)
      
      par_info <- subset(pars_info, pars_info$name == par_name)
      if (is.na(par_max[[par_name]])) {
        xmax <- max(par_info$max)
      } else {
        xmax <- par_max[[par_name]]
      }
      
      xmin <- min(par_info$min)
      
      if (par_name %in% names(par_labels)) {
        xlab <- par_labels[[par_name]]
      } else {
        if (grepl("^beta", par_name)) {
          k <- as.numeric(gsub("beta", "", par_name))
          xlab <- paste0(par_name, " (", beta_date[k], ")")
        } else {
          xlab <- par_name
        }
      }
      
      plot(0, 0, type = "n",
           ylab = "",
           xlab = xlab,
           xlim = c(xmin, xmax),
           ylim = ylim,
           yaxt = "n"
      )
      
      jitter <- 0.5
      distrits <- names(par)
      hp <- subset(hps, hps$name == par_name)
      if (is.na(hp$district[1])) {
        hp <- hp[rep(1, length(districts)), ]
        hp$district <- districts
        col <- "red"
      } else {
        col <- "grey20"
      }
      rownames(hp) <- hp$district
      hp <- hp[districts, ] # sort in correct order
      if (hp$type[1] == "beta") {
        shape1 <- hp$beta_shape1
        shape2 <- hp$beta_shape2
        if (!(all(shape1 == 1) && all(shape2 == 1))) {
          prior <- mapply(qbeta,
                          shape1 = shape1,
                          shape2 = shape2,
                          MoreArgs = list(p = c(0.025, 0.975)),
                          SIMPLIFY = TRUE)
          
          
          segments(x0 = prior[1, ],
                   y0 = at - jitter, y1 = at + jitter,
                   col = col_line, lty = 2, lwd = 1, lend = 2)
          segments(x0 = prior[2, ],
                   y0 = at - jitter, y1 = at + jitter,
                   col = col_line, lty = 2, lwd = 1, lend = 2)
        }
      }
      if (hp$type[1] == "gamma") {
        shape <- hp$gamma_shape
        scale <- hp$gamma_scale
        prior <- mapply(qgamma,
                        shape = shape,
                        scale = scale,
                        MoreArgs = list(p = c(0.025, 0.975)),
                        SIMPLIFY = TRUE)
        
        
        segments(x0 = prior[1, ],
                 y0 = at - jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
        segments(x0 = prior[2, ],
                 y0 = at - jitter, y1 = at + jitter,
                 col = col_line, lty = 2, lwd = 1, lend = 2)
      }
      
      mapply(FUN = plot_ci_bar, res = par, at = at, width = 0.1,
             col = col)
    }
  }
  
  plot_start_date <- function() {
    numeric_start_date <- extract_sample(par_name = "start_date")
    if ("start_date" %in% names(par_labels)) {
      xlab <- par_labels$start_date
    } else {
      xlab <- "start_date"
    }
    plot(x = sircovid::sircovid_date("2020-01-01"),
         y = 1,
         type = "n",
         ylab = "",
         xlab = xlab,
         xlim = sircovid::sircovid_date(c("2020-01-15", "2020-03-01")),
         ylim = ylim,
         yaxt = "n")
    mapply(FUN = plot_ci_bar, res = numeric_start_date, at = at, width = 0.1)
  }
  
  hps <- dat$parameters$prior
  pars_info <- dat$parameters$info
  
  op <- par(bty = "n", mar = c(3, 0, 1, 0), mgp = c(2, 0.75, 0),
            oma = c(2, 0, 3, 3))
  on.exit(par(op))
  
  if (plot_type == "all") {
    pars_to_plot <- par_names
  } else if (plot_type == "betas") {
    pars_to_plot <- beta_names
  } else if (plot_type == "non_betas") {
    pars_to_plot <- setdiff(par_names, beta_names)
  }
  
  if ("start_date" %in% pars_to_plot) {
    pars_to_plot <- c("start_date", pars_to_plot[pars_to_plot != "start_date"])
  }
  
  npar <- length(pars_to_plot)
  nrow <- nrow
  plot_per_row <- ceiling(npar / nrow)
  colwidth <- 64 / plot_per_row
  
  reps <- rep(c(3.5, rep(colwidth, plot_per_row)), nrow)
  
  layout(mat = matrix(rep(seq_along(reps), reps),
                      nrow = nrow, byrow = TRUE),
         heights = rep(3, nrow),
         widths = c(4, rep(1, plot_per_row * colwidth))
  )
  at <- seq_len(n_districts)
  
  for (i in seq_len(nrow)) {
    plot_axis()
    par(mar = c(3, 0, 1, 0.5))
    pars_row <- seq((i - 1) * plot_per_row + 1, min(npar, i * plot_per_row))
    mapply(plot_par, par_name = pars_to_plot[pars_row])
  }
  
  mtext(text = paste("Inferred epidemic parameters for Zambia districts at", date),
        side = 3, line = 0, outer = TRUE, cex = 1.1)
}


plot_ci_bar <- function(res, at, width = 1,
                        min = 0.025, max = 0.975, col = "grey20",
                        segments = FALSE, pt_col = NULL, horiz = TRUE, ...) {
  cols  <- c("grey80", col)
  qs <- quantile(res,
                 probs = seq(min, max, by = 0.005),
                 na.rm = TRUE)
  
  palette <- grDevices::colorRampPalette(cols)
  if (segments) {
    segments(y0 = at, x0 = min(res), x1 = max(res), col = cols[2])
    points(y = rep(at, 2), x = range(res), col = cols[2], pch = "-")
  }
  ci_bands(quantiles = cbind(qs, qs),
           y = at + c(-1, 1) * width,
           palette = palette, leg = FALSE, horiz = horiz)
  
  if (is.null(pt_col)) pt_col <- col
  if (horiz) {
    points(y = at, x = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  } else {
    points(x = at, y = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  }
}



plot_serology <- function(dat, data, all = FALSE, ylim = 1) {
  
  samples <- dat$fit$samples
  
  which <- c("all", "15_19", "20_29",
             "30_39", "40_49", "50_plus")
  sero_cols <- grep("sero", colnames(data), value = TRUE)
  dates <- ZamCovid:::numeric_date_as_date(samples$trajectories$date)[-1]
  data <- data %>%
    select(date = date_string, sero_cols) %>%
    filter(as.Date(date) %in% dates)
  
  values <- c("cadetblue", "darkgoldenrod3", "chartreuse4")
  breaks <- c("Model", "Data (fitted)", "Data (not fitted)")
  labels <- c("All", "15 to 19",
              "20 to 29", "30 to 39", "40 to 49",
              "50 and older")
  
  if (all) {
    which <- which[1]
    labels <- labels[1]
  }
  
  data_plot <- NULL
  for (w in which) {
    keep <- grep(w, colnames(data), value = TRUE)
    name <- paste0("sero_pos_", w)
    tmp <- calc_seropos(data[, keep], dates, name)
    rownames(tmp) <- NULL
    data_plot <- rbind(data_plot, tmp)
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
    n <- p[[paste0("N_tot_", gsub("sero_pos_", "", t))]]
    sero_pos <- states[t, , -1]
    tmp <- ((sens *  sero_pos + (1 - spec) * (n - sero_pos)) / n) - 0.01
    ret <- traj_to_long(t, tmp, dates) %>% select(!data)
    df <- rbind(df, ret)
  }
 
  df <- left_join(df, data_plot, by = c("date", "state")) %>%
    mutate(col = case_when(state == "sero_pos_over15" ~ "Data (not fitted)",
                           TRUE ~ "Data (fitted)"))
  df$state <- factor(df$state, levels = unique(df$state), labels = labels)
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data_mean, col = col), size = 0.9, alpha = 0.5) +
    geom_errorbar(aes(ymin = data_lb, ymax = data_ub, 
                        col = col), linewidth = 0.3, alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, ylim),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Seropositivity") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
}


plot_pcr_positivity <- function(dat, data, all = FALSE, ylim = 1) {
  
  samples <- dat$fit$samples
  
  dates <- ZamCovid:::numeric_date_as_date(samples$trajectories$date)[-1]
  
  data <- data %>%
    select(date = date_string, pcr_pos_all, pcr_tot_all) %>%
    filter(as.Date(date) %in% dates)
  
  values <- c("cadetblue", "darkgoldenrod3", "chartreuse4")
  breaks <- c("Model", "Data (fitted)", "Data (not fitted)")
  

  data_plot <- calc_seropos(data, dates, "pcr_pos_all")
  rownames(data_plot) <- NULL


  states <- samples$trajectories$state
  nms <- rownames(states)
  
  p <- samples$predict$transform(samples$pars[1, ])
  sens <- p$pcr_sensitivity
  spec <- p$pcr_specificity
  
  pcr_pos <- states["pcr_pos_all", , -1]
  n <- p[["N_tot_all"]]
  tmp <- ((sens *  pcr_pos + (1 - spec) * (n - pcr_pos)) / n) - 0.01
  
  df <- traj_to_long("pcr_pos_all", tmp, dates) %>%
    select(!data) %>%
    left_join(., data_plot, by = c("date", "state")) %>%
    mutate(col = "Data (fitted)")
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data_mean, col = col), size = 0.9, alpha = 0.5) +
    geom_errorbar(aes(ymin = data_lb, ymax = data_ub, 
                      col = col), linewidth = 0.3, alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, ylim),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_color_manual(values = values, breaks = breaks) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "PCR positivity") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
}


plot_deaths <- function(dat, data, week_only = FALSE, legend = NULL) {
  
  if(is.null(legend)) {legend <- "none"}
  
  which <- paste0("deaths_", c("all", "hosp"))
  values <- c("cadetblue", "darkgoldenrod3", "chartreuse4")
  breaks <- c("Model", "Data (fitted)", "Data (not fitted)")
  
  data_plot <- data %>%
    select(date = date_string,
           deaths_all_inc = deaths_all,
           deaths_hosp_inc = deaths_hosp) %>%
    mutate(deaths_comm_inc = deaths_all_inc - deaths_hosp_inc,
           date = as.Date(date))
  
           
  states <- dat$fit$samples$trajectories$state
  nms <- rownames(states)
  
  df <- traj_to_long("deaths_all_inc",
                     states["deaths_all_inc", , -1],
                     data_plot$date,
                     data_plot$deaths_all_inc)
  
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
    labs(x = "", y = "Daily deaths") +
    theme_minimal() +
    theme(legend.position = legend,
          legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_blank())
  
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
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = "Model"), alpha = 0.4, show.legend = FALSE) +
    geom_point(aes(y = data, col = "Data (not fitted)"), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_color_manual(values = values, breaks = breaks,
                       guide =
                         guide_legend(override.aes = list(linetype = c(1, NA),
                                                          size = rep(2, 2),
                                                          shape = c(NA, 16)))) +
    scale_fill_manual(values = values, breaks = breaks) +
    labs(x = "", y = "Weekly deaths (not fitted)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
  
  
  if (week_only) {
    weekly + labs(y = "Weekly deaths (all-cause)", title = "") +
      theme(legend.position = legend,
            legend.title = element_blank())
  } else {
    daily / weekly
  }
}


plot_deaths_disag <- function(dat, plot_age = TRUE) {
  
  which <- c("deaths_all_inc", "base_death_inc", "deaths_covid_inc",
             paste0("D_", seq(0, 75, 5)))
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$samples$trajectories$date)[-1]
  data <- rep(NA_real_, length(dates))
  
  states <- dat$fit$samples$trajectories$state
  
  df <- NULL
  for (t in which) {
    
    if (t == "deaths_covid_inc") {
      tmp <- states["deaths_comm_inc", , -1] + states["deaths_hosp_inc", , -1]
    } else {
      tmp <- states[t, , -1]
    }
    ret <- traj_to_long(t, tmp, dates, data) %>% select(!data)
    df <- rbind(df, ret)
  }
  
  
  deaths_agg <- df %>%
    filter(state %in% c("base_death_inc", "deaths_covid_inc")) %>%
    select(date, state, mean) %>%
    group_by(date, state) %>%
    summarise(n = sum(mean)) %>%
    mutate(percent = n / sum(n)) %>%
    filter(state == "deaths_covid_inc")
  
  agg <- ggplot(deaths_agg, aes(x = date)) +
    geom_area(aes(y = percent, fill = state), alpha = 0.75) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    scale_fill_manual(values = "brown4") +
    labs(x = "", y = "Covid-19 deaths (% all deaths)") +
    theme_minimal() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))
  
  deaths_age <- df %>%
    filter(!state %in% c("deaths_all_inc", "base_death_inc", "deaths_covid_inc")) %>%
    select(date, state, mean) %>%
    group_by(date, state) %>%
    summarise(n = sum(mean)) %>%
    mutate(percent = n / sum(n))
  
  deaths_age$state <- factor(deaths_age$state,
                             levels = paste0("D_", seq(0, 75, 5)),
                             labels = c(paste0(seq(0, 70, 5),
                                               " to ", seq(4, 74, 5)),
                                        "75+"))
  
  
  age <- ggplot(deaths_age, aes(x = date)) +
    geom_area(aes(y = percent, fill = state), alpha = 0.75) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format()) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_fill_viridis_d(direction = 1) +
    labs(x = "", y = "Deaths by age (% daily Covid-19 deaths)") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
  
  if (plot_age) {
    ((agg + theme(axis.text.x = element_blank())) / age) +
      plot_layout(heights = c(0.5, 1))
  } else {agg}
  
  
}


plot_deaths_disag_inc <- function(dat, data_full) {
  
  dates <- ZamCovid:::numeric_date_as_date(dat$fit$samples$trajectories$date)[-1]
  which <- c("deaths_0_14", "deaths_15_39", "deaths_40_59", "deaths_60_plus")
  
  data <- data_full[, c("date", which)] %>%
    mutate(date = as.Date(date)) %>%
    filter(date %in% dates)
  
  states <- dat$fit$samples$trajectories$state
  df <- NULL
  for (t in which) {
    
    if (t == which[1]) {
      tmp <- states[paste0("D_", seq(0, 10, 5)), , -1]
    } else if (t == which[2]) {
      tmp <- states[paste0("D_", seq(15, 35, 5)), , -1]
    } else if (t == which[3]) {
      tmp <- states[paste0("D_", seq(40, 55, 5)), , -1]
    } else if (t == which[4]) {
      tmp <- states[paste0("D_", seq(60, 75, 5)), , -1]
    }
    
    tmp <- apply(tmp, MARGIN = c(2, 3), sum)
    # tmp <- rbind(rep(0, dim(tmp)[2]), apply(tmp, MARGIN = 2, diff))
    # browser()
    ret <- traj_to_long(t, tmp, dates, data[, t]) %>%
      filter(date >= as.Date("2020-03-15"))
    ret$data <- cumsum(ret$data)
    
    df <- rbind(df, ret)
  }
  
  base_deaths <- data.frame(
    date = dates,
    deaths = colMeans(states["base_death_inc", , -1])) %>%
    filter(date >= as.Date("2020-03-15"))
  
  base_deaths <- df %>% 
    select(date, state, data) %>%
    group_by(date, state) %>%
    summarise(n = sum(data)) %>%
    mutate(percent = n / sum(n)) %>%
    ungroup() %>%
    select(date, name = state, value = percent) %>%
    pivot_wider() %>%
    left_join(base_deaths, .) %>%
    mutate_each(funs(cumsum(. * deaths)), starts_with("deaths_")) %>%
    # mutate_each(funs(. * deaths), starts_with("deaths_")) %>%
    select(!deaths) %>%
    pivot_longer(!date, names_to = "state", values_to = "baseline")
  
  df <- left_join(df, base_deaths)
  # browser()
  # x <- left_join(df, base_deaths) %>%
  #   mutate(week = lubridate::floor_date(date, "week", week_start = 3)) %>%
  #   select(!date) %>%
  #   pivot_longer(!c(week, state)) %>%
  #   group_by(week, state, name) %>%
  #   summarise(value = sum(value)) %>%
  #   ungroup() %>%
  #   pivot_wider()
  
  ggplot(df, aes(x = date)) +
    geom_line(aes(y = mean + baseline, col = state)) +
    geom_ribbon(aes(ymin = lb + baseline, ymax = ub + baseline, fill = state), alpha = 0.3) +
    geom_point(aes(y = data, col = state), alpha = 0.3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y",
                 limits = c(as.Date("2020-01-01"), NA)) +
    scale_fill_viridis_d(direction = 1) +
    scale_color_viridis_d(direction = 1) +
    facet_wrap(~state, scales = "free_y") +
    labs(x = "", y = "Cumulative deaths by age (all-cause)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10))
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
    theme(legend.position = "right",
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
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.015),
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
          axis.text.x = element_blank())
  
  p_reinf <- ggplot(reinf, aes(x = date)) +
    geom_ribbon(aes(ymax = 1, ymin = mean / inf$mean, fill = "1st infections"), alpha = 0.7) +
    geom_ribbon(aes(ymin = 0, ymax = mean / inf$mean, fill = "Re-infections"), alpha = 0.7) +
    scale_fill_manual(values = c("royalblue4", "orange2")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "% daily infections", x = "") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.line = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12))

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




plot_ci_bar <- function(res, at, width = 1,
                        min = 0.025, max = 0.975, col = "grey20",
                        segments = FALSE, pt_col = NULL, horiz = TRUE, ...) {
  cols  <- c("grey80", col)
  qs <- quantile(res,
                 probs = seq(min, max, by = 0.005),
                 na.rm = TRUE)
  
  palette <- grDevices::colorRampPalette(cols)
  if (segments) {
    segments(y0 = at, x0 = min(res), x1 = max(res), col = cols[2])
    points(y = rep(at, 2), x = range(res), col = cols[2], pch = "-")
  }
  ci_bands(quantiles = cbind(qs, qs),
           y = at + c(-1, 1) * width,
           palette = palette, leg = FALSE, horiz = horiz)
  
  if (is.null(pt_col)) pt_col <- col
  if (horiz) {
    points(y = at, x = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  } else {
    points(x = at, y = mean(res), col = "white", pch = 23, bg = pt_col, ...)
  }
}


ci_bands <- function (quantiles, y, palette = NULL, cols = NULL, leg = TRUE, 
                      leg_y = 0, leg_x = 1, horiz = TRUE, ...) 
{
  yy <- c(y, rev(y))
  yy <- c(yy, yy[1])
  n_bands <- (nrow(quantiles) - 1)/2 + 1
  if (!is.null(palette)) {
    cols <- do.call(what = palette, args = list(n = n_bands))
  }
  for (band in seq_len(n_bands)) {
    x1 <- quantiles[band, ]
    x2 <- quantiles[nrow(quantiles) + 1 - band, ]
    x2 <- rev(x2)
    x2 <- c(x2, x1[1])
    if (horiz) {
      polygon(y = yy, x = c(x1, x2), col = cols[band], 
              border = NA)
    }
    else {
      polygon(x = yy, y = c(x1, x2), col = cols[band], 
              border = NA)
    }
  }
  if (leg) {
    leg_cols <- which(row.names(quantiles) %in% leg)
    leg <- c(row.names(quantiles)[1], seq(5, 50, 5), "%")
    leg[seq(2, 10, 2)] <- ""
    legend(y = leg_y, x = leg_x, pch = 15, col = cols[leg_cols], 
           legend = leg, border = NA, bty = "n", ...)
  }
}

