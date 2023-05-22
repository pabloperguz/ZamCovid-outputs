ggtheme <- theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.7))
theme_set(ggtheme)


plot_trajectory_scenarios <- function(scenarios, which) {
  
  # Unlist objects of interest into tibbles for plotting
  traj <- unlist_scenarios(scenarios, "trajectories")
  rt <- unlist_scenarios(scenarios, "rt")
  ifr <- unlist_scenarios(scenarios, "severity")
  
  
  # Some graphic parameters
  scenario_cols <- c("grey60",
                     "goldenrod", "goldenrod4",
                     "darkolivegreen3", "darkolivegreen4",
                     "royalblue1", "royalblue4",
                     "slateblue1", "slateblue4",
                     "darkorchid1", "darkorchid4")
  
  scenario_labels <- c(
    'Central parameters',
    'Fit to deaths only', 'Fit to serology only',
    'Higher IFR', "Lower IFR",
    'Fast imm. waning', 'Slow imm. waning',
    'Higher sero. sensitivity', 'Lower sero. sensitivity',
    'Fast seroreversion', 'Slow seroreversion')
  
  state_nms <- c("Daily deaths (all cause)", "Seropositive (15+)")
  
  names(scenario_cols) <- unique(traj$scenario)
  names(scenario_labels) <- unique(traj$scenario)
  names(state_nms) <- unique(traj$state)
  cols <- scenario_cols[c("central", which)]
  
  
  # Plot trajectories
  traj_plot <- NULL
  for (s in names(state_nms)) {
    ymax <- ifelse(s == names(state_nms)[1], 20, 0.5)
    p <- traj %>%
      filter(state == s, scenario %in% c("central", which)) %>%
      ggplot(., aes(date)) +
      geom_line(aes(y = mean, col = scenario)) +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = scenario), alpha = 0.2) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
      scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
      labs(y = state_nms[s], x = "") +
      scale_color_manual(values = cols, breaks = c("central", which)) +
      scale_fill_manual(values = cols, breaks = c("central", which)) +
      facet_grid(cols = vars(scenario),
                 labeller = labeller(scenario = scenario_labels),
                 scales = "free_y") +
      theme(axis.text.x = element_blank())
    
    if (s == names(state_nms)[2]) {
      p <- p +
        theme(strip.text = element_blank()) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    }
    
    traj_plot[[s]] <- p
  }

  
  ifr_p <- ifr %>%
    filter(scenario %in% c("central", which)) %>%
    ggplot(., aes(date)) +
    geom_line(aes(y = mean, col = scenario)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = scenario), alpha = 0.2) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02),
                       labels = scales::percent_format(accuracy = 0.1)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "Infection fatality ratio", x = "") +
    scale_color_manual(values = cols, breaks = c("central", which)) +
    scale_fill_manual(values = cols, breaks = c("central", which)) +
    facet_grid(cols = vars(scenario),
               labeller = labeller(scenario = scenario_labels)) +
    theme(strip.text = element_blank())
  
  rt_p <- rt %>%
    filter(scenario %in% c("central", which)) %>%
    ggplot(., aes(date)) +
    geom_line(aes(y = mean, col = scenario, linetype = state)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = scenario), alpha = 0.2) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b-%y") +
    labs(y = "Reproduction numer (Rt)", x = "") +
    scale_color_manual(values = cols, breaks = c("central", which)) +
    scale_fill_manual(values = cols, breaks = c("central", which)) +
    facet_grid(cols = vars(scenario),
               labeller = labeller(scenario = scenario_labels)) +
    theme(strip.text = element_blank(),
          axis.text.x = element_blank())
  
  ((traj_plot[[1]] / traj_plot[[2]]) / rt_p / ifr_p) +
    plot_layout(heights = c(0.5, 0.5, 1, 1))
}

