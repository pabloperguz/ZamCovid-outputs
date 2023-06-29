plot_data <- function(df) {
  
  df <- df %>%
    mutate(district = stringr::str_to_title(district))
  
  
  # Serology
  sero <- df %>%
    select(date, district, sero_pos_all, sero_tot_all) %>%
    filter(!is.na(sero_tot_all))
  
  sero <- cbind(sero, Hmisc::binconf(sero$sero_pos_all, sero$sero_tot_all)) %>%
    rename(mean = PointEst, lb = Lower, ub = Upper)
  
  sero <- ggplot(sero, aes(date)) +
    geom_pointrange(aes(y = mean, ymin = lb, ymax = ub)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.2),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(limits = c(as.Date("2020-06-15"), as.Date("2020-08-01"))) +
    labs(y = "% seropositive", x = "") +
    facet_wrap(~district, ncol = 1) +
    theme_minimal() +
    theme(axis.line = element_line())
  
  # PCR-positivity
  pcr <- df %>%
    select(date, district, pcr_pos_all, pcr_tot_all) %>%
    filter(!is.na(pcr_tot_all))
  
  pcr <- cbind(pcr, Hmisc::binconf(pcr$pcr_pos_all, pcr$pcr_tot_all)) %>%
    rename(mean = PointEst, lb = Lower, ub = Upper)
  
  pcr <- ggplot(pcr, aes(date)) +
    geom_pointrange(aes(y = mean, ymin = lb, ymax = ub)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.2),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(limits = c(as.Date("2020-06-15"), as.Date("2020-08-01"))) +
    labs(y = "PCR positivity", x = "") +
    facet_wrap(~district, ncol = 1) +
    theme_minimal() +
    theme(axis.line = element_line())
    
  # Deaths
  deaths <- df %>%
    select(date, district, deaths = deaths_all) %>%
    mutate(date = lubridate::floor_date(date, "month"),
           source = "pandemic") %>%
    group_by(date, district, source) %>%
    summarise(deaths = sum(deaths, na.rm = TRUE))
  
  historic <- historic_deaths %>%
    pivot_longer(!c(district, month)) %>%
    mutate(date = as.Date(paste0(name, "-", month, "-01"), format = "%Y-%b-%d"),
           district = stringr::str_to_sentence(district),
           source = "historic") %>%
    select(date, district, source, deaths = value)
  
  deaths <- rbind(deaths, historic)
  
  deaths <- ggplot(deaths, aes(date, deaths, col = source)) +
    geom_point() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                       limits = c(0, NA)) +
    labs(y = "Monthly deaths", x = "") +
    facet_wrap(~district, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(axis.line = element_line(),
          legend.title = element_blank())
  
  sero | pcr | deaths
}
