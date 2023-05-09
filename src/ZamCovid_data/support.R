plot_data <- function(df) {
  
  plot <- NULL
  for (i in c("deaths", "pcr", "sero")) {
    
    tmp <- df %>% pivot_longer(!date) %>% filter(str_detect(name, i))
    
    if (i == "sero") {
      age_groups <- gsub("sero_positive_", "",
                         grep("sero_positive_", unique(tmp$name), value = TRUE))
      for (j in age_groups) {
        ret <- tmp %>% filter(str_detect(name, j)) %>% pivot_wider()
        colnames(ret) <- gsub(paste0("_", j), "", colnames(ret))
        nm <- paste0("seropos_", j)
        
        ret <- ret %>%  
          mutate(seropos = sero_positive / sero_total) %>%
          select(date, seropos) %>%
          `colnames<-`(c("date", nm)) %>%
          pivot_longer(!date)
        
        tmp <- rbind(tmp, ret)
      }
      tmp <- tmp %>%
        filter(!str_detect(name, "sero_positive_"),
               !str_detect(name, "sero_total_"))
    }
    nm <- if (i == "deaths") {
      "Daily deaths (all-cause mortuary)"
    } else if (i == "pcr") {
      "Daily PCR-positive cases (health facility testing)"
    } else if (i == "sero") {
      "Weekly seropositive (TREATS COVID)"
    }
    
    plot_all <- tmp %>% filter(str_detect(name, "_all")) %>%
      ggplot(., aes(date, value)) +
      geom_point(size = 0.7) + labs(x = "", y = "", title = nm)
    plot_age <- tmp %>% filter(!str_detect(name, "_all")) %>%
      ggplot(., aes(date, value, col = name)) +
      geom_point(size = 0.7) + labs(x = "", y = "")
    
    apply_theme <- function(gg) {
      gg +
        facet_grid(cols = vars(name), scales = "free_y") +
        scale_y_continuous(expand = c(0, NA)) +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
        theme_minimal() +
        theme(axis.line = element_line(),
              legend.position = "none",
              axis.text.x = element_text(angle = 45, size = 9, vjust = 0.7))
    }
    
    
    plot[[i]] <- apply_theme(plot_all) + apply_theme(plot_age) +
      plot_layout(widths = c(0.33, 1))
    
  }
  plot[[1]] / plot[[2]] / plot[[3]]
}


get_baseline_deaths <- function(...) {
  
  read_csv("data/kabwe_historic_deaths.csv") %>%
    mutate(month = factor(month, levels = month.abb)) %>%
    pivot_longer(!month, names_to = "year", values_to = "deaths") %>%
    mutate(year = as.integer(year),
           month_name = month,
           month = factor(match(month, month.abb)),
           site = "historic") %>%
    filter(year < 2020)
  
}
