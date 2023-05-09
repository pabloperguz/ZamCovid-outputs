current_ll <- deaths %>%
  mutate(month = factor(format(date, "%b"), levels = unique(historic_deaths$month)),
         year = factor(format(date, "%Y"))) %>%
  select(month, value = deaths_all, year) %>%
  group_by(month, year) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>%
  mutate(dataset = "linelist")

current_agg <- historic_deaths %>%
  pivot_longer(!month, names_to = "year") %>%
  filter(year %in% c("2020", "2021")) %>%
  group_by(month) %>%
  mutate(month = factor(month, levels = unique(historic_deaths$month))) %>%
  mutate(dataset = "aggregate") %>%
  select(month, year, value, dataset)

current_ll %>%
  rbind(., current_agg) %>%
  ggplot(., aes(month, value, fill = dataset, group = dataset)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.9) +
  facet_wrap(~year) +
  labs(y = "Monthly deaths", x = "") +
  scale_fill_manual(values = c("blue3", "red3")) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
  theme(axis.line = element_line())




historic <- historic_deaths %>%
  pivot_longer(!month, names_to = "year") %>%
  group_by(month) %>%
  mutate(month = factor(month, levels = unique(historic_deaths$month))) %>%
  select(month, year, value)

ggplot(historic, aes(month, value, col = year, group = year)) +
  geom_point(size = 3) +
  geom_line(linetype = 3) +
  labs(y = "Monthly deaths", x = "") +
  labs(title = "Just aggregate dataset") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
  theme(axis.line = element_line(),
        legend.position = "none") +
historic %>%
  filter(year %in% c("2017", "2018", "2019")) %>%
  rbind(., current_ll) %>%
  ggplot(., aes(month, value, col = year, group = year)) +
  geom_point(size = 3) +
  geom_line(linetype = 3) +
  labs(y = "Monthly deaths", x = "") +
  labs(title = "Aggregate (2017-19) and linelist (2020-21)") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
  theme(axis.line = element_line())



historic_poiss <- historic_deaths %>%
  pivot_longer(!month, names_to = "year") %>%
  filter(year %in% c("2017", "2018", "2019")) %>%
  group_by(month) %>%
  mutate(month = factor(month, levels = unique(historic_deaths$month))) %>%
  summarise(value = ceiling(mean(value))) %>%
  mutate(year = "Historic Poisson") %>%
  select(month, year, value) %>%
  rowwise() %>%
  mutate(lb = poisson.test(value)$conf.int[[1]],
         ub = poisson.test(value)$conf.int[[2]])

current_ll %>%
  mutate(lb = NA_real_,
         ub = NA_real_) %>%
  rbind(., historic_poiss) %>%
  ggplot(., aes(month, value, col = year, group = year)) +
  geom_point(size = 3) +
  geom_line(linetype = 3) +
  geom_errorbar(aes(ymin = lb, ymax = ub)) +
  labs(y = "Monthly deaths", x = "") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
  theme(axis.line = element_line())



## Calculate excess deaths using linear model

pacman::p_load(tidyverse, readxl, lubridate, ConvCalendar, 
               janitor, broom, glue, scales, ciTools, survey,
               update = FALSE, install = FALSE)

ggplottheme <- theme_classic() + 
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(color = "#cbcbcb"),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA),
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    axis.ticks = element_blank(),
    strip.background = element_blank()
  )
theme_set(ggplottheme)



alpha <- 0.95
deaths_for_model <- historic_deaths %>%
  mutate(month_name = month,
         month = match(month, month.abb)) %>%
  pivot_longer(!c(month, month_name), names_to = "year", values_to = "deaths") %>%
  mutate(year = as.numeric(year),
         deaths = floor(deaths * 1.937550))

model_deaths <- deaths_for_model %>%
  filter(year < 2020) %>%
  mutate(month = factor(month),
         region = "kabwe") %>%
  nest_by(region) %>%
  mutate(model = list(lm(deaths ~ year + month, data = data))) %>%
  mutate(sigma = broom::glance(model)$sigma,
         nobs = broom::glance(model)$nobs,
         df = model$df,
         t_df = abs(qt((1 - alpha) / 2, df)))

expected_deaths <- model_deaths %>%
  summarise(df = model$df,
            broom::augment(model,
                           newdata = tibble(month = c(1:12) %>% as.factor(),
                                            year = 2020),
                           interval = "prediction",
                           se_fit = TRUE)) %>%
  mutate(month_num = as.numeric(month)) %>%
  rename(expected = .fitted,
         lb = .lower,
         ub = .upper) %>%
  mutate(t_df = abs(qt((1 - alpha) / 2, df)),
         se_fit_pred = (ub - expected) / t_df)


agg_2020 <- current_agg[current_agg$year == "2020", ]
ll_2020 <- current_ll[current_ll$year == "2020", ]
expected_deaths %>%
  mutate(month = factor(month.abb, levels = month.abb)) %>%
  select(month, year, expected, lb, ub) %>%
  ggplot(., aes(month)) +
  geom_pointrange(aes(y = expected, ymin = lb, ymax = ub, col = "Expected"),
                  linetype = 3) +
  # geom_point(aes(y = agg_2020$value, col = "Data (agg)"), size = 2) +
  geom_point(aes(y = ll_2020$value, col = "Data (linelist)"), size = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  labs(y = "Monthly deaths 2020", x = "") +
  # scale_color_manual(values = c("blue3", "red3", "green4"),
  #                    breaks = c("Data (linelist)", "Data (agg)", "Expected")) +
  theme(legend.title = element_blank())

# count <- 0:4
# deaths <- c(144, 91, 32, 11, 2)
# 
# n <- sum(deaths)
# x <- sum(count * deaths)
# 
# lambda <- x/n
# 
# DescTools::PoissonCI(x=x, n=n, method = c("exact","score", "wald", "byar"))
# 
# exp <- dpois(0:4, lambda) * n
# 
# barplot(rbind(deaths, exp * n/sum(exp)), names=0:4, beside=TRUE,
#         col=c("red4", "blue4"), main = "Deaths", xlab = "count")
# legend("topright", legend=c("observed","expected"), fill=c("red4", "blue4"),
#        bg="white")
# 
# 
# ## SMR, Welsh Nickel workers
# DescTools::PoissonCI(x=137, n=24.19893)
# 
## My data
# n <- sum(historic$value)
# x <- sum(1:(nrow(historic)) * historic$value)
# lambda <- x / n
# 
# DescTools::PoissonCI(x, n, method = c("exact","score", "wald", "byar"))
# 
# exp <- dpois(1:(nrow(historic)), lambda) * n
# 
# barplot(rbind(historic$value, exp * n / sum(exp)), names=0:11, beside=TRUE,
#         col=c("red4", "blue4"), main = "Deaths Kabwe", xlab = "count")
# legend("topright", legend=c("observed","expected"), fill=c("red4", "blue4"),
#        bg="white")
# 
# 
# t <- poisson.test(historic$value[1])
# 
# t$conf.int[[1]]
# t$conf.int[[2]]
# 




